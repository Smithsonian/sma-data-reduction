from pyuvdata.telescopes import get_telescope
from pyuvdata.uvcal import UVCal
import pyuvdata.utils as uvutils
import h5py
import numpy as np
from astropy.time import Time
import os
import shutil

def export_compass_gains(filename, outpath='.'):
    cal_obj = UVCal()
    cal_obj._set_gain()
    cal_obj._set_wide_band()

    cal_obj.ref_antenna_name = "unknown"
    cal_obj.gain_convention = "divide"
    cal_obj.telescope = get_telescope("SMA")
    cal_obj.telescope.Nants = 8
    cal_obj.telescope.x_orientation = "east"

    cal_obj.history = "Created by COMPASS v"
    cal_obj.gain_convention = "divide"
    cal_obj.sky_catalog = "COMPASS"
    cal_obj.gain_scale = "Jy"
    cal_obj.observer = "SMA"
    cal_obj._set_sky()
    cal_obj.Nants_data = 8
    cal_obj.Nfreqs = 1
    cal_obj.antenna_names = ["Ant%i" % idx for idx in range(1, 9)]
    cal_obj.antenna_numbers = np.arange(1, 9)
    cal_obj.ant_array = np.arange(1, cal_obj.Nants_data + 1)
    lon = cal_obj.telescope.location.lon.rad

    with h5py.File(filename, "r") as compass_gains:
        obsid = (''.join(np.array(compass_gains['obsID'][0].view('U1'))))
        gain_corr = compass_gains['gainApEff'][...]

        cal_obj.gain_array = (gain_corr**0.5) * (
            compass_gains['gainWinReal'][...] + (1j*compass_gains['gainWinImag'][...])
        )
        cal_obj.flag_array = np.isnan(abs(cal_obj.gain_array))
        cal_obj.quality_array = np.where(cal_obj.flag_array, 0.0, 1.0)
        cal_obj.Ntimes = len(cal_obj.gain_array)
        cal_obj.Nspws = cal_obj.gain_array.shape[2]
        cal_obj.freq_range = np.zeros((cal_obj.Nspws, 2))
        cal_obj.freq_range[:, 0] = np.min(compass_gains['freqArr'], axis=1)
        cal_obj.freq_range[:, 1] = np.max(compass_gains['freqArr'], axis=1)
        compass_version = (''.join(np.array(compass_gains['compassVersion'][0].view('U1'))))
        sou_id_arr = compass_gains['gainWinSouID'][0].astype(int)
        sou_list = ((''.join(np.array(compass_gains['gainUniqueSou'][0].view('U1')))).split(';'))
        sou_ra = compass_gains['gainSouRA'][0]
        sou_dec = compass_gains['gainSouDec'][0]
        sou_epoch = compass_gains['gainSouEpoch'][0]
        mjd_arr = compass_gains['gainWinTimes'][...]
        cal_obj.time_array = Time(np.mean(mjd_arr, axis=0), scale="tt", format="mjd").utc.jd
        cal_obj.integration_time = np.diff(mjd_arr, axis=0)[0] * 86400
        cal_obj.set_lsts_from_time_array()
        cal_obj.antenna_positions = uvutils.ECEF_from_rotECEF(
            compass_gains['antPos'][...], lon
        )
        same_freq = bool(compass_gains['sameFreq'][0,0])
        ref_ant = compass_gains['gainWinRefAnt'][...]
        sb_arr = compass_gains['sbArr'][...]
        rx1_arr = compass_gains['rx1Arr'][...]
        rx2_arr = compass_gains['rx2Arr'][...]
        sb_arr = compass_gains['sbArr'][...]
        win_arr = compass_gains['winArr'][...]

        use_gains = np.array(compass_gains['gainWinUse'][0], dtype=bool)
        use_pha = np.array(compass_gains['gainWinPha'][0], dtype=bool)
        use_amp = np.array(compass_gains['gainWinAmp'][0], dtype=bool)

    # If we have the same freq, order the windows in a good way
    if same_freq:
        gain_shape = (cal_obj.Ntimes, cal_obj.Nants_data, 2, cal_obj.Nspws // 2)
        cal_obj.gain_array = np.swapaxes(cal_obj.gain_array.reshape(gain_shape), 2, 3)
        cal_obj.flag_array = np.swapaxes(cal_obj.flag_array.reshape(gain_shape), 2, 3)
        cal_obj.quality_array = np.swapaxes(cal_obj.quality_array.reshape(gain_shape), 2, 3)
        cal_obj.Nspws = cal_obj.Nspws // 2
        cal_obj.Njones = 2
        cal_obj.jones_array = [-5, -6]  # Fix for pol obs
        cal_obj.spw_array = np.array(
            win_arr.reshape(2,-1)[0] * ((-1)**(1 - sb_arr.reshape(2,-1)[0])), dtype=int
        )
        cal_obj.freq_range = cal_obj.freq_range.reshape(2, -1, 2)[0]
    else:
        gain_shape = (cal_obj.Ntimes, cal_obj.Nants_data, 1, cal_obj.Nspws)
        cal_obj.gain_array = np.swapaxes(np.repeat(cal_obj.gain_array.reshape(gain_shape), 2, axis=2), 2, 3)
        cal_obj.flag_array = np.swapaxes(np.repeat(cal_obj.flag_array.reshape(gain_shape), 2, axis=2), 2, 3)
        cal_obj.quality_array = np.swapaxes(np.repeat(cal_obj.quality_array.reshape(gain_shape), 2, axis=2), 2, 3)
        cal_obj.Njones = 2 # Default to split_pol and fix later
        cal_obj.jones_array = [-5, -6]  # Mark unknown for now
        cal_obj.spw_array = np.array(
            win_arr * ((-1)**(1 - sb_arr)) + (512 * rx1_arr), dtype=int
        )[0]


    # Reorganize the data to be the right shape/dimensions
    cal_obj.gain_array = np.conj(np.transpose(cal_obj.gain_array, [1, 2, 0, 3]))
    cal_obj.flag_array = np.transpose(cal_obj.flag_array, [1, 2, 0, 3])
    cal_obj.quality_array = np.transpose(cal_obj.quality_array, [1, 2, 0, 3])
    cal_obj.history += compass_version

    # Source handling
    cal_obj.phase_center_id_array = sou_id_arr.astype(int)

    for idx in range(len(sou_list)):
        cal_obj._add_phase_center(
            cat_name=sou_list[idx],
            cat_type="sidereal",
            cat_lon=sou_ra[idx],
            cat_lat=sou_dec[idx],
            cat_frame="icrs",
            cat_epoch=sou_epoch[idx],
            cat_id=(idx + 1),
            info_source="file",
        )

    # Make sure that things work correctly.
    cal_obj.check()

    # First do the amplitude solns
    gain_solns_obj = cal_obj.copy()
    gain_solns_obj._select_by_index(
        time_inds=np.where(use_amp & use_gains)[0],
        spw_inds=None,
        ant_inds=None,
        freq_inds=None,
        jones_inds=None,
        history_update_string=' Selected amplitude solns.',
    )
    gain_solns_obj.gain_array = np.abs(gain_solns_obj.gain_array)
    gain_solns_obj.write_ms_cal(
        os.path.join(outpath, f'{obsid}_amp_solns.ms'), clobber=True
    )

    # Next do the phase solns
    gain_solns_obj = cal_obj.copy()
    gain_solns_obj._select_by_index(
        time_inds=np.where(use_pha & use_gains)[0],
        spw_inds=None,
        ant_inds=None,
        freq_inds=None,
        jones_inds=None,
        history_update_string=' Selected phase solns.',
    )
    gain_solns_obj.gain_array /= np.abs(gain_solns_obj.gain_array)
    gain_solns_obj.write_ms_cal(
        os.path.join(outpath, f'{obsid}_pha_solns.ms'), clobber=True
    )

    # Finally, do the self-cal sources
    for cat_id in np.unique(sou_id_arr[~use_gains]):
        selfcal_gains = cal_obj.select(phase_center_ids=cat_id, inplace=False)
        # For now, just treat this as phase-only
        # selfcal_gains.gain_array /= np.abs(selfcal_gains.gain_array)
        sou_name = selfcal_gains.phase_center_catalog[cat_id]['cat_name']
        selfcal_gains.write_ms_cal(
            os.path.join(outpath, f'{obsid}_{sou_name}_selfcal_solns.ms'), clobber=True
        )
    
    # Finally, copy the matlab file over to the main directory for the spectral solns
    _ = shutil.copy(filename, outpath)

    print("Gains conversion to MS format complete!")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Convert MIR to MS.')

    parser.add_argument('-f', '--filename',
                        type=str,
                        help='MIR filename (comma separated when reading multiple files)')

    parser.add_argument('-o', '--outpath',
                        type=str,
                        help='Output directory for MS')

    args = parser.parse_args()

    export_compass_gains(filename=args.filename, outpath=args.outpath)
