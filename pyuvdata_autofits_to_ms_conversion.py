'''
Convert MIR autofits export to a CASA measurement set (MS).

>>> python pyuvdata_ms_conversion.py -f [autofits_folder] -o [output_folder]


'''

import time
import tracemalloc


def convert_mir_to_ms(datafolder,
                      outputdir='./',
                      verbose=True,
                      track_mem_usage=True):

    from pyuvdata import UVData
    from pyuvdata import get_telescope
    from numpy import median

    sma_tel = get_telescope("SMA")

    # Check if multiple filenames were given.

    uv_data = UVData()

    if track_mem_usage:
        tracemalloc.start()

    t0 = time.perf_counter()

    filenames = sorted(datafolder.rglob('*.UVF'))

    if len(filenames) == 0:
        raise ValueError('No autofits files found!')

    for ii, filename in enumerate(filenames):
        temp_uv = UVData.from_file(
            filename, file_type="uvfits", run_check_acceptability=False
        )

        # Get the object in the right shape
        temp_uv.use_future_array_shapes()
        temp_uv._set_flex_spw()

        # Change up the SPWs to make sense
        temp_uv.spw_array[:] = ii
        temp_uv.flex_spw_id_array[:] = ii

        # Fix the telescope positions, since they're recorded wrong
        temp_uv.telescope_location = sma_tel.telescope_location

        # There seems to be a minor issue w/ how UVWs are calculated in the data -- this
        # step basically calculates a scaling factor for the UVWs to correct for this.
        ant_dict = {ant: idx for idx, ant in enumerate(temp_uv.antenna_numbers)}
        ant_1_pos = temp_uv.antenna_positions[[ant_dict[idx] for idx in temp_uv.ant_1_array]]
        ant_2_pos = temp_uv.antenna_positions[[ant_dict[idx] for idx in temp_uv.ant_2_array]]
        uvw_scale_fac = median(((ant_1_pos - ant_2_pos)**2.0).sum(axis=1) /
            ((temp_uv.uvw_array**2.0).sum(axis=1))
        )**0.5

        temp_uv.uvw_array *= uvw_scale_fac

        if ii == 0:
            uv_data = temp_uv
        else:
            uv_data += temp_uv

    t3 = time.perf_counter()
    if verbose:
        print(f"Total read in time: {(t3 - t0) / 60:.4f} min")

    if track_mem_usage:
        current, peak = tracemalloc.get_traced_memory()
        print(f"MIR read: Current memory during usage is {current / 10**6:.1f}MB; Peak was {peak / 10**6:.1f}MB")
        tracemalloc.stop()

    # Now write out the file! Note that setting clobber=True will
    # overwrite any identically named measurement set.
    if verbose:
        print("Writing out MS.")

    if track_mem_usage:
        tracemalloc.start()

    # Grab the basename

    if len(filenames) == 1:
        filename_base = filenames[0].name
    else:
        # I'm assuming this is always 1 track and the date is therefore the same.
        date_str = filenames[0].name.split("_")[0]
        time_str = filenames[0].name.split("_")[1]
        for filename in filenames[1:]:
            this_time_str = filename.name.split("_")[1]
            time_str += f"_{this_time_str}"
        filename_base = f"{date_str}_{time_str}"

    out_filename = outputdir / filename_base
    uv_data.write_ms(f"{out_filename}_bin{rechunk}.ms", clobber=True,
                     run_check=run_check)
    t4 = time.perf_counter()
    if verbose:
        print(f"Done write in: {(t4 - t3) / 60:.4f} min")

    if verbose:
        print(f"Read time: {(t3 - t0) / 60:.4f} min")
        print(f"Write time: {(t4 - t3) / 60:.4f} min")
        print(f"Total time: {(t4 - t0) / 60:.4f} min")

    if track_mem_usage:
        current, peak = tracemalloc.get_traced_memory()
        print(f"MS write: Current memory during usage is {current / 10**6:.1f}MB; Peak was {peak / 10**6:.1f}MB")
        tracemalloc.stop()


if __name__ == "__main__":

    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description='Convert MIR to MS.')

    parser.add_argument('-f', '--filename',
                        type=str,
                        help='Autofits folder')

    parser.add_argument('-o', '--outpath',
                        type=str,
                        help='Output directory for MS')

    args = parser.parse_args()

    convert_mir_to_ms(Path(args.filename),
        outputdir=Path(args.outpath),
    )

