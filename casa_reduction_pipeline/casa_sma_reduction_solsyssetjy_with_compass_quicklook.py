
# Based on Todd Hunter's CASA calibration script, modified by Charlie Qi 05/2021
# The measurement set is generated with CASA 5.7.2-4.

import os
import numpy as np
import sys
from datetime import datetime

# Quicklook SMA export functions:
from quicklook_sma.utilities import read_config
# Additional QA plotting routines
from quicklook_sma import make_qa_tables, make_all_caltable_txt
# Quicklook imaging
from quicklook_sma.quicklook_imaging import quicklook_continuum_imaging


# CASA 6
from casatasks import (listobs, clearcal, flagmanager, flagdata,
                       setjy, bandpass, gaincal, applycal, blcal, fluxscale)
from casaplotms import plotms
from casatools import table, msmetadata


from casatools import logsink

casalog = logsink()


this_date = datetime.now()

config_filename = sys.argv[-1]

if not os.path.exists(config_filename):
    # raise ValueError(f"Cannot fine config filename {config_filename}")
    raise ValueError("Cannot fine config filename {0}".format(config_filename))

import configparser

config_file = configparser.ConfigParser()
config_file.read(config_filename)

sma_pipe_config = config_file

###################
myvis = sma_pipe_config.get('SMA-Pipe', 'myvis')
###################

manual_flag_file = sma_pipe_config.get('SMA-Pipe', 'manual_flag_file')
# Set to None if not filename is given
if len(manual_flag_file) == 0:
    manual_flag_file = None
else:
    if not os.path.exists(manual_flag_file):
        raise ValueError(f"Cannot find the specified flagging file {manual_flag_file}")

# Reset the corrected column when re-running the script.
restart_pipeline = sma_pipe_config.getboolean('SMA-Pipe', 'restart_pipeline')

# Interactive run
interactive_on = sma_pipe_config.getboolean('SMA-Pipe', 'interactive_on')


flux = sma_pipe_config.get('SMA-Pipe', 'flux')
# This script will use the SS models and so doesn't need the manually
# set brightness.
# flux_stokesI = float(sma_pipe_config.get('SMA-Pipe', 'flux_stokesI'))


bpcal = sma_pipe_config.get('SMA-Pipe', 'bpcal')
pcal1 = sma_pipe_config.get('SMA-Pipe', 'pcal1')
pcal2 = sma_pipe_config.get('SMA-Pipe', 'pcal2')
science_fields = sma_pipe_config.get('SMA-Pipe', 'science_fields')
is_mosaic = sma_pipe_config.getboolean('SMA-Pipe', 'is_mosaic')

if len(pcal2) > 0:
    bothpcal = ",".join([pcal1, pcal2])
else:
    bothpcal= pcal1

if flux == bpcal or flux in bothpcal:
    calfields= ",".join([bpcal, bothpcal])
else:
    calfields= ",".join([bpcal, bothpcal, flux])


if is_mosaic:
    # Find all matching target names
    tb.open("{0}/FIELD".format(myvis))
    field_names = tb.getcol('NAME')
    tb.close()

    science_match = science_fields.strip("*")
    science_field_list = []

    for field in field_names:
        if science_match in field:
            science_field_list.append(field)

    science_fields = ",".join(science_field_list)

bpchans = sma_pipe_config.get('SMA-Pipe', 'bpchans')
calchans = sma_pipe_config.get('SMA-Pipe', 'calchans')
spwrange = sma_pipe_config.get('SMA-Pipe', 'spwrange')
nspws = int(sma_pipe_config.get('SMA-Pipe', 'nspws'))
bpscan = sma_pipe_config.get('SMA-Pipe', 'bpscan')

rechunk = int(sma_pipe_config.get('SMA-Pipe', 'rechunk'))
edgechan_frac = float(sma_pipe_config.get('SMA-Pipe', 'edgechan_frac'))

chan_num = int(16384 / rechunk)
edge_chan_low = int(np.floor(chan_num * edgechan_frac))
edge_chan_high = int(np.floor(chan_num * (1. - edgechan_frac)))

# Flag edges of the SPWs
edgechan = "{0}:0~{1};{2}~{3}".format(spwrange, edge_chan_low,
                                      edge_chan_high, chan_num-1)

try:
    extra_chunk56_edgeflag = sma_pipe_config.getboolean('SMA-Pipe',
                                                        'extra_chunk56_edgeflag')
except Exception:
    extra_chunk56_edgeflag = False

# Add a bit of extra flagging between chunks 5/6 where there's a steeper drop in sensitivity
if extra_chunk56_edgeflag:
    # add an extra ~10% flagging between those two only.
    edge_frac_extra = 0.12
    edge_chan_high_56 = int(np.floor(chan_num * (1. - edge_frac_extra)))

    # This accounts for the reversal in channel vs frequency ordering
    chunk56_lsb = f"0~1:{edge_chan_high_56}~{chan_num-1}"
    chunk56_usb = f"10~11:{edge_chan_high_56}~{chan_num-1}"

    # Add  these to the edgechan selection
    edgechan += f", {chunk56_lsb}, {chunk56_usb}"

refant = sma_pipe_config.get('SMA-Pipe', 'refant')
bpsolint= sma_pipe_config.get('SMA-Pipe', 'bpsolint')

# Minimum time for solution intervals.
# Defaults to 30s (the normal SMA integration time)
try:
    min_solint = sma_pipe_config.get('SMA-Pipe', 'min_solint')
except Exception:
    min_solint = '30s'


###############
# Setup output directories
###############
caltab_plot_path = sma_pipe_config.get('SMA-Pipe', 'caltab_plot_path')

if not os.path.exists(caltab_plot_path):
    os.mkdir(caltab_plot_path)

plot_path = sma_pipe_config.get('SMA-Pipe', 'plot_path')

if not os.path.exists(plot_path):
    os.mkdir(plot_path)

scitarget_plot_path = sma_pipe_config.get('SMA-Pipe', 'scitarget_plot_path')

if not os.path.exists(scitarget_plot_path):
    os.mkdir(scitarget_plot_path)


###############
# Restart and clear previous calibration/model columns
# Restore original flags.
###############

if restart_pipeline:
    # Overwrite the CORRECTED column.
    clearcal(vis=myvis)
    try:
        flagmanager(vis=myvis, mode='restore', versionname='original')
    except RuntimeError:
        casalog.post("Unable to find original flag version")
        print("Unable to find original flag version")
else:
    # Backup a version of the original flags
    flagversions = flagmanager(vis=myvis, mode='list')
    flagversion_names = [flagversions[key]['name']
                         for key in flagversions.keys()
                         if key != "MS"]
    if "original" not in flagversion_names:
        flagmanager(vis=myvis, mode='save', versionname='original')

###############
# List MS file contents
###############

listobs(myvis)

###############
# Edge and birdie flagging
###############

# Ensure zeros are flagged
flagdata(vis=myvis, mode='clip', clipzeros=True, flagbackup=False)


# Flag edges and birdies
# flagdata(vis=myvis, mode='manual', spw=edgechan, flagbackup=False)
# flagmanager(vis=myvis, mode='save', versionname='flag_edges')


# Use plotms to view the bright quasars interactively.
# We are looking for large outliers in either channels or time.

# Read in flag commands from a flagging txt file
if manual_flag_file is not None:
    flagdata(vis=myvis, mode='list',
            inpfile=manual_flag_file,
            flagbackup=False)

    flagmanager(vis=myvis, mode='save', versionname='manual_flagging')

else:
    print("No manual flagging file given.")

# From here, apply flags. Then replot to check that we caught all the spikes.

if interactive_on:
    # First look for outliers in amp vs channel
    for myfield in list(set(calfields.split(","))):
        plotms(vis=myvis, xaxis='channel',
            yaxis='amp',field=myfield, avgtime='1e8', avgscan=False,
            coloraxis='ant1',iteraxis='spw', ydatacolumn='data',
            gridrows=4, gridcols=3, yselfscale=True, showgui=True)
        # input(f"Done adding freq. flags for {myfield}?")
        input("Done adding freq. flags for {0}?".format(myfield))

    # Next look for time outliers
    for myfield in list(set(calfields.split(","))):
        plotms(vis=myvis, xaxis='time',
            yaxis='amp',field=myfield, avgchannel='1e8',
            coloraxis='ant1',iteraxis='spw', ydatacolumn='data',
            gridrows=4, gridcols=3, yselfscale=True, showgui=True)
        # input(f"Done adding time flags for {myfield}?")
        input("Done adding time flags for {0}?".format(myfield))

    # input("Stop here to run new flagging commands.")
    input("Stop here to run new flagging commands.")

    print("Exiting after interactive examination. Re-run with interactive_on=False "
          "and set the manual_flagging_file in the config file to run the pipeline.")

    sys.exit(0)


###############
# PRIORCALS
###############

# Check for an antenna pos correction table:
antpos_table = f"{myvis[:-3]}.antpos"
if os.path.exists(antpos_table):
    priorcals = [antpos_table]
else:
    priorcals = []


###############
# SETJY
###############

setjy(vis=myvis, field=flux, spw='',
      scalebychan=True,
      standard='Butler-JPL-Horizons 2012',
      usescratch=False)



###############
# GAINS
###############

# per-int phase with combining spws across sidebands.
gain_phase_int_table = '{0}.intphase_combinespw.gcal'.format(myvis)
if os.path.exists(gain_phase_int_table):
    os.system(f'rm -rf {gain_phase_int_table}')

# Solve per sideband
# LSB
gaincal(vis=myvis,caltable=gain_phase_int_table,
        field=calfields,refant=refant,
        combine='spw',spw='0~5',
        calmode='p',solint=min_solint,minsnr=2.0,minblperant=3,
        gaintable=priorcals)
# USB
gaincal(vis=myvis,caltable=gain_phase_int_table, append=True,
        field=calfields, refant=refant,
        combine='spw',spw='6~11',
        calmode='p',solint=min_solint,minsnr=2.0,minblperant=3,
        gaintable=priorcals)


plotms(vis=gain_phase_int_table,xaxis='time',
       yaxis='phase',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       xconnector='line', timeconnector=True,
       plotfile=f"{caltab_plot_path}/{gain_phase_int_table}.png",
       gridrows=3, gridcols=3,
       yselfscale=True, showgui=False, overwrite=True, dpi=400)

# per-scan phase with combining spws
gain_phase_scan_table = '{0}.scanphase_combinespw.gcal'.format(myvis)
if os.path.exists(gain_phase_scan_table):
    os.system(f'rm -rf {gain_phase_scan_table}')

# Solve per sideband
# LSB
gaincal(vis=myvis,caltable=gain_phase_scan_table,
        field=calfields,refant=refant,
        combine='spw',spw='0~5',
        calmode='p',solint='300s',minsnr=2.0,minblperant=3,
        gaintable=priorcals)
# USB
gaincal(vis=myvis,caltable=gain_phase_scan_table,append=True,
        field=calfields,refant=refant,
        combine='spw',spw='6~11',
        calmode='p',solint='300s',minsnr=2.0,minblperant=3,
        gaintable=priorcals)

plotms(vis=gain_phase_scan_table,xaxis='time',
       yaxis='phase',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       xconnector='line', timeconnector=True,
       plotfile=f"{caltab_plot_path}/{gain_phase_scan_table}.png",
       gridrows=3, gridcols=3,
       yselfscale=True, showgui=False, overwrite=True, dpi=400)

# Amplitude gain calibration

gain_amp_scan_table = '{0}.amp.gcal'.format(myvis)
if os.path.exists(gain_amp_scan_table):
    os.system(f'rm -rf {gain_amp_scan_table}')

gaincal(vis=myvis,caltable=gain_amp_scan_table,
        field=calfields,spw=calchans,refant=refant,
        combine='',
        spwmap=[[0,0,0,0,0,0,6,6,6,6,6,6]],
        calmode='a', solint='10min', minsnr=3.0, minblperant=3,
        gaintable=[gain_phase_int_table] + priorcals)


plotms(vis=gain_amp_scan_table,xaxis='time',
       yaxis='amp',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       xconnector='line', timeconnector=True,
       plotfile=f"{caltab_plot_path}/{gain_amp_scan_table}_time.png",
       gridrows=3, gridcols=3,
       yselfscale=True, showgui=False, overwrite=True, dpi=400)

plotms(vis=gain_amp_scan_table,
       xaxis='freq', yaxis='amp',
       coloraxis='spw',iteraxis='antenna', ydatacolumn='data',
       xconnector='line', timeconnector=True,
       plotfile=f"{caltab_plot_path}/{gain_amp_scan_table}_freq.png",
       gridrows=3, gridcols=3,
       yselfscale=True, showgui=False, overwrite=True, dpi=400)


# Bootstrap fluxes from our flux cal (see setjy above)
# to the gain calibrators
fluxboot_table = '{0}.flux.cal'.format(myvis)
if os.path.exists(fluxboot_table):
    os.system(f'rm -rf {fluxboot_table}')


if flux == bpcal:
    transfer_fields = [bothpcal]
elif flux == pcal1:
    transfer_fields = [bpcal, pcal2]
elif flux == pcal2:
    transfer_fields = [bpcal, pcal1]
else:
    transfer_fields = [bpcal, bothpcal]


fluxresults = fluxscale(vis=myvis,
                        caltable=gain_amp_scan_table,
                        refspwmap=[-1],
                        transfer=transfer_fields,
                        fluxtable=fluxboot_table,
                        reference=flux,
                        fitorder=1,
                        incremental=True)

fluxdict_numpyfile = "{0}.fluxresults.npy".format(myvis)
if os.path.exists(fluxdict_numpyfile):
    os.remove(fluxdict_numpyfile)

np.save(fluxdict_numpyfile,
        fluxresults, allow_pickle=True)


# Import functions for summarizing/parsing the fluxscale output.
# from utils.fluxscale_fit_plots import fluxscale_to_tables, plot_flux_fits

# Inspect the fit solution and visualize.
def fluxscale_to_tables(fluxresults):
    '''
    Make a fit solution table and flux per SPW table with measured +
    predicted fluxes.
    '''

    from astropy.table import Table, Column

    freqs = Column(fluxresults['freq'], name='freq')
    spws = Column(fluxresults['spwID'], name='spw')

    num_spws = len(spws)

    # numfield = num_keys - (freq, spwID, spwName)
    num_fields = len(fluxresults.keys()) - 3

    field_ids = [key for key in fluxresults.keys() if key.isdigit()]

    # Make a column of table names for the data table and fit table
    field_ids_fit = Column(field_ids, name='field_idx')
    field_names_fit = Column([fluxresults[idx]['fieldName'] for idx in field_ids_fit],
                              name='field')

    all_fields_ids = []
    all_field_spws = []
    all_field_freqs = []
    for field_id in field_ids:
        all_fields_ids.extend([field_id] * num_spws)
        all_field_spws.extend(list(fluxresults['spwID']))
        all_field_freqs.extend(list(fluxresults['freq']))
    field_ids_data = Column(all_fields_ids, name='field_idx')
    field_spws_data = Column(all_field_spws, name='spw')
    field_freqs_data = Column(all_field_freqs, name='freq')
    field_names_data = Column([fluxresults[idx]['fieldName'] for idx in all_fields_ids],
                              name='field')

    for nn, field_id in enumerate(field_ids):

        # field_name = fluxresults[field_id]['fieldName']

        fit_flux = fluxresults[field_id]['fitFluxd']
        fit_fluxerr = fluxresults[field_id]['fitFluxdErr']
        fit_flux_reffreq = fluxresults[field_id]['fitRefFreq']

        # 𝑙𝑜𝑔(𝑆𝜈)=𝑎𝑜+𝑎1∗(𝑙𝑜𝑔(𝜈/𝜈0))+𝑎2∗(𝑙𝑜𝑔(𝜈/𝜈0))∗∗2
        fit_pars = fluxresults[field_id]['spidx']
        fit_errs = fluxresults[field_id]['spidxerr']

        num_params = len(fit_pars)

        fit_vals = np.array([fit_flux, fit_fluxerr, fit_flux_reffreq])
        fit_vals = np.append(fit_vals, fit_pars)
        fit_vals = np.append(fit_vals, fit_errs)

        # 4 * [IQUV], 4 * [IQUV err]
        data_vals = np.zeros((len(spws), 8))

        for ii, this_spw in enumerate(spws):
            data_vals[ii, :] = np.append(fluxresults[field_id][str(this_spw)]['fluxd'],
                                         fluxresults[field_id][str(this_spw)]['fluxdErr'])

        if nn == 0:
            all_data_vals = data_vals
            all_fit_vals = fit_vals
        else:
            all_data_vals = np.vstack([all_data_vals, data_vals])
            all_fit_vals = np.vstack([all_fit_vals, fit_vals])

    # Now construct the final tables:
    # Fit table
    fit_cols = ['fitFluxd', 'fitFluxdErr', 'fitRefFreq']
    fit_cols.extend([f'a{num}' for num in range(num_params)])
    fit_cols.extend([f'a{num}_err' for num in range(num_params)])

    fit_table = Table(data=all_fit_vals, names=fit_cols)
    fit_table.add_column(field_ids_fit)
    fit_table.add_column(field_names_fit)

    # Data table
    data_cols = ['I', 'Q', 'U', 'V', 'I_err', 'Q_err', 'U_err', 'V_err']
    data_table = Table(data=all_data_vals, names=data_cols)
    data_table.add_column(field_ids_data)
    data_table.add_column(field_names_data)
    data_table.add_column(field_spws_data)
    data_table.add_column(field_freqs_data)

    return data_table, fit_table


def log_fit_func(freq, a0, a1, a2, ref_freq):
    return a0 + a1 * np.log10(freq / ref_freq) + a2 * np.log10(freq / ref_freq)**2


def return_flux_models(fit_table):
    '''
    𝑙𝑜𝑔(𝑆𝜈)=𝑎𝑜+𝑎1∗(𝑙𝑜𝑔(𝜈/𝜈0))+𝑎2∗(𝑙𝑜𝑔(𝜈/𝜈0))∗∗2
    '''

    from functools import partial

    flux_models = {}

    for ii, field in enumerate(fit_table['field']):

        a0 = fit_table['a0'][ii]
        a1 = fit_table['a1'][ii]
        ref_freq = fit_table['fitRefFreq'][ii]
        flux0 = fit_table['fitFluxd'][ii]

        if 'a2' in fit_table.colnames:
            a2 = fit_table['a2'][ii]
        else:
            a2 = 0.0

        flux_models[field] = partial(log_fit_func, a0=a0, a1=a1,
                                     a2=a2, ref_freq=ref_freq)

    return flux_models


def return_flux_fit_residuals(fit_table, data_table, flux_models=None):

    from astropy.table import Table, Column

    residuals = {}

    if flux_models is None:
        flux_models = return_flux_models(fit_table)

    for ii, field in enumerate(fit_table['field']):

        this_model_func = flux_models[field]

        mask = data_table['field'] == field

        freqs = data_table['freq'][mask]
        spws = data_table['spw'][mask]

        # Could generalize this at some point.
        stokes_I = data_table['I'][mask]
        stokes_I_err = data_table['I_err'][mask]

        resids = stokes_I - 10**this_model_func(freqs)

        # -1.0 values indicates that SPW was not used in the fit
        resids[stokes_I == -1.] = 0.

        # Calc the N-sigma deviation from the model
        sigma_resids = resids / stokes_I_err

        residual_table = Table([spws, freqs,
                                Column(resids, name="I_resid"),
                                Column(sigma_resids, name="sigma_resid")])

        residuals[field] = residual_table

    return residuals


def plot_flux_fits(fit_table, data_table,
                   flux_models=None,
                   resids_dict=None,
                   output_folder="fluxfit_plots"):
    '''
    Makes a plot per field.
    '''

    import matplotlib.pyplot as plt

    if flux_models is None:
        flux_models = return_flux_models(fit_table)

    if resids_dict is None:
        resids_dict = return_flux_fit_residuals(fit_table, data_table,
                                                flux_models=flux_models)

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for this_field in fit_table['field']:

        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True,
                                figsize=(12, 12))

        freqs = resids_dict[this_field]['freq']
        stokes_I = data_table['I'][data_table['field'] == this_field]
        stokes_I_err = data_table['I_err'][data_table['field'] == this_field]

        valid_data = stokes_I != -1.

        stokes_I_resid = resids_dict[this_field]['I_resid']

        this_model_func = flux_models[this_field]
        stokes_I_model = 10**this_model_func(freqs)

        axs[0].errorbar(freqs[valid_data] / 1e9,
                        stokes_I[valid_data],
                        yerr=stokes_I_err[valid_data])

        a0 = fit_table['a0'][fit_table['field'] == this_field][0]
        a1 = fit_table['a1'][fit_table['field'] == this_field][0]
        ref_freq_GHz = fit_table['fitRefFreq'][fit_table['field'] == this_field][0] / 1.e9
        model_label = f"{a0:.2f} + {a1:.2f} * (freq / {ref_freq_GHz:.1f})"

        axs[0].plot(freqs / 1e9,
                    stokes_I_model,
                    label=model_label)

        axs[0].legend(frameon=True, fontsize=10)

        # Residual
        axs[1].errorbar(freqs[valid_data] / 1e9,
                        stokes_I_resid[valid_data],
                        yerr=stokes_I_err[valid_data])

        axs[1].axhline(0., color='k', linewidth=2, linestyle='--', zorder=1)

        # Set symmetric y axis for residual
        resid_max = np.abs(stokes_I_resid).max()
        axs[1].set_ylim([-1.25 * resid_max, 1.25 * resid_max])

        axs[1].set_xlabel("Freq (GHz)", fontsize=10)

        axs[0].set_ylabel("Amplitude (Jy)", fontsize=10)
        axs[1].set_ylabel("Residual (Jy)", fontsize=10)

        fig.savefig(f"{output_folder}/{this_field}_fluxscale_fit.png", dpi=300)

        plt.close()


# Make and save the table results:
flux_data_table, flux_fit_table = fluxscale_to_tables(fluxresults)

flux_data_table_name = "{0}.fluxscale_data.csv".format(myvis)
if os.path.exists(flux_data_table_name):
    os.remove(flux_data_table_name)
flux_data_table.write(flux_data_table_name, overwrite=True)

flux_fits_table_name = "{0}.fluxscale_fits.csv".format(myvis)
if os.path.exists(flux_fits_table_name):
    os.remove(flux_fits_table_name)
flux_fit_table.write(flux_fits_table_name, overwrite=True)

# Make summary plots:
plot_flux_fits(flux_fit_table, flux_data_table)

# Set the model column based on the fits:
for ii, this_field in enumerate(flux_fit_table['field']):

    this_index = [flux_fit_table['a1'][ii]]
    if 'a2' in flux_fit_table.colnames:
        this_index.append(flux_fit_table['a2'][ii])

    # NOTE: only sets the Stokes I for now.
    setjy(vis=myvis, field=this_field,
          spw='',
          fluxdensity=[10**float(flux_fit_table['a0'][ii]), 0, 0, 0],
          spix=this_index,
          reffreq=f"{float(flux_fit_table['fitRefFreq'][ii])/1.e9}GHz",
          scalebychan=True,
          standard='manual',
          usescratch=False)


# Now rederive the amp gains set by the model column,.

gain_amp_scan_final_table = '{0}.amp_final.gcal'.format(myvis)
if os.path.exists(gain_amp_scan_final_table):
    os.system(f'rm -rf {gain_amp_scan_final_table}')

gaincal(vis=myvis,caltable=gain_amp_scan_final_table,
        field=calfields,spw=calchans,refant=refant,
        combine='',
        spwmap=[[0,0,0,0,0,0,6,6,6,6,6,6], []],
        calmode='a', solint='300s', gaintype='G',
        minsnr=2.0, minblperant=3,
        solnorm=False,
        gaintable=[gain_phase_int_table] + priorcals)


plotms(vis=gain_amp_scan_final_table,xaxis='time',
       yaxis='amp',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       xconnector='line', timeconnector=True,
       plotfile=f"{caltab_plot_path}/{gain_amp_scan_final_table}.png",
       gridrows=3, gridcols=3,
       yselfscale=True, showgui=False, overwrite=True, dpi=400)

plotms(vis=gain_amp_scan_table,
       xaxis='freq', yaxis='amp',
       coloraxis='spw',iteraxis='antenna', ydatacolumn='data',
       xconnector='line', timeconnector=True,
       plotfile=f"{caltab_plot_path}/{gain_amp_scan_final_table}_freq.png",
       gridrows=3, gridcols=3,
       yselfscale=True, showgui=False, overwrite=True, dpi=400)


####################################################################
# Applycal
####################################################################

flagmanager(vis=myvis,mode='save',versionname='beforeapplycal')

## BP cal apply:

# With perSB mapping for fluxscaling
# This WORKS and should be consistent with how MIR applies
# flux calibration
applycal(vis=myvis,field=bpcal,
         spw=spwrange,
         gaintable=[gain_phase_int_table,
                    gain_amp_scan_final_table] + priorcals,
         interp=['linear,linear',
                 'nearest,linear'],
         spwmap=[[0,0,0,0,0,0,6,6,6,6,6,6], []],
         gainfield=[bpcal, bpcal],
         flagbackup=False,
         calwt=False)



for pcal in bothpcal.split(","):

    applycal(vis=myvis,field=pcal,
            spw=spwrange,
            gaintable=[gain_phase_int_table,
                       gain_amp_scan_final_table] + priorcals,
            interp=['linearPD,linear',
                    'nearest,linear'],
            spwmap=[[0,0,0,0,0,0,6,6,6,6,6,6], []],
            gainfield=[pcal, pcal],
            flagbackup=False, calwt=False)


# If the flux overlaps with the other calibrations, no need to
# reapply the calibration.
if flux != bpcal and flux not in bothpcal:
    # Flux calibration:
    applycal(vis=myvis,field=flux,
            spw=spwrange,
            gaintable=[gain_phase_int_table,
                       gain_amp_scan_final_table] + priorcals,
            interp=['nearestPD','nearest'],
            spwmap=[[0,0,0,0,0,0,6,6,6,6,6,6], []],
            gainfield=[flux, flux],
            flagbackup=False, calwt=False)



# Science calibration:

applycal(vis=myvis,field=science_fields,
        spw=spwrange,
        gaintable=[gain_phase_scan_table,
                   gain_amp_scan_final_table] + priorcals,
        interp=['linearPD,linear',
                'linear,linear'],
        spwmap=[[0,0,0,0,0,0,6,6,6,6,6,6], []],
        gainfield=[bothpcal, bothpcal],
        flagbackup=False, calwt=False,
        applymode='calflagstrict')

flagmanager(vis=myvis,mode='save',versionname='afterapplycal')


# Export the calibrated science targets
target_vis = "{0}.target".format(myvis)
if os.path.exists(target_vis):
    os.system("rm -r {}".format(target_vis))

split(vis=myvis, outputvis=target_vis,
      field=science_fields, datacolumn='CORRECTED',
      keepflags=False)

# import sys
# sys.exit(0)

# Export summary products:
timestring = this_date.strftime("%Y%m%d_%H%M")

products_folder = f"products_{timestring}"

if not os.path.exists(products_folder):

    os.mkdir(products_folder)

this_logfile = casalog.logfile()
# Copy the current log file to the products folder:
os.system(f"cp {this_logfile} {products_folder}/casa_reduction.log")

# Copy the config file used for the pipeline run to the products folder:
os.system(f"cp {config_filename} {products_folder}/")

# If given, copy the flags file used to the products folder:
if manual_flag_file is not None:
    os.system(f"cp {manual_flag_file} {products_folder}/manual_flags.txt")

# Copy the fluxscale vals, fits and plots:
os.system(f"cp {flux_data_table_name} {products_folder}/")
os.system(f"cp {flux_fits_table_name} {products_folder}/")
os.system(f"cp -r fluxfit_plots {products_folder}/")

# Copy THIS SCRIPT in to the products so it's clear what was run
# for the reduction
this_scriptname = sys.argv[-2]
os.system(f"cp {this_scriptname} {products_folder}/casa_reduction_script.py")

# Gather up calibration tables and final flag version to enable restoration:
# TODO: add this step here!


# Command line inputs.

sma_config = read_config(config_filename)

casalog.post(f"Making quicklook products for: {sma_config['myvis']}")


# --------------------------------
# Make quicklook images of targets
# --------------------------------
run_quicklook = True

# Run dirty imaging only for a quicklook
if run_quicklook:

    # Dirty images per sideband per target.
    quicklook_continuum_imaging(config_filename,
                                image_type='target',
                                niter=0, nsigma=5.,
                                output_folder="quicklook_imaging")


    # Gain and bandpass cals. No imaging of the flux cal by default.
    # It's helpful to clean for a few iterations on point source
    # calibrators.
    quicklook_continuum_imaging(config_filename,
                                image_type='calibrator',
                                niter=0, nsigma=5.,
                                output_folder="quicklook_calibrator_imaging")

    os.system("mv {0} {1}".format('quicklook_imaging', products_folder))
    os.system("mv {0} {1}".format('quicklook_calibrator_imaging', products_folder))

# ----------------------------
# Now make additional QA plots:
# -----------------------------

# Calibration table:
make_all_caltable_txt(config_filename)

# chans_to_show : int
# Number of channels to keep for visualizing in plots. Default is to average down
# to 128 per chunk/SPW. CHOOSING LARGER VALUES WILL RESULT IN LARGE DATA FILES!
chans_to_show = 128

this_config = read_config(config_filename)

# Calculate the number of channels from the given rechunk factor
chans_in_ms = 16384 / int(this_config['rechunk'])
chans_to_avg = chans_in_ms / chans_to_show
print(f"Averaging channels by {chans_to_avg} from {chans_in_ms} to {chans_to_show}")
casalog.post(f"Averaging channels by {chans_to_avg} from {chans_in_ms} to {chans_to_show}")

chans_to_avg = int(chans_to_avg)

# Per field outputs:
# Avg over all channels over time
# Avg over
make_qa_tables(config_filename,
                output_folder='scan_plots_txt',
                outtype='txt',
                overwrite=False,
                chanavg_vs_time=16384,
                chanavg_vs_chan=chans_to_avg)

# make_all_flagsummary_data(myvis, output_folder='perfield_flagfraction_txt')

# Move these folders to the products folder.
os.system("mv {0} {1}".format('final_caltable_txt', products_folder))
os.system("mv {0} {1}".format('scan_plots_txt', products_folder))
# os.system("cp -r {0} {1}".format('perfield_flagfraction_txt', products_folder))


casalog.post("Finished! To create interactive figures, run QAPlotter in the products"
             " directory.")
