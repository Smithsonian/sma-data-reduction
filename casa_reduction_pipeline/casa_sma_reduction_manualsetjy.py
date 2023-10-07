
# Based on Todd Hunter's CASA calibration script, modified by Charlie Qi 05/2021
# The measurement set is generated with CASA 5.7.2-4.

import os
import numpy as np
import sys
from datetime import datetime

# CASA 6
from casatasks import (listobs, clearcal, flagmanager, flagdata,
                       setjy, bandpass, gaincal, applycal, blcal, fluxscale)
from casaplotms import plotms
from casatools import table, msmetadata


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

if not os.path.exists(myvis):
    raise FileNotFoundError(f"Unable to find MS: {myvis}")

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
flux_stokesI = float(sma_pipe_config.get('SMA-Pipe', 'flux_stokesI'))


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
flagdata(vis=myvis, mode='manual', spw=edgechan, flagbackup=False)
flagmanager(vis=myvis, mode='save', versionname='flag_edges')


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

# Get the avg freq.
msmd.open(myvis)
scan_idx = msmd.scannumbers()[0]
all_spws = msmd.spwsforscan(scan_idx)
mean_freq = np.mean([np.mean(msmd.chanfreqs(spw_idx)) for spw_idx in all_spws])
msmd.close()

setjy(vis=myvis, field=flux, spw='',
      fluxdensity=[flux_stokesI, 0, 0, 0],
      reffreq=f"{mean_freq/1e9}GHz",
      scalebychan=True,
      standard='manual',
      usescratch=False)


###############
# BANDPASS
###############

# phase-only selfcal
phaseshortgaincal_table = '{0}.bpself.gcal'.format(myvis)
if os.path.exists(phaseshortgaincal_table):
    os.system(f'rm -rf {phaseshortgaincal_table}')

gaincal(vis=myvis,caltable=phaseshortgaincal_table,
        field=bpcal,spw=bpchans,refant=refant, scan=bpscan,
        calmode='p',solint=min_solint,minsnr=2.0, minblperant=3,
        gaintable=priorcals)

plotms(vis=phaseshortgaincal_table,
       xaxis='time',
       yaxis='phase',
       coloraxis='spw',
       iteraxis='antenna',
       ydatacolumn='data',
       plotfile=f"{caltab_plot_path}/{phaseshortgaincal_table}.png",
       gridrows=3, gridcols=3,
       yselfscale=True,
       xconnector='line', timeconnector=True,
       showgui=False, overwrite=True, dpi=400)

# Solnorm True bandpass
bandpass_table = '{0}.bandpass.solnorm_true.bcal'.format(myvis)
if os.path.exists(bandpass_table):
    os.system(f'rm -rf {bandpass_table}')

# smooth some channels CHECK SMOOTHING WINDOW
bandpass(vis=myvis,caltable=bandpass_table,
         bandtype='B', scan=bpscan,
         field=bpcal, spw=spwrange,
         combine='scan,field',
         refant=refant,
         solint=bpsolint, solnorm=True, minblperant=3,
         fillgaps=10,  # If some channels are flagged above, interpolate over in the BP
         gaintable=[phaseshortgaincal_table] + priorcals)

# Plot bandpass phase and amplitude solutions.
plotms(vis=bandpass_table,xaxis='freq',
       yaxis='phase',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       plotfile=f"{caltab_plot_path}/{bandpass_table}.phase.png",
       gridrows=3, gridcols=3,
       yselfscale=True,
       showgui=False, overwrite=True, dpi=400)
plotms(vis=bandpass_table,xaxis='freq',
       yaxis='amp',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       plotfile=f"{caltab_plot_path}/{bandpass_table}.amp.png",
       gridrows=3, gridcols=3,
       yselfscale=True,
       showgui=False, overwrite=True, dpi=400)


# Make ap selfcal table applying the bandpass
ampphaseshortgaincal_table = '{0}.bpself.ap.gcal'.format(myvis)
if os.path.exists(ampphaseshortgaincal_table):
    os.system(f'rm -rf {ampphaseshortgaincal_table}')

gaincal(vis=myvis, caltable=ampphaseshortgaincal_table,
        field=bpcal, spw=bpchans, refant=refant, scan=bpscan,
        calmode='ap', solint=min_solint, minsnr=2.0, minblperant=3,
        gaintable=[bandpass_table] + priorcals)

plotms(vis=ampphaseshortgaincal_table,
       xaxis='time',
       yaxis='phase',
       coloraxis='spw',
       iteraxis='antenna',
       ydatacolumn='data',
       xconnector='line', timeconnector=True,
       plotfile=f"{caltab_plot_path}/{ampphaseshortgaincal_table}.png",
       gridrows=3, gridcols=3,
       yselfscale=True,
       showgui=False, overwrite=True, dpi=400)
plotms(vis=ampphaseshortgaincal_table,
       xaxis='time',
       yaxis='amp',
       coloraxis='spw',
       iteraxis='antenna',
       ydatacolumn='data',
       xconnector='line', timeconnector=True,
       plotfile=f"{caltab_plot_path}/{ampphaseshortgaincal_table}.png",
       gridrows=3, gridcols=3,
       yselfscale=True,
       showgui=False, overwrite=True, dpi=400)



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
        gaintable=[bandpass_table] + priorcals)
# USB
gaincal(vis=myvis,caltable=gain_phase_int_table, append=True,
        field=calfields, refant=refant,
        combine='spw',spw='6~11',
        calmode='p',solint=min_solint,minsnr=2.0,minblperant=3,
        gaintable=[bandpass_table] + priorcals)


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
        gaintable=[bandpass_table] + priorcals)
# USB
gaincal(vis=myvis,caltable=gain_phase_scan_table,append=True,
        field=calfields,refant=refant,
        combine='spw',spw='6~11',
        calmode='p',solint='300s',minsnr=2.0,minblperant=3,
        gaintable=[bandpass_table] + priorcals)

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
        spwmap=[[],[0,0,0,0,0,0,6,6,6,6,6,6]],
        calmode='a', solint='10min', minsnr=3.0, minblperant=3,
        gaintable=[bandpass_table,
                   gain_phase_int_table] + priorcals)


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
from utils.fluxscale_fit_plots import fluxscale_to_tables, plot_flux_fits

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
        spwmap=[[],[0,0,0,0,0,0,6,6,6,6,6,6], []],
        calmode='a', solint='300s', gaintype='G',
        minsnr=2.0, minblperant=3,
        solnorm=False,
        gaintable=[bandpass_table,
                   gain_phase_int_table] + priorcals)


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
         gaintable=[bandpass_table,
                    gain_phase_int_table,
                    gain_amp_scan_final_table] + priorcals,
         interp=['linear,linearflag',
                 'linear,linear',
                 'nearest,linear'],
         spwmap=[[], [0,0,0,0,0,0,6,6,6,6,6,6], []],
         gainfield=[bpcal, bpcal, bpcal],
         flagbackup=False,
         calwt=False)



for pcal in bothpcal.split(","):

    applycal(vis=myvis,field=pcal,
            spw=spwrange,
            gaintable=[bandpass_table,
                       gain_phase_int_table,
                       gain_amp_scan_final_table] + priorcals,
            interp=['linear,linearflag',
                    'linearPD,linear',
                    'nearest,linear'],
            spwmap=[[], [0,0,0,0,0,0,6,6,6,6,6,6], []],
            gainfield=[bpcal, pcal, pcal],
            flagbackup=False, calwt=False)


# If the flux overlaps with the other calibrations, no need to
# reapply the calibration.
if flux != bpcal and flux not in bothpcal:
    # Flux calibration:
    applycal(vis=myvis,field=flux,
            spw=spwrange,
            gaintable=[bandpass_table,
                        gain_phase_int_table,
                        gain_amp_scan_final_table] + priorcals,
            interp=['nearest','nearestPD','nearest'],
            spwmap=[[], [0,0,0,0,0,0,6,6,6,6,6,6], []],
            gainfield=[bpcal, flux, flux],
            flagbackup=False, calwt=False)



# Science calibration:

applycal(vis=myvis,field=science_fields,
        spw=spwrange,
        gaintable=[bandpass_table,
                   gain_phase_scan_table,
                   gain_amp_scan_final_table] + priorcals,
        interp=['linear,linearflag',
                'linearPD,linear',
                'linear,linear'],
        spwmap=[[], [0,0,0,0,0,0,6,6,6,6,6,6], []],
        gainfield=[bpcal, bothpcal, bothpcal],
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

