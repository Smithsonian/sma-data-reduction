'''
Reduction of SMA SWARM data using CASA 5 or 6
---------------------------------------------

Based on Todd Hunter's CASA calibration script,
modified by Chunhua Qi May-2021

Dec-2021: added flux cal (Qi)
Dec-2021: added summary plots option (Koch)
Nov-2022: notes added (Koch)
Nov-2022: add zeros flagging (Koch)

Notes:
-----

* This script will work with outputs from mir2ms (conversion in IDL) or pyuvdata
 (conversion using python). In both cases, system temperature corrections are pre-applied to the data
  so the raw data column is in Jy units. See here:
  https://lweb.cfa.harvard.edu/rtdc/SMAdata/process/casa/convertcasa/
* Using mir2ms: mir2ms will create an measurement set (MS) file _per receiver_.
  This script needs to be run on each receiver MS separately.
* Using pyuvdata: pyuvdata will combine both receivers into a single MS.
   * If both receivers are tuned together: pyuvdata correctly assigns the correlation and
     this script will calibrate both receivers appropriately (no calibration solution combines
     correlations).
   * If the receivers are tuned separately, the MS will have 2x the number of spectral windows
     (SPWs), each with a single correlation. This script requires that:
     (i) you split this MS into the SPWs for each receiver (12 SPWs for data taken in 2021-);
     (ii) OR you can modify this script to correctly assign SPWs for the gain calibration.
     For (ii), this requires changing the SPW mapping in the gain phase/amp calibration section,
     and the applycal calls. The gain calibration combines all SPWs in a sideband (contiguous 8-12
     GHz frequency coverage in 2~GHz SPWs). In this script, we assume SPWs 0~5 are one sideband and
     SPWs 6~11 are the other sideband.

Default flagging
----------------
The following are automatically flagged by the script:
  * Zero amplitudes flagged: these are bad scans that are usually caught, but a flagdata command ensures
    they are flagged here.
  * SPW edge flagging: the outer edges of each 2 GHz SPW are flagged where the passband loses sensitvity.
    The default is to flag the outer 2.5% of channels and can be altered using the `chan_frac` setting
    below. Key to this setting is the factor that the SMA data have been rechunked by, which pre-averages
    channels when fine spectral resolution is not needed in prior steps (see above on pyuvdata and here:
    https://lweb.cfa.harvard.edu/rtdc/SMAdata/process/rechunk/). This can be changed with `rechunk_factor`
    below.


Manual flagging
---------------
Additional flagging of poor data can be added to the script. Example are included (but commented out) below
(see flagdata calls with mode='manual').

The most common type of flag for SMA
observations are 'birdies' that appear at harmonics from the local oscillator. These are typically
symmetric in the upper and lower sidebands, and appear as massive amplitude outliers. These few channels
will need to be flagged by specifying the SPW and channel number for the flagdata task. The commented out
examples show this.

To find data needing manual flagging, we currently recommend visualizing the data with plotms (see example
plotms commands with the example flagging lines below). These will summarize the calibrator data in
amplitude vs. frequency and time.

Bandpass
--------
The bandpass solution is currently solved per-channel (bandtype='B'). If the SNR is not high
enough, a polynomial model (bandtype='BPOLY') make be a better choice.

We encourage users to carefully examine bandpass solutions and to try different parameters for improvement.


Flux calibration
----------------
The SMA uses solar system objects for absolute flux calibration. This uses the setjy
CASA task with model solutions from Butler-JPL-Horizons 2012. If there is no flux calibrator in the
observation, the flux levels can be set manually in the setjy command.

Gain calibration
----------------
This script produces gain calibration solutions by combining each sideband together to
maximize sensitivity. The script assumes that one sideband is contained in SPWs 0~5 and the other in
SPW 6~11. SWARM data taken before chunks 5 and 6 were included will have 4 SPWs (each with 2 GHz coverage)
and this numbering will need to be changed. Note the warning above for split tuning observations above.

Short solution time intervals
-----------------------------
This script defaults to a short time interval of a single integration, which for
the SMA is 30 seconds by default. If the integration times are shorter, you may need to moderately increase
this time interval to around 30 seconds (or based visually from examining how the phase fluctuates on few
integration timescales; the long bandpass scan may be useful to visualize this). Places where `solint='int'`
are where this short solution time interval can be changed.

Long solution time intervals
-----------------------------
The default long solution interval for the gain phase is per-scan and 10 min/per-scan for the gain amplitude
(whichever is the shortest).

Summary and reviewing calibration
---------------------------------
This script produces three sets of plots to summarize the calibration and output data:

1. Throughout the script, plotms is called to visualize the calibration solution tables. To automate,
   the `plotfile` line in the plotms calls can be uncommented out to instead save these plots as png
   files.
2. The end of the script will loop through various properties (amp/phase vs. frequency, time, uv-distance)
   of the calibrator sources and produce several plots in a `calibrator_plots` folder. This will take awhile
   to run (~30 min, depending on the data volume) but allows for assessing the quality of the calibration for
   known point-like calibrators (bandpass, gain) and extended sources (flux w/ Solar System objects).
3. Similarly, the script will produce several plots in a `target_plots` folder summarizing amplitude vs.
   frequency, time, and uv-distance for the science targets. Science targets are often fainter and so these
   plots help assess how well the calibration was applied to the targets based on their consistency
   with noise.

'''


import os

# CASA 5
from tasks import (listobs, delmod, flagdata, plotms,
                   setjy, bandpass, gaincal, applycal, blcal, fluxscale)

# CASA 6
# from casatasks import (listobs, delmod, flagdata, setjy, bandpass, gaincal,
#                        applycal, blcal, fluxscale)
# from casaplotms import plotms

# Measurement set name
###################
vis = 'r64.rx230.ms'
###################

# Run listobs to identify the calibrator and target field names
# to input below.
output = listobs(vis)
#, listfile=vis+'.listobs', overwrite=True)



# Setting fields and spectral windows
# Flux cal (usually a Solar system object)
flux = 'URANUS'
# Bandpass cal
bpcal = '3C84'
# Gain/phase calibrator (usually 2)
bothpcal = '0102+584,0112+227,0120-270,0217+017,0423-013,0958+655,3C273,3C345'

calfields=bpcal+','+bothpcal+','+flux

# Specify the field names
# If you have a mosaic, using a wildcard will most likely work (e.g. M51*)
science='ACTGALAXY,LSI61303,NGC1600,NGC3894,NGC4472,SIS22-G4-T1,SIS22-G4-T2,SIS22-G4-T3,SIS22-G4-T4,SIS22-G4-T5,SIS22-G4-T6'

# Spectral window range to consider during the calibration steps.
bpchans='0~11'
calchans='0~11'
spwrange='0~11' # i.e. '1~24' or '1~48'

# Scan number of long bandpass scan (typically 10s of min)
bpscan='1201'

# Set the reference antenna. Check the observing log
# to ensure the antenna was in the array for the observation!
# By definition, this sets phase = 0.
refant='AN06'

# Setting edge channels and bandpass smoothing channels
# All SMA observations using SWARM (>2017-2018) have 16384 channels
# per SPW.
# Data sets are often rechunked to average down the number of channels
# as this can vastly reduce the data volume.
# If you're unsure of the rechunk factor, use the # channels per SPW
# in the listobs output to calculate:
# rechunk = 16384 / Nchan
rechunk_factor = 64

# We will clip off the edges of each SPW where the response drops off
# This will remove 2.5% on the upper and lower sides
chan_frac = 0.025
num_chan = 16384 // rechunk_factor
low_chan_edge = int(num_chan * chan_frac)
high_chan_edge = int(num_chan * (1 - chan_frac))

# Edge channels to flag out from the above calculation
edgechan='0~11:0~{0};{1}~{2}'.format(low_chan_edge, high_chan_edge, num_chan)

# This sets two things for calculating the bandpass solution:
# time duration, channels to smooth over
# Time duration 'inf' combines all times together (standard for max S/N)
# Channels to smooth over forces coherency between some adjacent channels
# and can boost the S/N if needed.
# Use a lower # channels for larger rechunk factors.
bpsolint='inf,2ch'


# Setting summary plot output
# CASA uses plotms to view averaged visibility data for examination
# This can be tedious to view many plots so the above setting enables
# saving a large set of png plots that can be reviewed after the script runs.
# NOTE: this can take hours to run! It's helpful to run overnight or in the background.
enable_summary_plots = True


###############
# Initial flagging.
###############

# Ensure zeros are flagged
flagdata(vis=vis, mode='clip', clipzeros=True, flagbackup=False)

# Flag spw edges (see above)
flagdata(vis=vis,spw=edgechan, flagbackup=False)
# This backs up a flagging version that could later be restored.
flagmanager(vis=vis, mode='save', versionname='flag_edges')

# Next, use plotms interactively to identify "birdies"
# There are 1-2 channels with large amplitude due to harmonics
# in the local oscillator frequency.
# List the flagging commands below

# Uncomment below to use plotms to loop through each
# calibrator to identify the spw and channel number with a spike

# First look for outliers in amp vs channel
# for myfield in calfields.split(","):
#     plotms(vis=vis, xaxis='channel',
#         yaxis='amp',field=myfield, avgtime='1e8', avgscan=False,
#         coloraxis='ant1',iteraxis='spw', ydatacolumn='data',
#         gridrows=4, gridcols=3, yselfscale=True, showgui=True)
#     input(f"Done adding freq. flags for {myfield}?")  # FOR python 3 (CASA 6)
# #     raw_input("Done adding freq. flags for {0}?".format(myfield))  # For python 2 (CASA 5)

# # Next look for time outliers
# for myfield in calfields.split(","):
#     plotms(vis=vis, xaxis='time',
#         yaxis='amp',field=myfield, avgchannel='1e8',
#         coloraxis='ant1',iteraxis='spw', ydatacolumn='data',
#         gridrows=4, gridcols=3, yselfscale=True, showgui=True)
#     input(f"Done adding time flags for {myfield}?")  # FOR python 3 (CASA 6)
# #     raw_input("Done adding time flags for {0}?".format(myfield))  # For python 2 (CASA 5)

# Flag birdies: add SPW and channel number following the examples here
# These persist typically over all fields. By not specifying a field name,
# we flag over all of them.

# flagdata(vis=vis, mode='manual', spw='2:36')
# flagdata(vis=vis, mode='manual', spw='4:78')
# flagdata(vis=vis, mode='manual', spw='10:76')
# flagdata(vis=vis, mode='manual', spw='7:25')
# flagdata(vis=vis, mode='manual', spw='8:34')

# Example of flagging a time period with poor data:
# flagdata(vis=vis, mode='manual', field='3C345',
#          timerange='14:55:00~14:59:00')


# After calibration, additional flags can be added below here
# and the script re-run to recalibrate the data with all flags applied.

#flagdata(vis=vis,spw='0:8189~8193,4:8189~8193,2:7427~7429,6:6644~6648;7426~7428')
#flagdata(vis=vis,spw='0:10384,2:9710~9760',antenna='AN03,AN08')
#flagdata(vis=vis,spw='2:6620~6680,6:6620~6680,2:7420~7440,6:7420~7440',antenna='AN03,AN08')

#flagdata(vis=vis,spw='',timerange='06:36:00~06:50:00')
#flagdata(vis=vis,spw='',timerange='07:28:30~07:35:00')
#flagdata(vis=vis,spw='',scan='60,116')

####################################################################
# CALIBRATION
####################################################################

# To clear explicit setting of quasar fluxes below (if rerunning script)
delmod(vis,otf=True,scr=True)

# Use a standard for setting the absolute flux scale.
# Use updated Solar System models
setjy(vis=vis, field=flux, standard='Butler-JPL-Horizons 2012',
      scalebychan=False)



# Phase-only selfcal on the BP
# We're attempting to correct for short time fluctuations on the scale
# of individual time integrations. This will help boost the S/N when solving
# for the bandpass
os.system('rm -rf *bpself.gcal')
gaincal(vis=vis,caltable=vis+'.bpself.gcal',
        field=bpcal,spw=bpchans,refant=refant,scan=bpscan,
        calmode='p',solint='int',minsnr=3.0,minblperant=3)

# Plot phase vs time per antenna colored by SPW
# Check that the phase variations are smooth in time.
# Offsets between antennas are OK.
plotms(vis=vis+'.bpself.gcal',
       xaxis='time',
       yaxis='phase',
       coloraxis='spw',
       iteraxis='antenna',
       ydatacolumn='data',
#        plotfile="{0}.png".format(vis+'.bpself.gcal'),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)

# NOTE: in all plotms calls, uncomment the plotfile input to save
# plots of each calibration table.

# Solnorm True bandpass
os.system('rm -rf *solnorm_true.bcal')
# Smooth some channels to help produce a smooth solution over channels
# This will produce 1 solution per SPW in both amplitude and phase.
# `bpsolint` (defined at the top), sets 2 things: time-average,channel-average
# In this case we have bpsolint="inf,2ch" meaning we will integrate together ALL times
# and average over 2 channels.
bandpass(vis=vis,caltable=vis+'.bandpass.solnorm_true.bcal',
         bandtype='B',scan=bpscan,
         field=bpcal,spw=spwrange,combine='scan',refant=refant,
         solint=bpsolint,solnorm=True,minblperant=3,
         gaintable=[vis+'.bpself.gcal'])


# Plot bandpass phase and amplitude solutions vs. freq.
# There will be bumps+wiggles but it should be smoothly varying
# Note if there are large outliers in phase or amp! We may need to add flagging
# and re-do if necessary.

# Each panel shows 1 antenna and the colors are the different SPWs
plotms(vis=vis+'.bandpass.solnorm_true.bcal',xaxis='freq',
       yaxis='phase', coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
#        plotfile="{0}.phase.png".format(vis+'.bandpass.solnorm_true.bcal'),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)

plotms(vis=vis+'.bandpass.solnorm_true.bcal',xaxis='freq',
       yaxis='amp', coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
#        plotfile="{0}.amp.png".format(vis+'.bandpass.solnorm_true.bcal'),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)

# We can apply the bandpass table to the bandpass in the MS.
# This allows us to check whether the amp vs. freq is now roughly
# flat across all SPWs.
applycal(vis=vis,field=bpcal,
         spw=spwrange,
         gaintable=[vis+'.bandpass.solnorm_true.bcal'],
         interp='nearest',
         spwmap=[],
         gainfield=bpcal,
         flagbackup=False,calwt=False)

# Make amp+phase selfcal table applying the bandpass
# This is re-doing the self-cal on the bandpass to setup for an extra step
# with "blcal" (per-baseline calibrations) that the SMA data often needs.
os.system('rm -rf *bpself.ap.gcal')
gaincal(vis=vis,caltable=vis+'.bpself.ap.gcal',
        field=bpcal,spw=bpchans,refant=refant,scan=bpscan,
        calmode='ap',solint='int',minsnr=3.0,minblperant=3,
        gaintable=[vis+'.bandpass.solnorm_true.bcal'])

# Same as the last phase self-cal on the bandpass, we want things to
# be smoothly varying in time.
# Plot phase vs time per antenna colored by SPW
plotms(vis=vis+'.bpself.ap.gcal',
       xaxis='time',
       yaxis='phase',
       coloraxis='spw',
       iteraxis='antenna',
       ydatacolumn='data',
#        plotfile="{0}.png".format(vis+'.bpself.ap.gcal'),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)

# Plot amp vs time per antenna colored by SPW
plotms(vis=vis+'.bpself.ap.gcal',
       xaxis='time',
       yaxis='amp',
       coloraxis='spw',
       iteraxis='antenna',
       ydatacolumn='data',
#        plotfile="{0}.png".format(vis+'.bpself.ap.gcal'),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)

# blcal
os.system('rm -rf *blcal')

# SMA data from the above steps can have offsets between different SPW.
# For consistent amplitude calibration, we can enforce per-baseline solutions.
# amplitude calibration (vs. per-antenna that all other steps use).
blcal(vis=vis, caltable=vis+'.blcal',
      field=bpcal, spw=spwrange,solint='inf',scan=bpscan,
      solnorm=False,freqdep=False,
      gaintable=[vis+'.bpself.ap.gcal',vis+'.bandpass.solnorm_true.bcal'],
      calmode='a')

# Plot amp vs freq per antenna with points colored by the SPW
plotms(vis=vis+'.blcal',xaxis='freq',
       yaxis='amp',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
#        plotfile="{0}.amp.png".format(vis+'.blcal'),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)


###############
# GAINS
###############

# Last steps of the calibration: using the gain calibrators.

# We start with a short time interval (1 time integration) period.
# These will be applied in the next step to smooth out short time variations
# (similar to what we did at the start for the bandpass) when solving for the
# amplitude gains.

# int phase with combining spws per side-band

os.system('rm -rf *intphase_combinespw.gcal')

# Solve per sideband
# Each receiver has a lower and upper sideband in frequency.
# We will use the entire sideband (so 6 SPWs) for each solution.

# LSB
# The lower side-band in our data is spw 0 to 5 (0~5 to select all)
gaincal(vis=vis,caltable=vis+'.intphase_combinespw.gcal',
        field=calfields,refant=refant,
        combine='spw',spw='0~5',
        calmode='p',solint='int',minsnr=3.0,minblperant=3,
        gaintable=[vis+'.bandpass.solnorm_true.bcal',
                   vis+'.blcal'])

# USB
# The upper side-band in our data is spw 6 to 11 (6~11 to select all)
# NOTE: we can append to the same table using append=True.
gaincal(vis=vis,caltable=vis+'.intphase_combinespw.gcal',append=True,
        field=calfields,refant=refant,
        combine='spw',spw='6~11',
        calmode='p',solint='int',minsnr=3.0,minblperant=3,
        gaintable=[vis+'.bandpass.solnorm_true.bcal',
                   vis+'.blcal'])

# Plot phase vs time per antenna
# As before, we're checking for large phase jumps within a single scan.
# Small variations are fine as the S/N is lower over a single time interval.
plotms(vis=vis+'.intphase_combinespw.gcal',xaxis='time',
       yaxis='phase',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
#        plotfile="{0}.png".format(vis+'.intphase_combinespw.gcal'),
       gridrows=4, gridcols=2,
       yselfscale=True, showgui=True, overwrite=True)

# Next we repeat the same, but now averaging over an entire scan on
# the gain calibrators.
# This should boost the S/N since we're integrating over more data.
# Same idea as above with applying to the LSB and USB.

# scan phase with combining spws (long solint)
os.system('rm -rf *scanphase_combinespw.gcal')
# LSB
gaincal(vis=vis,caltable=vis+'.scanphase_combinespw.gcal',
        field=calfields,refant=refant,
        combine='spw',spw='0~5',
        calmode='p',solint='inf',minsnr=3.0,minblperant=3,
        gaintable=[vis+'.bandpass.solnorm_true.bcal',
                   vis+'.blcal'])
# USB
gaincal(vis=vis,caltable=vis+'.scanphase_combinespw.gcal',append=True,
        field=calfields,refant=refant,
        combine='spw',spw='6~11',
        calmode='p',solint='inf',minsnr=3.0,minblperant=3,
        gaintable=[vis+'.bandpass.solnorm_true.bcal',
                   vis+'.blcal'])

# Plot phase vs time per antenna colored by SPW
# We're look for small phase variations between gain calibrator scans.
# Large jumps in phase may mean that, even with the calibration applied,
# we may not get a coherent solution on our targets! In that case, we may
# need to flag the time period during the large phase jump.
plotms(vis=vis+'.scanphase_combinespw.gcal',xaxis='time',
       yaxis='phase',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
#        plotfile="{0}.png".format(vis+'.scanphase_combinespw.gcal'),
       gridrows=4, gridcols=2,
       yselfscale=True, showgui=True, overwrite=True)

# Gain amplitude calibration.
# We ALSO do this per sideband, and this is the repeated
# [0,0,0...,6,6,6,6] you see below.
# The 0s correspond to spws 0-5 (the LSB) and the 6s to 6-11 (the USB).
# We also apply the scanphase_combinespw.gcal table
# (per-integration phase from above) to correct for rapid phase
# variations that will help boost the S/N for the amplitude solutions.

os.system('rm -rf *amp.gcal')
gaincal(vis=vis,caltable=vis+'.amp.gcal',
        field=calfields,spw=calchans,refant=refant,
        combine='',
        spwmap=[[],[],[0,0,0,0,0,0,6,6,6,6,6,6]],
        calmode='ap',solint='10min',minsnr=3.0,minblperant=3,
        gaintable=[vis+'.bandpass.solnorm_true.bcal',
                   vis+'.blcal',
                   vis+'.scanphase_combinespw.gcal'])

# Plot amp vs time per antenna colored by SPW
# As with the phase, we want smooth variations in time.
plotms(vis=vis+'.amp.gcal', xaxis='time',
       yaxis='amp',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
#        plotfile="{0}.png".format(vis+'.amp.gcal'),
       gridrows=4, gridcols=2,
       yselfscale=True, showgui=True, overwrite=True)


# Lastly, we need to use our flux calibrator model (set with setjy at the
# beginning) to our other calibrator sources so they have the correct flux
# density.
# NOTE: the data already have the Tsys correction applied so they will be
# close to the actual flux density already, but underestimated by
# ~20%(ish) due to unmodelled efficiency losses.

os.system('rm -rf *flux.cal')
fluxresults = {}

fluxresults[vis]=fluxscale(vis=vis,caltable=vis+'.amp.gcal',
                           refspwmap=[-1],transfer=[bpcal,bothpcal],
                           fluxtable=vis+'.flux.cal',reference=flux)
print(fluxresults[vis])

# While you can explore the log or flux table to see the calibrator fluxes per SPW,
# the dictionary in fluxresults can also be handy to review.
# One way to retain the dictionary is to save it as a numpy binary file:
import numpy as np
np.save(vis+".fluxscale.npy", fluxresults, allow_pickle=True)

# To load, use:
# fluxresults = np.load(vis+".fluxscale.npy").item()

####################################################################
# Applycal
####################################################################

# We have all the calibration tables we need. Now we apply the necessary
# ones to the data to create our calibrated science data.

# Make a backup so we can compare before/after calibration (optional)
flagmanager(vis=vis,mode='save',versionname='beforeapplycal')

# Apply the calibrations to the bandpass calibrators:
# NOTE: same usage for per-sideband gain phase solutions as above
# (list of 0s and 6s). We also set the type of interpolation to convert
# from the calibration tables to our sources. In most cases we use
# "nearest" to use the gain calibration nearest in time.
applycal(vis=vis,field=bpcal,
         spw=spwrange,
         gaintable=[vis+'.bandpass.solnorm_true.bcal',
                    vis+'.blcal',
                    vis+'.scanphase_combinespw.gcal',
                    vis+'.flux.cal'],
         interp=['nearest','nearest','nearestPD','nearest'],
         spwmap=[[],[],[0,0,0,0,0,0,6,6,6,6,6,6],[]],
         gainfield=[bpcal,bpcal,bpcal,bpcal],flagbackup=False,calwt=False)

# Same thing, but now looping over the gain calibrators.
for this_pcal in bothpcal.split(","):

    applycal(vis=vis,field=this_pcal,
            spw=spwrange,
            gaintable=[vis+'.bandpass.solnorm_true.bcal',
                       vis+'.blcal',
                       vis+'.scanphase_combinespw.gcal',
                       vis+'.flux.cal'],
            interp=['nearest','nearest','nearestPD','nearest'],
            spwmap=[[],[],[0,0,0,0,0,0,6,6,6,6,6,6],[]],
            # NOTE: here we need to use the bandpass source for the BP
            # related calibration tables
            # and the phase cal for the phase calibrator tables.
            gainfield=[bpcal, bpcal, this_pcal, this_pcal],
            flagbackup=False, calwt=False)

# Lastly, we apply the calibration to our science targets. This can be
# done all at once because the inputs to gainfield does not change
# (but we use all the gain calibrators for the gain tables)
applycal(vis=vis,field=science,
         spw=spwrange,
         gaintable=[vis+'.bandpass.solnorm_true.bcal',
                    vis+'.blcal',
                    vis+'.scanphase_combinespw.gcal',
                    vis+'.flux.cal'],
         interp=['nearest','nearest','linearPD','linear'],
         spwmap=[[],[],[0,0,0,0,0,0,6,6,6,6,6,6],[]],
         gainfield=[bpcal,bpcal,bothpcal,bothpcal],
         flagbackup=False,calwt=False)


# Additional optional plotting:
# plotms(vis=vis,xaxis='Time',yaxis='phase',
#        field=pcal1,avgchannel='999999',
#        avgfield=True,coloraxis='baseline',iteraxis='antenna',ydatacolumn='corrected',
#        plotfile=vis+'.aftercal.gainpha.png',gridrows=3,gridcols=2)

# plotms(vis=vis,xaxis='freq',
#          yaxis='amp',field='0',avgtime='1e8',avgscan=True,
#          coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
#          plotfile=vis+'.beforecal.specamp.png',gridrows=3,gridcols=2)

# plotms(vis=vis,xaxis='freq',
#          yaxis='amp',field='0',avgtime='1e8',avgscan=True,
#          coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
#          plotfile=vis+'.aftercal.specamp.png',gridrows=3,gridcols=2)

# plotms(vis=vis,xaxis='freq',
#          yaxis='phase',field='0',avgtime='1e8',avgscan=True,
#          coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
#          plotfile=vis+'.beforecal.specpha.png',gridrows=3,gridcols=2)

# plotms(vis=vis,xaxis='freq',
#          yaxis='phase',field='0',avgtime='1e8',avgscan=True,
#          coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
#          plotfile=vis+'.aftercal.specpha.png',gridrows=3,gridcols=2)


#plotms(vis=vis,yaxis='amp',xaxis='time',avgchannel='16384',spw='',coloraxis='field',iteraxis='spw',
#       ydatacolumn='corrected')

#plotms(vis=vis,yaxis='amp',xaxis='freq',avgtime='1e8',avgscan=True,field='3',coloraxis='spw',iteraxis='spw',
#       ydatacolumn='corrected')


# Additional field summary plots added by Eric:

####################################################################
# Examine calibrated data
####################################################################

# Loop through and save plots of the calibrators in:
# 1. Amp vs Freq (time averaged) per antenna
# 2. Amp vs Time (freq averaged) per antenna
# 3. Amp vs uv-dist (time + freq averaged) per SPW
# 4. Phase vs Freq (time averaged) per antenna
# 5. Phase vs Time (freq averaged) per antenna
# 6. Phase vs uv-dist (time + freq averaged) per SPW
# 7. Real vs Imag (time + freq averaged)
# 8. Amp vs baseline (time + freq averaged)
# 9. Phase vs baseline (time + freq averaged)

# For comparison, we will make two of each plot: before and after calibration.

# NOTE: this can take a long time to create all plots! The code is written to
# run non-interactively with png plots created to more rapidly review and
# identify issues that can be then be examined interactively in plotms.


if enable_summary_plots:

    plot_path = "calibrator_plots"
    if not os.path.exists(plot_path):
        os.mkdir(plot_path)

    scitarget_plot_path = "target_plots"
    if not os.path.exists(scitarget_plot_path):
        os.mkdir(scitarget_plot_path)

    # Calibrator plots:

    myvis=vis
    for myfield in calfields.split(","):

        # Amp/Phase vs freq
        # After calibration:
        plotms(vis=myvis,xaxis='freq',
                yaxis='amp',field=myfield,avgtime='1e8',avgscan=True,
                coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.aftercal.amp_freq.png".format(plot_path,
                                                                    myvis,
                                                                    myfield.lower()))
        plotms(vis=myvis,xaxis='freq',
                yaxis='phase',field=myfield,avgtime='1e8',avgscan=True,
                coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.aftercal.pha_freq.png".format(plot_path,
                                                                    myvis,
                                                                    myfield.lower()))

        # Before calibration:
        plotms(vis=myvis,xaxis='freq',
                yaxis='amp',field=myfield,avgtime='1e8',avgscan=True,
                coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.beforecal.amp_freq.png".format(plot_path,
                                                                    myvis,
                                                                    myfield.lower()))
        plotms(vis=myvis,xaxis='freq',
                yaxis='phase',field=myfield,avgtime='1e8',avgscan=True,
                coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.beforecal.pha_freq.png".format(plot_path,
                                                                    myvis,
                                                                    myfield.lower()))

        # Amp/Phase vs time
        # After calibration:
        plotms(vis=myvis,xaxis='time',
                yaxis='amp', field=myfield, avgchannel='1e8', avgscan=False,
                coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.aftercal.amp_time.png".format(plot_path,
                                                                    myvis,
                                                                    myfield.lower()))
        plotms(vis=myvis,xaxis='time',
                yaxis='phase',field=myfield, avgchannel='1e8', avgscan=False,
                coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.aftercal.pha_time.png".format(plot_path,
                                                                    myvis,
                                                                    myfield.lower()))

        # Before calibration:
        plotms(vis=myvis,xaxis='time',
                yaxis='phase',field=myfield, avgchannel='1e8', avgscan=False,
                coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.beforecal.amp_time.png".format(plot_path,
                                                                    myvis,
                                                                    myfield.lower()))
        plotms(vis=myvis,xaxis='time',
                yaxis='phase',field=myfield, avgchannel='1e8', avgscan=False,
                coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.beforecal.pha_time.png".format(plot_path,
                                                                    myvis,
                                                                    myfield.lower()))

        # Amp/Phase vs uvdist
        plotms(vis=myvis,xaxis='uvdist',
                yaxis='amp', field=myfield, avgchannel='8', avgscan=False, avgtime='1e8',
                coloraxis='ant1',iteraxis='spw',ydatacolumn='corrected',
                gridrows=4, gridcols=3, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.aftercal.amp_uvdist.png".format(plot_path,
                                                                    myvis,
                                                                    myfield.lower()))
        plotms(vis=myvis,xaxis='uvdist',
                yaxis='phase',field=myfield, avgchannel='8', avgscan=False, avgtime='1e8',
                coloraxis='ant1',iteraxis='spw',ydatacolumn='corrected',
                gridrows=4, gridcols=3, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.aftercal.pha_uvdist.png".format(plot_path,
                                                                    myvis,
                                                                    myfield.lower()))

        # Real vs Imag
        plotms(vis=myvis,xaxis='imag',
                yaxis='real',field=myfield, avgchannel='8', avgscan=False, avgtime='1e8',
                coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                plotfile="{0}/{1}.{2}.aftercal.amp_pha.png".format(plot_path,
                                                                myvis,
                                                                myfield.lower()))

        # Amp/Phase vs baseline
        plotms(vis=myvis,xaxis='baseline',
                yaxis='amp', field=myfield, avgchannel='1e8', avgscan=False, avgtime='1e8',
                coloraxis='spw',iteraxis=None, ydatacolumn='corrected',
                gridrows=1, gridcols=1, showgui=False, overwrite=True,
                plotfile="{0}/{1}.{2}.aftercal.amp_baseline.png".format(plot_path,
                                                                        myvis,
                                                                        myfield.lower()))
        plotms(vis=myvis,xaxis='baseline',
                coloraxis='spw',iteraxis=None, ydatacolumn='corrected',
                gridrows=1, gridcols=1, showgui=False, overwrite=True,
                plotfile="{0}/{1}.{2}.aftercal.pha_baseline.png".format(plot_path,
                                                                        myvis,
                                                                        myfield.lower()))

    # Target plots:

    # Now examine our science targets. Only amplitude is typically helpful here
    # as complex sources structure will not have easily interpretable phases!


    for myfield in science.split(","):

        # Amp vs freq
        plotms(vis=myvis,xaxis='freq',
                yaxis='amp',field=myfield,avgtime='1e8',avgscan=True,
                coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.aftercal.amp_freq.png".format(scitarget_plot_path,
                                                                    myvis,
                                                                    myfield.lower()))

        # Amp vs time
        plotms(vis=myvis,xaxis='time',
                yaxis='amp', field=myfield, avgchannel='1e8', avgscan=False,
                coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.aftercal.amp_time.png".format(scitarget_plot_path,
                                                                    myvis,
                                                                    myfield.lower()))

        # Amp vs uvdist
        plotms(vis=myvis,xaxis='uvdist',
                yaxis='amp', field=myfield, avgchannel='1e8', avgscan=False, avgtime='1e8',
                coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
                gridrows=4, gridcols=2, showgui=False, overwrite=True,
                yselfscale=True,
                plotfile="{0}/{1}.{2}.aftercal.amp_uvdist.png".format(scitarget_plot_path,
                                                                    myvis,
                                                                    myfield.lower()))
