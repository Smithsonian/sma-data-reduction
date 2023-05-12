# Based on Todd Hunter's CASA calibration script,
# modified by Charlie Qi 05/2021
# revised on 12/2021 for flux cal
# The measurement set is generated with CASA 5.7.2-4.

import os

# CASA 5
from tasks import (listobs, delmod, flagdata, plotms,
                   setjy, bandpass, gaincal, applycal, blcal, fluxscale)

# CASA 6
# from casatasks import (listobs, delmod, flagdata, setjy, bandpass, gaincal,
#                        applycal, blcal, fluxscale)
# from casaplotms import plotms

###################
vis='210528.rx230.ms'
###################

listobs(vis,listfile=vis+'.listobs')

#Fields: 8
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    A    3C279               12:56:11.166000 -05.47.21.52400 J2000   0          16740
#  1         VESTA               11:06:02.977798 +14.45.46.32823 J2000   1           3600
#  2    A    1517-243            15:17:41.813000 -24.22.19.47500 J2000   2           7200
#  3    A    1626-298            16:26:06.021000 -29.51.26.97100 J2000   3          11232
#  4         IMLUP               15:56:09.206712 -37.56.06.12616 J2000   4         105240
#  5         CALLISTO            22:14:09.885093 -11.51.21.91465 J2000   5           1800
#  6         MWC349A             20:32:45.540000 +40.39.36.61100 J2000   6           1800
#  7    A    3C454.3             22:53:57.748000 +16.08.53.56300 J2000   7          21780
#Spectral Windows:  (12 unique spectral windows and 1 unique polarization setups)
#  SpwID  Name   #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz)  Corrs
#  0      none   16384   TOPO  203046.960       139.648  -2288000.0 201903.0299   RR
#  1      none   16384   TOPO  198746.993       139.648   2288000.0 199890.9235   RR
#  2      none   16384   TOPO  199047.026       139.648  -2288000.0 197903.0961   RR
#  3      none   16384   TOPO  194747.060       139.648   2288000.0 195890.9898   RR
#  4      none   16384   TOPO  195047.093       139.648  -2288000.0 193903.1624   RR
#  5      none   16384   TOPO  190747.126       139.648   2288000.0 191891.0560   RR
#  6      none   16384   TOPO  210746.794       139.648   2288000.0 211890.7244   RR
#  7      none   16384   TOPO  215046.762       139.648  -2288000.0 213902.8317   RR
#  8      none   16384   TOPO  214746.728       139.648   2288000.0 215890.6579   RR
#  9      none   16384   TOPO  219046.696       139.648  -2288000.0 217902.7657   RR
#  10     none   16384   TOPO  218746.661       139.648   2288000.0 219890.5914   RR
#  11     none   16384   TOPO  223046.630       139.648  -2288000.0 221902.6997   RR


flux='5'
bpcal='0,7'
pcal1='2'
pcal2='3'
bothpcal=pcal1
calfields=bpcal+','+bothpcal+','+flux
science='4' # include all targets;
bpchans='0~11'
calchans='0~11'
spwrange='0~11' # i.e. '1~24' or '1~48'
nspws=7 # NB: this needs to be one more than number of chunks to cover psuedo-cont spw 0
bpscan='1,701' # scan number of long BP scan

refant='AN02'

edgechan='0~11:0~1039;15344~16383' # edge channel to flag out
bpsolint='inf,64ch' # bandpass smoothing window per 16 channel

###############
# BANDPASS
###############

# Flag edges and birdies
flagdata(vis=vis,spw=edgechan)

# Flag birdies
#flagdata(vis=vis,spw='0:8189~8193,4:8189~8193,2:7427~7429,6:6644~6648;7426~7428')
#flagdata(vis=vis,spw='0:10384,2:9710~9760',antenna='AN03,AN08')
#flagdata(vis=vis,spw='2:6620~6680,6:6620~6680,2:7420~7440,6:7420~7440',antenna='AN03,AN08')

# Aftercal
#flagdata(vis=vis,spw='',timerange='06:36:00~06:50:00')
#flagdata(vis=vis,spw='',timerange='07:28:30~07:35:00')
#flagdata(vis=vis,spw='',scan='60,116')

####################################################################
# CALIBRATION
####################################################################

# To clear explicit setting of quasar fluxes below (if rerunning script)
delmod(vis,otf=True,scr=True)

# use updated Solar System models
setjy(vis=vis,field=flux,standard='Butler-JPL-Horizons 2012',scalebychan=False)



# phase-only selfcal
os.system('rm -rf *bpself.gcal')
gaincal(vis=vis,caltable=vis+'.bpself.gcal',
        field=bpcal,spw=bpchans,refant=refant,scan=bpscan,
        calmode='p',solint='int',minsnr=3.0,minblperant=3)

# Plot phase vs time per antenna colored by SPW
plotms(vis=vis+'.bpself.gcal,
       xaxis='time',
       yaxis='phase',
       coloraxis='spw',
       iteraxis='antenna',
       ydatacolumn='data',
       plotfile="{0}.png".format(vis+'.bpself.gcal),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)


# Solnorm True bandpass
os.system('rm -rf *solnorm_true.bcal')
# smooth some channels CHECK SMOOTHING WINDOW
bandpass(vis=vis,caltable=vis+'.bandpass.solnorm_true.bcal',
         bandtype='B',scan=bpscan,
         field=bpcal,spw=spwrange,combine='scan',refant=refant,
         solint=bpsolint,solnorm=True,minblperant=3,
         gaintable=[vis+'.bpself.gcal'])
# B can't be combined with combine='spw'

#plotbandpass(caltable=vis+'.bandpass.solnorm_true.bcal',
#             field=bpcal,xaxis='chan',yaxis='amp',
#             interactive=False,subplot=42,vis=vis,figfile='figures/'+vis+'.bandpass.amp')

# Plot bandpass phase and amplitude solutions vs. freq.
# Each panel shows 1 antenna and the colors are the different SPWs
plotms(vis=vis+'.bandpass.solnorm_true.bcal',xaxis='freq',
       yaxis='phase', coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       plotfile="{0}.phase.png".format(vis+'.bandpass.solnorm_true.bcal'),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)
plotms(vis=vis+'.bandpass.solnorm_true.bcal',xaxis='freq',
       yaxis='amp', coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       plotfile="{0}.amp.png".format(vis+'.bandpass.solnorm_true.bcal'),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)


# try 'BPOLY' ?
applycal(vis=vis,field=bpcal,
         spw=spwrange,
         gaintable=[vis+'.bandpass.solnorm_true.bcal'],
         interp='nearest',
         spwmap=[],
         gainfield=bpcal,
         flagbackup=False,calwt=False)

# Make ap selfcal table applying the bandpass
os.system('rm -rf *bpself.ap.gcal')
gaincal(vis=vis,caltable=vis+'.bpself.ap.gcal',
        field=bpcal,spw=bpchans,refant=refant,scan=bpscan,
        calmode='ap',solint='int',minsnr=3.0,minblperant=3,
        gaintable=[vis+'.bandpass.solnorm_true.bcal'])

# Plot phase vs time per antenna colored by SPW
plotms(vis=vis+'.bpself.ap.gcal,
       xaxis='time',
       yaxis='phase',
       coloraxis='spw',
       iteraxis='antenna',
       ydatacolumn='data',
       plotfile="{0}.png".format(vis+'.bpself.ap.gcal),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)
# Plot amp vs time per antenna colored by SPW
plotms(vis=vis+'.bpself.ap.gcal,
       xaxis='time',
       yaxis='amp',
       coloraxis='spw',
       iteraxis='antenna',
       ydatacolumn='data',
       plotfile="{0}.png".format(vis+'.bpself.ap.gcal),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)

# blcal
os.system('rm -rf *blcal')

blcal(vis=vis, caltable=vis+'.blcal',
      field=bpcal, spw=spwrange,solint='inf',scan=bpscan,
      solnorm=False,freqdep=False,
      gaintable=[vis+'.bpself.ap.gcal',vis+'.bandpass.solnorm_true.bcal'],
      calmode='a')

# Plot amp vs freq per antenna with points colored by the SPW
plotms(vis=vis+'.blcal',xaxis='freq',
       yaxis='amp',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       plotfile="{0}/{1}.amp.png".format(vis+'.blcal'),
       gridrows=4, gridcols=2,
       yselfscale=True,
       showgui=True, overwrite=True)

# Not sure about blcal table, but we needed to do this with the old correlator
# to avoid spw to spw offsets, so we continue to use it.

###############
# Phase up table for phase combine
###############


###############
# GAINS
###############

# int phase with combining spws per side-band

os.system('rm -rf *intphase_combinespw.gcal')

# Solve per sideband
# LSB
gaincal(vis=vis,caltable=vis+'.intphase_combinespw.gcal',
        field=calfields,refant=refant,
        combine='spw',spw='0~5',
        calmode='p',solint='int',minsnr=3.0,minblperant=3,
        gaintable=[vis+'.bandpass.solnorm_true.bcal',
                   vis+'.blcal'])

# USB
# NOTE: we can append to the same table using append=True
gaincal(vis=vis,caltable=vis+'.intphase_combinespw.gcal',append=True,
        field=calfields,refant=refant,
        combine='spw',spw='6~11',
        calmode='p',solint='int',minsnr=3.0,minblperant=3,
        gaintable=[vis+'.bandpass.solnorm_true.bcal',
                   vis+'.blcal'])

# Plot phase vs time per antenna
plotms(vis=vis+'.intphase_combinespw.gcal',xaxis='time',
       yaxis='phase',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       plotfile="{0}.png".format(vis+'.intphase_combinespw.gcal'),
       gridrows=4, gridcols=2,
       yselfscale=True, showgui=True, overwrite=True)

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
plotms(vis=vis+'.scanphase_combinespw.gcal',xaxis='time',
       yaxis='phase',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       plotfile="{0}.png".format(vis+'.scanphase_combinespw.gcal'),
       gridrows=4, gridcols=2,
       yselfscale=True, showgui=True, overwrite=True)

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
plotms(vis=vis+'.amp.gcal', xaxis='time',
       yaxis='amp',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       plotfile="{0}.png".format(vis+'.amp.gcal),
       gridrows=4, gridcols=2,
       yselfscale=True, showgui=True, overwrite=True)

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

# final applycals with bootstrapped flux cal table


flagmanager(vis=vis,mode='save',versionname='beforeapplycal')


applycal(vis=vis,field=bpcal,
         spw=spwrange,
         gaintable=[vis+'.bandpass.solnorm_true.bcal',
                    vis+'.blcal',
                    vis+'.scanphase_combinespw.gcal',
                    vis+'.flux.cal'],
         interp=['nearest','nearest','nearestPD','nearest'],
         spwmap=[[],[],[0,0,0,0,0,0,6,6,6,6,6,6],[]],
         gainfield=[bpcal,bpcal,bpcal,bpcal],flagbackup=False,calwt=False)


applycal(vis=vis,field=pcal1,
         spw=spwrange,
         gaintable=[vis+'.bandpass.solnorm_true.bcal',
                    vis+'.blcal',
                    vis+'.scanphase_combinespw.gcal',
                    vis+'.flux.cal'],
         interp=['nearest','nearest','nearestPD','nearest'],
         spwmap=[[],[],[0,0,0,0,0,0,6,6,6,6,6,6],[]],
         gainfield=[bpcal,bpcal,pcal1,pcal1],flagbackup=False,calwt=False)



applycal(vis=vis,field=pcal2,
         spw=spwrange,
         gaintable=[vis+'.bandpass.solnorm_true.bcal',
                    vis+'.blcal',
                    vis+'.scanphase_combinespw.gcal',
                    vis+'.flux.cal'],
         interp=['nearest','nearest','nearestPD','nearest'],
         spwmap=[[],[],[0,0,0,0,0,0,6,6,6,6,6,6],[]],
         gainfield=[bpcal,bpcal,pcal1,pcal1],flagbackup=False,calwt=False)


applycal(vis=vis,field=science,
         spw=spwrange,
         gaintable=[vis+'.bandpass.solnorm_true.bcal',
                    vis+'.blcal',
                    vis+'.scanphase_combinespw.gcal',
                    vis+'.flux.cal'],
         interp=['nearest','nearest','linearPD','linear'],
         spwmap=[[],[],[0,0,0,0,0,0,6,6,6,6,6,6],[]],
         gainfield=[bpcal,bpcal,bothpcal,bothpcal],flagbackup=False,calwt=False)





plotms(vis=vis,xaxis='Time',yaxis='phase',
       field=pcal1,avgchannel='999999',
       avgfield=True,coloraxis='baseline',iteraxis='antenna',ydatacolumn='corrected',
       plotfile=vis+'.aftercal.gainpha.png',gridrows=3,gridcols=2)



plotms(vis=vis,xaxis='freq',
         yaxis='amp',field='0',avgtime='1e8',avgscan=True,
         coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
         plotfile=vis+'.beforecal.specamp.png',gridrows=3,gridcols=2)


plotms(vis=vis,xaxis='freq',
         yaxis='amp',field='0',avgtime='1e8',avgscan=True,
         coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
         plotfile=vis+'.aftercal.specamp.png',gridrows=3,gridcols=2)

plotms(vis=vis,xaxis='freq',
         yaxis='phase',field='0',avgtime='1e8',avgscan=True,
         coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
         plotfile=vis+'.beforecal.specpha.png',gridrows=3,gridcols=2)


plotms(vis=vis,xaxis='freq',
         yaxis='phase',field='0',avgtime='1e8',avgscan=True,
         coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
         plotfile=vis+'.aftercal.specpha.png',gridrows=3,gridcols=2)



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

enable_summary_plots = False

if enable_summary_plots:

    plot_path = "calibrator_plots"
    if not os.path.exists(plot_path):
        os.mkdir(plot_path)

    scitarget_plot_path = "target_plots"
    if not os.path.exists(scitarget_plot_path):
        os.mkdir(scitarget_plot_path)

    # Calibrator plots:

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
                plotfile="{0}/{1}.{2}'.aftercal.amp_baseline.png".format(plot_path,
                                                                        myvis,
                                                                        myfield.lower()))
        plotms(vis=myvis,xaxis='baseline',
                coloraxis='spw',iteraxis=None, ydatacolumn='corrected',
                gridrows=1, gridcols=1, showgui=False, overwrite=True,
                plotfile="{0}/{1}.{2}'.aftercal.pha_baseline.png".format(plot_path,
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

