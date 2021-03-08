#Follow the online instructions to get data from MIR-printed UVFITS into CASA:
import sys 
sys.path.append('data/HD127821/.') 
import MIRFITStoCASA 
fullvis='data/HD127821/HD127821_track6.ms' #Change name to something appropriate
allNames = []
for sou in ['HD127821']: #Change name to exact source name from MIR (beginning characters of UVFITS files)
	for rx in ['230','240']: #Select applicable receivers
		for sb in ['L','U']: #Select applicable sidebands
			for i in ['1','2','3','4']: #Select applicable chunks
				name = sou+"_"+sb+"_S"+i+"_RX"+rx
				print("------converting "+name+" ....")
				MIRFITStoCASA.MIRFITStoCASA(UVFITSname=name+'.UVFITS', MSname=name+'.ms')  
				allNames.append(name+'.ms')
				
#Check CASA MS files containing visibilities are there and well:
!ls -d ./*.ms

#Now concatenate all chunks, sidebands and receivers into one big visibility MS dataset
concat(vis=allNames,concatvis=fullvis,timesort=True) #(Safely ignore warnings about negative or zero total bandwidth)

#Now inspect structure of data in the active CASA log (integration times, spectral windows and channels within them, sources, antennas etc)
listobs(fullvis)

#Then flag (bad) edge channels, typically for nflagedge you want to use the same value as you did for 'ntrim' in MIR. ntotal can be found under '#Chans' in listobs
nflagedge=130
ntotal=2048
#Execute flagging
flagdata(vis=fullvis, mode='manualflag', spw='*:0~'+str(nflagedge)+';'+str(ntotal-nflagedge)+'~'+str(ntotal-1))

#Inspect visibilities with plotms task. Did we flag enough, are there any outliers? 
#To plot time-averaged amplitudes as a function of channel, iterating over each spectral window (spw, i.e. chunk):
plotms(vis=fullvis, xaxis='channel', yaxis='amp', avgtime='1e20', avgscan=True, iteraxis='spw')  #Remember to check every spw!
#To plot full u-v coverage:
plotms(vis=fullvis, xaxis='uwave', yaxis='vwave', coloraxis='Baseline')
#To plot Real(V) after averaging in time and frequency as a function of u-v distance (radius in u-v space):
plotms(vis=fullvis, xaxis='uvdist', yaxis='Real', avgtime='1e20', avgscan=True, avgchannel='1e20')



#######################
#CONTINUUM
#######################

#All looks good. For continuum imaging, duplicate the MS to produce a dataset that we will use solely for continuum:
contvis=fullvis[:-3]+'_cont.ms'
split(vis=fullvis, outputvis=contvis, keepflags=False, datacolumn='data')

#Then use plotms, e.g.:
plotms(vis=contvis, xaxis='Frequency', yaxis='amp', freqframe='LSRK', avgtime='1e20', avgscan=True, coloraxis='spw')
#To search for and manually (within plotms) flag strong lines. WATCH for line vs atmospheric feature! If you know a weak line may be present, you can flag around the frequency where you expect it to be.
#Note that you can check for atmospheric features locations using the online SMA Passband Visualizer.

#Great, now average in frequency so that imaging is faster. For the SMA in subcompact configuration, it's safe enough to average all channels to 1 channel per SPW:
contvisavg=contvis[:-3]+'_chanavg.ms'
width=[2048,2048,2048,2048,2048,2048,2048, 2048] #This is the factor you want to average by, for each chunk. The number of channels you end up with is ntotal/width per chunk. Make sure ntotal/width ends up being an integer!
split(vis=contvis, outputvis=contvisavg, keepflags=False, datacolumn='data', width=width)
#Check that this worked in listobs:
listobs(contvisavg)

#Now, if you had several observations of this target, you would repeat the above for each observations and concatenate here.
#I am just going to supply you with the concatenated dataset containing 6 different SMA tracks in both compact and subcompact configurations:
contvisavg_alltracks='data/HD127821_calibratedvis_cont.ms'



#Time to start the actual imaging process with tclean. Use 'inp tclean' to verify inputs before running, and 'help tclean' to clarify what is needed as input.
#File names
imagename='continuum/HD127821/HD127821_natural_cont' #Output images
vis=contvisavg_alltracks #Input visibility data

#Imaging parameters
interactive=True #Shows residual image and waits after every major cycle iteration. Number of minor cycle iterations will be set manually in viewer.
niter=10000000 #Number of iterations before stopping deconvolution (global stopping criterion). Set to a very large number, as we will stop manually before this in interactive mode.
cell=['0.1arcsec'] #Pixel size of image in arcsec. Needs to be a small fraction (<1/3) of expected beam size (few arcsec in this case).
imsize=[512,512] #Typically about the size of the primary beam (55" for SMA at 1.3mm), measured in number of pixels. Better if it is a power of 2, as it makes the FFT algorithm more efficient.
weighting='natural' #Weighting of the visibilities before imaging. Natural gives larger beam, lowest noise (hence max SNR). Uniform gives smallest beam, highest noise. Briggs is something in between depending on robust parameter.
robust=0.5 #Only needed if choosing 'briggs' weighting.
uvtaper=[''] #Do you want to further taper the data in u-v space to lower the resolution below what is delivered by natural weighting? If so, specify tapering in arcsec.
#Continuum parameters
specmode='mfs' #For continuum imaging, use multi-frequency synthesis mode.

#gridding and CLEAN algorithm choice. All of these, to begin with, are standard inputs, so we do not actually need to input them, but we will here for completeness.
gridder='standard'
deconvolver='hogbom' #Classic Hogbom 1974, but modified with major and minor cycles. 
gain=0.1 #Standard loop gain gamma which is multiplied by the maximum intensity of the residual image and added to the model.
threshold=0.0 #Do not use a global threshold for the residual image max in Jy to stop iterations, stop manually in interactive mode.
nsigma=0.0 #Do not use a global threshold for the residual image max in nsigma*rms to stop iterations, stop manually in interactive mode.
cycleniter=-1 #Max number of minor cycle iterations per major cycle. Set to -1 initially as we will decide iteratively in interactive mode.
cyclefactor=1.0 #used to determine minor cycle threshold. Factor multiplied by the maximum dirty beam sidelobe level to calculate when to trigger major cycle.
minpsffraction=0.05 #used to determine minor cycle threshold. If max dirty beam sidelobe level is less than this, use 5% as a threshold to trigger major cycle. Lower boundary for major cycle trigger.
maxpsffraction=0.8 #used to determine minor cycle threshold. If max dirty beam sidelobe level is more than this, use 80% as a threshold to trigger major cycle. Upper boundary for major cycle trigger.

#Remove images if they already exist
import os
os.system('rm -rf ./'+imagename+'.*')

#Run tclean
tclean(vis=vis, imagename=imagename, interactive=interactive, niter=niter, cell=cell, imsize=imsize, weighting=weighting, robust=robust, uvtaper=uvtaper, gridder=gridder, deconvolver=deconvolver, gain=gain, threshold=threshold, nsigma=nsigma, cycleniter=cycleniter, cyclefactor=cyclefactor, minpsffraction=minpsffraction, maxpsffraction=maxpsffraction, specmode=specmode)
#Every time viewer comes up (initially to show you the produced dirty image after gridding, weighting, and inverting, and subsequently to show you residual image after every major cycle), you want to select a region (create it, then double click within it) within which you believe real source signal to lie. CLEAN will find the maximum (over and over for cycleniter minor cycles) only within this region. The region can be modified as you CLEAN deeper. Stop (hit the X button) when the maximum in the residual image is a few times the RMS noise level in the image!
#Note, initially and at every major cycle, you have to select the number of *further* iterations in the 'max cycleniter' box of viewer. At the same time, you have to keep 'iterations left' >> 'max cycleniter', so make sure a very large number in that box the first time the viewer comes up. This 'iterations left' value won't affect anything, but it will quit your CLEAN process if <= cycleniter.

#Finally, use the CASA viewer to check out your image, which will be the .image file. Note that dirty beam (.psf file), primary beam (.pb file), residual (.res), model (.model) can all be accessed after CLEANing.
viewer




#######################
#LINE IMAGING
#######################

#Since there was no line detection in the previous dataset, we start from a different CASA-imported SMA dataset, where I have already split-ted out the chunk (spw) containing the CO J=2-1 line, at the full SMA resolution.
linevis='data/disk/disk_track1_CO.ms'

#Run a quick listobs to check it out
listobs(linevis)

#Then use plotms, e.g.:
plotms(vis=linevis, xaxis='Frequency', yaxis='amp', avgtime='1e20', avgscan=True)
#To search for the wanted line, and note frequency range where line is significantly above the continuum. If line is too weak, choose a generous range around the expected line frequency.
freqmin=230.530 #in GHz
freqmax=230.536 #in GHz
#Now carry out continuum subtraction in u-v space through uvcontsub task, where this is done for each time interval and baseline.
#The continuum is first fit at frequencies outside of the freqmin-to-freqmax range, and a polynomial of degree:
fitorder=1 
#is fitted and then subtracted from visibilities.
fitspw='*:'+str(freqmin)+'~'+str(freqmax)+'GHz'
excludechans=True #IMPORTANT to exclude rather than select frequency range in fitspw for continuum fitting.
#Then, run continuum subtraction
uvcontsub(vis=linevis, fitspw=fitspw, fitorder=1, excludechans=excludechans)

#Now, if you had several observations of this target, you would repeat the above for each observations and concatenate here (or first concatenate, then do the above)
#Here we just use this date so nothing to do:
linevisavg_alltracks=linevis+'.contsub'

#Time to start the actual imaging process with tclean. Use 'inp tclean' to verify inputs before running, and 'help tclean' to clarify what is needed as input.
#File names
imagename='line/disk/disk_briggs00_CO21' #Output images
vis=linevisavg_alltracks #Input visibility data

#Imaging parameters
interactive=True #Shows residual image and waits after every major cycle iteration. Number of minor cycle iterations will be set manually in viewer.
niter=10000000 #Number of iterations before stopping deconvolution (global stopping criterion). Set to a very large number, as we will stop manually before this in interactive mode.
cell=['0.1arcsec'] #Pixel size of image in arcsec. Needs to be a small fraction (<1/3) of expected beam size (few arcsec in this case).
imsize=[512,512] #Typically about the size of the primary beam (55" for SMA at 1.3mm), measured in number of pixels. Better if it is a power of 2, as it makes the FFT algorithm more efficient.
weighting='briggs' #Weighting of the visibilities before imaging. Natural gives larger beam, lowest noise (hence max SNR). Uniform gives smallest beam, highest noise. Briggs is something in between depending on robust parameter.
robust=0.0 #Only needed if choosing 'briggs' weighting.
uvtaper=[''] #Do you want to further taper the data in u-v space to lower the resolution below what is delivered by natural weighting? If so, specify tapering in arcsec.

#Line parameters
specmode='cube' #For line imaging, use cube mode.
width='0.4km/s' #Channel width chosen for output cube. CLEAN will carry out interpolation to resample the visibility data before imaging. '' for native (see listobs spw channel width). Use deltanu/nu_line=deltav/c to figure out velocity widths from frequency widths.
start='0km/s' #velocity of starting channel of cube. Make sure to cover whole line!
nchan=35 #number of channels in cube. Make sure to cover whole line!
outframe='LSRK' #output reference frame for velocities
restfreq='230.538e9' #Rest frequency of line of interest, in this case CO J=2-1.

#gridding and CLEAN algorithm choice. All of these, to begin with, are standard inputs, so we do not actually need to input them, but we will here for completeness.
gridder='standard'
deconvolver='hogbom' #Classic Hogbom 1974, but modified with major and minor cycles. 
gain=0.1 #Standard loop gain gamma which is multiplied by the maximum intensity of the residual image and added to the model.
threshold=0.0 #Do not use a global threshold for the residual image max in Jy to stop iterations, stop manually in interactive mode.
nsigma=0.0 #Do not use a global threshold for the residual image max in nsigma*rms to stop iterations, stop manually in interactive mode.
cycleniter=-1 #Max number of minor cycle iterations per major cycle. Set to -1 initially as we will decide iteratively in interactive mode.
cyclefactor=1.0 #used to determine minor cycle threshold. Factor multiplied by the maximum dirty beam sidelobe level to calculate when to trigger major cycle.
minpsffraction=0.05 #used to determine minor cycle threshold. If max dirty beam sidelobe level is less than this, use 5% as a threshold to trigger major cycle. Lower boundary for major cycle trigger.
maxpsffraction=0.8 #used to determine minor cycle threshold. If max dirty beam sidelobe level is more than this, use 80% as a threshold to trigger major cycle. Upper boundary for major cycle trigger.

#Remove images if they already exist
import os
os.system('rm -rf ./'+imagename+'.*')

#Run tclean
tclean(vis=vis, imagename=imagename, interactive=interactive, niter=niter, cell=cell, imsize=imsize, weighting=weighting, robust=robust, uvtaper=uvtaper, gridder=gridder, deconvolver=deconvolver, gain=gain, threshold=threshold, nsigma=nsigma, cycleniter=cycleniter, cyclefactor=cyclefactor, minpsffraction=minpsffraction, maxpsffraction=maxpsffraction, specmode=specmode, width=width, start=start, nchan=nchan, restfreq=restfreq, outframe=outframe)
#Repeat procedure done for the continuum, but now mask drawing and minor and major cycles are carried out for every channel of the cube! Ideally want to create a different mask for each channel, and remove masks (at different times for different channels) as emission is CLEANed down to the noise level.

#Finally, use the CASA viewer to check out your cube, which will be the .image file. Note that dirty beam (.psf file), primary beam (.pb file), residual (.res), model (.model) can all be accessed after CLEANing.
viewer

