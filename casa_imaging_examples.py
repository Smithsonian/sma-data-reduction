
##################################
# Imaging in CASA
##################################

# Script adapted from a tutorial written by Luca Matra

# Imaging in CASA uses the tclean task:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.imaging.tclean.html
# This does all of the steps for producing "cleaned" image:
# 1. Gridding the uv-data and taking an inverse-FFT to make the dirty image
# 2. Deconvolution of the dirty map until some stopping threshold is reached
# 3. Restoration of the cleaned map from the clean model components + residual field

# Below are examples of the two types of imaging we will use: continuum and spectral line.


# Add the name of your MS file here!
vis_name = "FILENAME.ms.target"


# Select your target name here.
# We need to tell CASA which source in the data that we want to make an image of
fieldname = 'Arp148'


#######################
#CONTINUUM
#######################

# Continuum imaging needs to avoid frequencies dominated by strong spectral line emission
# as this will bias the inferred continuum flux.

# plotms can be used to search for bright line emission

#To search for and manually (within plotms) flag strong lines.
# WATCH for line vs atmospheric feature! If you know a weak line may be present,
# you can flag around the frequency where you expect it to be.
# Note that you can check for atmospheric features locations using the online
# SMA Passband Visualizer.

# For the tuning of this data set, it include the lines 12CO, 13CO, and C18(2-1).
# To identify these lines, we can plot an averaged spectrum to pick out lines
# bright enough to show up in the visibilities.:

# For example, to identify the 12CO(2-1) line we need to know:
# 1. the rest frequency of 12CO(2-1) (230.538 GHz)
# 2. the line-of-sight velocity range of our science target
# 3. Identify the corresponding SPW and channel range these frequencies correspond to
#    using the plotms command below.

plotms(vis=vis_name, field=fieldname,
       xaxis='Frequency', yaxis='amp',
	   freqframe='LSRK',
       avgtime='1e20', avgscan=True, coloraxis='spw',
	   )

# Use the format spw:chan1~chan2 (e.g. 7:800~1000)
spw_freqstring =

# Also record the min. max freq those correspond to (we will need this for the spectral line
# imaging).
# Input as a string or float: e.g. "230.538"
freqmin =
freqmax =

# We want to flag this range to avoid it being included with the continuum, but we WANT to keep this
# range for the spectral line imaging below.
# We'll make a copy of the MS here and then apply the flagging.

# We make a copy of the data here so that we can remove the
# spectra line emission so it does not bias the continuum image.
import os
os.system("cp -r {0} {0}.cont".format(vis_name))

vis_name_cont = "{}.cont".format(vis_name)

flagdata(vis=vis_name_cont, mode='manual', spw=spw_freqstring)

# Time to start the actual imaging process with tclean.

# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.imaging.tclean.html

# COPY AND PASTE THE PARAMETER SETTINGS INTO YOUR CASA TERMINAL

# File names
# Make the filename based on the name of the target we're imaging.
# 'natural' refers to how the data are gridded (see 'weighting' setting below).
imagename = f'{fieldname}_continuum' #Output images

# Imaging parameters
#Shows residual image and waits after every major cycle iteration.
# Number of minor cycle iterations will be set manually in viewer.
interactive = True
# Number of iterations before stopping deconvolution (global stopping criterion).
niter = 1000
#Pixel size of image in arcsec.
# Needs to be a small fraction (<1/3) of expected beam size (few arcsec in this case).
cell = '0.5arcsec'
# Typically about the size of the primary beam (55" for SMA at 1.3mm), measured in
# number of pixels. Better if it is a power of 2, as it makes the FFT algorithm more efficient.
imsize = 512
#Weighting of the visibilities before imaging. Natural gives larger beam, lowest noise (hence max SNR).
# Uniform gives smallest beam, highest noise. Briggs is something in between depending on robust parameter.
weighting = 'robust'
# Only needed if choosing 'briggs' weighting. Balances between natural (+2) and uniform (-2)
# in increments of 0.5.
robust = 0.5

#Continuum parameters
#For continuum imaging, use multi-frequency synthesis (mfs) mode.
# This combines data from all frequencies into our image, which maximizes the sensitivity
specmode = 'mfs'

# Gridding and CLEAN algorithm choice. All of these, to begin with, are standard inputs,
# so we do not actually need to input them, but we will here for completeness.
gridder = 'standard'
# Classic Hogbom 1974, but modified with major and minor cycles.
deconvolver = 'hogbom'
# Standard loop gain gamma which is multiplied by the maximum intensity of the residual image
# and added to the model.
gain = 0.1
# Sets a threshold in Jy/beam. e.g., "5.0mJy/beam". We don't set this as we will set nsigma to estimate it.
threshold = ""
# Set global threshold for the residual image max in nsigma*rms to stop iterations
nsigma = 3.0
# Max number of minor cycle iterations per major cycle.
# Set to -1 initially as we will decide iteratively in interactive mode.
cycleniter = -1
# Used to determine minor cycle threshold. Factor multiplied by the maximum dirty beam
# sidelobe level to calculate when to trigger major cycle.
cyclefactor = 1.0
# Used to determine minor cycle threshold. If max dirty beam sidelobe level is less than
# this, use 5% as a threshold to trigger major cycle. Lower boundary for major cycle trigger.
minpsffraction = 0.05
# Used to determine minor cycle threshold. If max dirty beam sidelobe level is more than this,
# use 80% as a threshold to trigger major cycle. Upper boundary for major cycle trigger.
maxpsffraction = 0.8
# Primary beam limit sets the size of the field where valid data is included in the field-of-view
# The primary beam size is set by the antenna size (7 m for SMA antennas).
# Roughly speaking, the noise level goes as 1 / pb decreasing radially outward.
pblimit = 0.25


#Remove images if they already exist
import os
os.system('rm -rf ./'+imagename+'.*')

#Run tclean
tclean(vis=vis_name_cont,
	   imagename=imagename,
	   field=fieldname,
	   interactive=interactive,
	   niter=niter,
	   cell=cell,
	   imsize=imsize,
	   weighting=weighting,
	   robust=robust,
	   gridder=gridder,
	   deconvolver=deconvolver,
	   gain=gain,
	   threshold=threshold,
	   nsigma=nsigma,
	   cycleniter=cycleniter,
	   cyclefactor=cyclefactor,
	   minpsffraction=minpsffraction,
	   maxpsffraction=maxpsffraction,
	   specmode=specmode)

# Every time viewer comes up (initially to show you the produced dirty image after gridding, weighting,
# and inverting, and subsequently to show you residual image after every major cycle), you want to select
# a region (create it, then double click within it) within which you believe real source signal to lie.
# CLEAN will find the maximum (over and over for cycleniter minor cycles) only within this region.
# The region can be modified as you CLEAN deeper. Stop (hit the X button) when the maximum in the residual
# image is a few times the RMS noise level in the image!
# Note, initially and at every major cycle, you have to select the number of *further* iterations
# in the 'max cycleniter' box of viewer. At the same time, you have to keep
# 'iterations left' >> 'max cycleniter', so make sure a very large number in that box the first time
# the viewer comes up. This 'iterations left' value won't affect anything, but it will quit your CLEAN
# process if <= cycleniter.

# Finally, use the CASA viewer to check out your image, which will be the .image file.
# Note that dirty beam (.psf file), primary beam (.pb file), residual (.res), model (.model) can
# all be accessed after CLEANing.

# Uncomment out to open the viewer GUI.
# viewer()


#######################
#LINE IMAGING
# Example for the 12CO(2-1) line
#######################


# Here we can use the vis_name that we defined above:
# vis_name

# Similar to avoiding the line emission for the continuum, we also want to correct for the
# continuum. Instead of remove a region, we will fit a model in the visibility space (or uv-space)
# and subtract it off.


#Now carry out continuum subtraction in u-v space through uvcontsub task, where this is done for each time interval and baseline.
# The continuum is first fit at frequencies outside of the freqmin-to-freqmax range, and a polynomial of degree:
fitorder = 0
# is fitted and then subtracted from visibilities.
# NOTE: the freq ranges were set above
fitspw = '*:'+str(freqmin)+'~'+str(freqmax)+'GHz'
# * will select all SPWs

#IMPORTANT to exclude rather than select frequency range in fitspw for continuum fitting.
excludechans = True
#Then, run continuum subtraction
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.manipulation.uvcontsub.html
uvcontsub(vis=vis_name, fitspw=fitspw, fitorder=fitorder,
		  excludechans=excludechans)

# NOTE: if this command fails, use "uvcontsub_old" instead:
# uvcontsub_old(vis=vis_name, fitspw=fitspw, fitorder=fitorder,
# 		  excludechans=excludechans)

# uvcontsub creates a new MS ending in ".contsub" (this cannot be changed)
vis_name_contsub = "{}.contsub".format(vis_name)
# Sanity check that the MS exists:
assert os.path.exists(vis_name_contsub)

# Time to start the actual imaging process with tclean.

##########
# Strategy to image our data:
##########
# NOTE: our targets have extended CO emission (i.e., not just point sources).
# Because of this, we want to only deconvolve the emission from the source,
# BUT there are sidelobes from the telescope's PSF that make this confusing.
# Instead of imaging in a single command, we can make it easier to identify
# the real emission by splitting the deconvolution ('cleaning') into two steps:

# Step 1. Deconvolve the emission starting without a clean mask to a brightness level
#    5x the noise
# Step 2. Then use interactive clean to draw a clean mask and deconvolve to a lower level
#    (3x the noise is a good level to reach).

# Why 2 steps? The strong and confusing sidelobe structures will be minimized by
# deconvolving the brightest sources. This will make it easier to identify the real emission
# in step 2. We want to avoid including the sidelobe structure in the clean model as this
# will give false emission in different areas of our map. Our goal is to reconstruct the
# true image on the sky.

# File names
imagename = f'{fieldname}_12co21' # Output images

# Imaging parameters
#Shows residual image and waits after every major cycle iteration.
# Number of minor cycle iterations will be set manually in viewer.
interactive = True
# Number of iterations before stopping deconvolution (global stopping criterion).
niter = 1000
#Pixel size of image in arcsec.
# Needs to be a small fraction (<1/3) of expected beam size (few arcsec in this case).
cell = '0.5arcsec'
# Typically about the size of the primary beam (55" for SMA at 1.3mm), measured in
# number of pixels. Better if it is a power of 2, as it makes the FFT algorithm more efficient.
imsize = 512
#Weighting of the visibilities before imaging. Natural gives larger beam, lowest noise (hence max SNR).
# Uniform gives smallest beam, highest noise. Briggs is something in between depending on robust parameter.
weighting = 'robust'
# Only needed if choosing 'briggs' weighting. Balances between natural (+2) and uniform (-2)
# in increments of 0.5.
robust = 0.5

# Line parameters
# Note that this is a lot different than the continuum case!
# This is where recording the velocity range above is useful! Choose the start, width and nchan to cover
# this velocity range.
specmode = 'cube' #For line imaging, use cube mode.
# Channel width chosen for output cube. CLEAN will carry out interpolation to resample the visibility data
# before imaging. '' for native (see listobs spw channel width). Use deltanu/nu_line = deltav/c to figure
# out velocity widths from frequency widths.
width = '30.0km/s' # I suggest this as a good starting point for the velocity channel width***

# NOTE: this will depend on your target!!
start = '0km/s' # velocity of starting channel of cube. Make sure to cover whole line!

nchan = 50 # number of channels in cube. Make sure to cover whole line!
outframe = 'LSRK' #output reference frame for velocities This is local-standard-of-rest kinematic
restfreq = '230.538GHz' #Rest frequency of line of interest, in this case CO J=2-1.


# Gridding and CLEAN algorithm choice. All of these, to begin with, are standard inputs,
# so we do not actually need to input them, but we will here for completeness.
gridder = 'standard'

# Our emission is extended, and so we use the multi-scale deconvolution algorithm.
# The standard clean methods assume the emission is a collection of point-sources,
# which is a poor approximation when we have extended emission.
deconvolver = 'multiscale'
# We also need to supply which scales to use for multiscale clean.
# These values are in pixel, and a good rule-of-thumb is [0, beam, 3*beam]
# where 0 is a point source, and the others are the scales of the resolution.
# The actual values we give are in pixel units. From `cell=0.5arcsec` and the
# sensitivity calculator we used for the proposal, we expect the resolution to be
# ~3-4arcsec. So in pixels, 1 beam is about 6 pixels across.
# In pixels, [0, beam, 3*beam] corresponds to:
scales = [0, 6, 18]

# Standard loop gain gamma which is multiplied by the maximum intensity of the residual image
# and added to the model.
gain = 0.1
# Sets a threshold in Jy/beam. e.g., "5.0mJy/beam". We don't set this as we will set nsigma to estimate it.
threshold = ""
# Set global threshold for the residual image max in nsigma*rms to stop iterations
nsigma = 5.0

# Max number of minor cycle iterations per major cycle.
# Set to -1 initially as we will decide iteratively in interactive mode.
cycleniter = -1
# Used to determine minor cycle threshold. Factor multiplied by the maximum dirty beam
# sidelobe level to calculate when to trigger major cycle.
cyclefactor = 1.0
# Used to determine minor cycle threshold. If max dirty beam sidelobe level is less than
# this, use 5% as a threshold to trigger major cycle. Lower boundary for major cycle trigger.
minpsffraction = 0.05
# Used to determine minor cycle threshold. If max dirty beam sidelobe level is more than this,
# use 80% as a threshold to trigger major cycle. Upper boundary for major cycle trigger.
maxpsffraction = 0.8
# Primary beam limit sets the size of the field where valid data is included in the field-of-view
# The primary beam size is set by the antenna size (7 m for SMA antennas).
# Roughly speaking, the noise level goes as 1 / pb decreasing radially outward.
pblimit = 0.25

#Remove images if they already exist
import os

if os.path.exists("{}.residual".format(imagename)):
	os.system('rm -r {}.*'.format(imagename))

# Step 1. Deconvolve the emission starting without a clean mask to a brightness level
#    5x the noise

#Run tclean
tclean(vis=vis_name_contsub,
	   imagename=imagename,
	   field=fieldname,
	   interactive=False,
	   niter=niter,
	   cell=cell,
	   imsize=imsize,
	   weighting=weighting,
	   robust=robust,
	   gridder=gridder,
	   deconvolver=deconvolver,
	   gain=gain,
	   threshold=threshold,
	   nsigma=nsigma,
	   cycleniter=cycleniter,
	   cyclefactor=cyclefactor,
	   minpsffraction=minpsffraction,
	   maxpsffraction=maxpsffraction,
	   specmode=specmode,
	   width=width,
	   start=start,
	   nchan=nchan,
	   restfreq=restfreq,
	   outframe=outframe)


# Step 2. Then use interactive clean to draw a clean mask and deconvolve to a lower level
#    (3x the noise is a good level to reach).

#Repeat procedure done for the continuum, but now mask drawing and minor and major cycles are carried out
# for every channel of the cube! Ideally want to create a different mask for each channel, and remove masks
# (at different times for different channels) as emission is CLEANed down to the noise level.

# Re-running tclean will restart from the previous call, including the initial
# deconvolution done in step 1.

nsigma = 3.

tclean(vis=vis_name_contsub,
	   imagename=imagename,
	   field=fieldname,
	   interactive=True,
	   niter=niter,
	   cell=cell,
	   imsize=imsize,
	   weighting=weighting,
	   robust=robust,
	   gridder=gridder,
	   deconvolver=deconvolver,
	   gain=gain,
	   threshold=threshold,
	   nsigma=nsigma,
	   cycleniter=cycleniter,
	   cyclefactor=cyclefactor,
	   minpsffraction=minpsffraction,
	   maxpsffraction=maxpsffraction,
	   specmode=specmode,
	   width=width,
	   start=start,
	   nchan=nchan,
	   restfreq=restfreq,
	   outframe=outframe,
	   calcres=False, calcpsf=False,  # These save some time to avoid recalculating saved products
	   )



#Finally, use the CASA viewer to check out your cube, which will be the .image file. Note that dirty beam
# (.psf file), primary beam (.pb file), residual (.res), model (.model) can all be accessed after CLEANing.
# viewer()
