import os
import sys

#************************************************************************
#************************************************************************
#*****                         DATA                                 *****
#************************************************************************
#************************************************************************

# Input root: folder where nightly folders are stored
input_root = '/data/mash/marusa/2m3data/wifes/'

# Output root: where to save nightly data
output_root = '/data/mash/marusa/2m3reduced/wifes/'

# Save to folders with this prefix
prefix=None

# Minimal number of each of the calibration frames. Default is 3.
calmin = 3

# Object list: reduce only these objects
#~ object_list_filename = '/data/mash/marusa/2m3data/wifes/reduce_these_objects.dat'
#~ object_list = ['RZ Mic']

# Run numbers of bad files that you don't want to include in the reduction
excluderun_filename = '/data/mash/marusa/2m3data/wifes/list_of_bad_exposures_that_we_shouldn_use.dat'

# List of bad calibration files
badcalib_filename = '/data/mash/marusa/2m3data/wifes/list_of_high_biases_pay_attention.dat' # None


#************************************************************************
#************************************************************************
#*****                USER REDUCTION DESIGN IS SET HERE             *****
#************************************************************************
#************************************************************************

band = 'r' # RedBand

# Don't show plots. Just save them into diagnostics folder so reduction doesn't need any interaction if not asked for. Verbose.

# Co-add images of the same object
coadd_images = True
# If false, then separate them in the metadata file. Take care with arcs.

# SET MULTITHREAD ?
multithread=False

# SET SKIP ALREADY DONE FILES ?
skip_done=True

proc_steps = [
    # All steps up to p08 were deleted.
    #------------------
    {'step':'extract_stars'  , 'run':True, 'suffix':None,
     'args':{'ytrim':4, 
             'type':'flux'}},
    {'step':'derive_calib'   , 'run':True, 'suffix':None,
     'args':{'plot_stars':False,
             'plot_sensf':False,
             'polydeg':10,
             'method':'smooth_SG', # 'poly' or 'smooth_SG'
             'boxcar':10, # smoothing for smooth_SG only
             'norm_stars':True}},
    {'step':'flux_calib'     , 'run':True, 'suffix':'09', 'args':{}},
    #------------------
    {'step':'extract_stars'  , 'run':True, 'suffix':None,
     'args':{'ytrim':4, 
             'type':'telluric'}},
    {'step':'derive_telluric', 'run':True, 'suffix':None,
     'args':{'plot':False}},
    {'step':'telluric_corr'  , 'run':True, 'suffix':'10', 'args':{}},
    #------------------
    {'step':'save_3dcube'    , 'run':True, 'suffix':'11', 'args':{}}
    #~ #------------------
    ]
