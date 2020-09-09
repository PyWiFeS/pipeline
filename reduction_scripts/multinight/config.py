"""
Parameters for data reduction of the WiFeS data.

Processing steps:
(1) Edit this file.
(2) Generate metadata.
(3) Run reduction.
(4) Process stellar.
(5) Cleanup.
"""

import os
import sys

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# root needed only for output
output_root = '/priv/mulga1/marusa/2m3reduced/wifes/'

# Do you want to reduce only specific objects? Names must match those in the fits file headers (OBJNAME).
objectnames=None
exclude_objectnames=['PDS 70']

# Do you want to reduce only images with specific binning?
ccdsum=None #'1 1' # '1 2' # binning; False or None
naxis2=None #2056 # False # 2056 for PDS 70

# Save to folders with this prefix
prefix=None#'ys'

# This thing with metadata_filename is actually not used yet.
if prefix is not None and len(prefix)>0:
    metadata_filename='%s_metadata'%prefix
else:
    metadata_filename='metadata'

#------------------------------------------------------------------------

generate_metadata={'prefix': prefix,
                    'CCDSUM': ccdsum,
                    'objectnames': objectnames,
                    'exclude_objectnames': exclude_objectnames,
                    'metadata_filename': metadata_filename,
                    'output_root': output_root,
                    }

#------------------------------------------------------------------------
#------------------------------------------------------------------------
#************************************************************************
#*****                USER REDUCTION DESIGN IS SET HERE             *****
#************************************************************************

# Don't show plots. Just save them into diagnostics folder so reduction doesn't need any interaction if not asked for. Verbose.

# Co-add images of the same object
coadd_images = True
# If false, then separate them in the metadata file. Take care with arcs.

# SET MULTITHREAD ?
multithread=False

# SET SKIP ALREADY DONE FILES ?
skip_done=True

proc_steps = [
    #------------------
    {'step':'overscan_sub'   , 'run':True, 'suffix':'00', 'args':{}},
    {'step':'bpm_repair'     , 'run':True, 'suffix':'01', 'args':{}},
    #------------------
    {'step':'superbias'      , 'run':True, 'suffix':None,
     'args':{'method':'row_med', 
             'plot':False, 
             'verbose':False}},
    {'step':'bias_sub'       , 'run':True, 'suffix':'02',
     'args':{'method':'subtract', 
             'plot':False, 
             'verbose':False}},
    #------------------
    {'step':'superflat'      , 'run':True, 'suffix':None,
     'args':{'source':'dome'}},
#~ #    {'step':'superflat'      , 'run':False, 'suffix':None,
#~ #     'args':{'source':'twi', 
#~ #             'scale':'median_nonzero'}},
    {'step':'slitlet_profile', 'run':True, 'suffix':None, 'args':{}},
    #------------------
    {'step':'flat_cleanup'   , 'run':True, 'suffix':None,
     #~ 'args':{'type':['dome','twi'], 
     'args':{'type':['dome'], 
             'verbose':True, 
             'plot':False,
             'buffer':4,
             'offsets':[0.4,0.4],
             'radius':10.0,
             'nsig_lim':3.0}},
    #------------------
    {'step':'superflat_mef'  , 'run':True, 'suffix':None,
     'args':{'source':'dome'}},
#~ #    {'step':'superflat_mef'  , 'run':False, 'suffix':None,
#~ #     'args':{'source':'twi'}},
    #------------------
    {'step':'slitlet_mef'    , 'run':True, 'suffix':'03',
     'args':{'ns':False}},
    #------------------
    {'step':'wave_soln'      , 'run':True, 'suffix':None,
     'args':{'verbose':True,
             'method' : 'optical',
             'doalphapfit' : True,
             'doplot' : False, #['step2'], # True, False, or ['step1','step2']
             'shift_method' : 'xcorr_all',
             'find_method' : 'mpfit',
             'dlam_cut_start':5.0,
             'multithread': multithread}},
#    {'step':'wire_soln'      , 'run':True, 'suffix':None, 'args':{}},
    {'step':'flat_response'  , 'run':True, 'suffix':None,
#     'args':{'mode':'all'}},
     'args':{'mode':'dome'}},
    #------------------
    {'step':'cosmic_rays'    , 'run':True, 'suffix':'04',
     'args':{'ns':False, 
             'multithread':multithread}},
    #------------------
    {'step':'sky_sub'        , 'run':True, 'suffix':'05',
     'args':{'ns':False}},
    #------------------
    {'step':'obs_coadd'      , 'run': True, 'suffix':'06',
     'args':{'method':'sum'}},
    #------------------
    {'step':'flatfield'      , 'run':True, 'suffix':'07', 'args':{}},
    #------------------             
    {'step':'cube_gen'       , 'run':True, 'suffix':'08',
     'args':{'multithread':multithread,
             'adr':True,
             'dw_set': (0.44 if band=='r' else 0.77),
             'wmin_set': (5400.0 if band=='r' else 3500.0), 
             'wmax_set': (7000.0 if band=='r' else 5700.0)}},     
    #------------------
    #~ {'step':'extract_stars'  , 'run':True, 'suffix':None,
     #~ 'args':{'ytrim':4, 
             #~ 'type':'flux'}},
    #~ {'step':'derive_calib'   , 'run':True, 'suffix':None,
     #~ 'args':{'plot_stars':True,
             #~ 'plot_sensf':True,
             #~ 'polydeg':10,
             #~ 'method':'smooth_SG', # 'poly' or 'smooth_SG'
             #~ 'boxcar':10, # smoothing for smooth_SG only
             #~ 'norm_stars':True}},
    #~ {'step':'flux_calib'     , 'run':True, 'suffix':'09', 'args':{}},
    #------------------
    #~ {'step':'extract_stars'  , 'run':True, 'suffix':None,
     #~ 'args':{'ytrim':4, 
             #~ 'type':'telluric'}},
    #~ {'step':'derive_telluric', 'run':True, 'suffix':None,
     #~ 'args':{'plot':True}},
    #~ {'step':'telluric_corr'  , 'run':True, 'suffix':'10', 'args':{}},
    #------------------
    #~ {'step':'save_3dcube'    , 'run':True, 'suffix':'11', 'args':{}}
    #~ #------------------
    ]
