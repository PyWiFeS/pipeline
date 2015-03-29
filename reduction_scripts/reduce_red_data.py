#! /usr/bin/env python

import sys
import os
import pickle
import pyfits
import pywifes
import gc

#--- required for Fred's updates -------
import datetime

#---------------------- Fred's update -----------------------------------
start_time = datetime.datetime.now()

#------------------------------------------------------------------------
# get name of metadata file
meta_fn = sys.argv[1]
f1 = open(meta_fn, 'r')
obs_metadata = pickle.load(f1)
f1.close()

# testing on Fred's data
proj_dir = '/priv/maipenrai/skymap/mjc/wifes/20130622/'
data_dir = proj_dir+'raw_data/'
out_dir  = proj_dir+'proc_data/'

calib_prefix = out_dir+'wifesR_20130622'
my_data_hdu=0

# SET MULTITHREAD
#multithread=False
multithread=True

# SET SKIP ALREADY DONE FILES
skip_done=False
#skip_done=True

#------------------------------------------------------------------------
#------------------------------------------------------------------------
#************************************************************************
#*****                USER REDUCTION DESIGN IS SET HERE             *****
#************************************************************************
#  sky subtraction, no image combination
proc_steps = [
    #------------------
    {'step':'overscan_sub'   , 'run':False, 'suffix':'00', 'args':{}},
    {'step':'bpm_repair'     , 'run':False, 'suffix':'01', 'args':{}},
    #------------------
    {'step':'superbias'      , 'run':False, 'suffix':None,
     'args':{'method':'row_med', 'plot':True, 
             'verbose':False}},
    {'step':'bias_sub'       , 'run':False, 'suffix':'02',
     'args':{'method':'subtract', 'plot':False, 
             'verbose':False}},
    #------------------
    {'step':'superflat'      , 'run':False, 'suffix':None,
     'args':{'source':'dome'}},
    {'step':'superflat'      , 'run':False, 'suffix':None,
     'args':{'source':'twi', 'scale':'median_nonzero'}},
    {'step':'slitlet_profile', 'run':False, 'suffix':None, 'args':{}},
    #------------------
    {'step':'flat_cleanup'   , 'run':False, 'suffix':None,
     'args':{'type':['dome','twi'], 
             'verbose':True, 
             'plot':True,
             'buffer':4,
             'offsets':[0.4,0.4],
             'radius':10.0,
             'nsig_lim':3.0}},
    #------------------
    {'step':'superflat_mef'  , 'run':False, 'suffix':None,
     'args':{'source':'dome'}},
    {'step':'superflat_mef'  , 'run':False, 'suffix':None,
     'args':{'source':'twi'}},
    #------------------
    {'step':'slitlet_mef'    , 'run':False, 'suffix':'03',
     'args':{'ns':False}},
    #------------------
    {'step':'wave_soln'      , 'run':False, 'suffix':None,
     'args':{'verbose':True,
             'method' : 'optical',
             'doalphapfit' : True,
             'doplot' : True,
             'shift_method' : 'xcorr_all',
             'find_method' : 'mpfit',
             'dlam_cut_start':5.0}},
    {'step':'wire_soln'      , 'run':False, 'suffix':None, 'args':{}},
    {'step':'flat_response'  , 'run':False, 'suffix':None,
     'args':{'mode':'all'}},
    #------------------
    {'step':'cosmic_rays'    , 'run':False, 'suffix':'04',
     'args':{'ns':False, 'multithread':multithread}},
    #------------------
    {'step':'sky_sub'        , 'run':False, 'suffix':'05',
     'args':{'ns':False}},
    #------------------
    {'step':'obs_coadd'      , 'run':False, 'suffix':'06',
     'args':{'method':'sum'}},
    #------------------
    {'step':'flatfield'      , 'run':False, 'suffix':'07', 'args':{}},
    #------------------
    {'step':'cube_gen'       , 'run':False, 'suffix':'08',
     'args':{'multithread':multithread,'adr':True,
             'wmin_set':5400.0}},
    #------------------
    {'step':'extract_stars'  , 'run':False, 'suffix':None,
     'args':{'ytrim':4, 'type':'flux'}},
    {'step':'derive_calib'   , 'run':False, 'suffix':None,
     'args':{'plot_stars':False,
             'plot_sensf':True,
             'polydeg':10,
             'method':'poly',
             'norm_stars':True}},
    {'step':'flux_calib'     , 'run':False, 'suffix':'09', 'args':{}},
    #------------------
    {'step':'extract_stars'  , 'run':False, 'suffix':None,
     'args':{'ytrim':4, 'type':'telluric'}},
    {'step':'derive_telluric', 'run':False, 'suffix':None,
     'args':{'plot':True}},
    {'step':'telluric_corr'  , 'run':False, 'suffix':'10', 'args':{}},
    # Fred's update -> save the cube in a standard shape
    #------------------
    {'step':'save_3dcube'    , 'run':False, 'suffix':'11', 'args':{}}
    #------------------
    ]

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# METADATA WRANGLING FUNCTION
def get_full_obs_list(metadata):
    full_obs_list = []
    base_fn_list = (
        metadata['bias']+
        metadata['arc']+
        metadata['wire']+
        #metadata['dark']+
        metadata['domeflat']+
        metadata['twiflat'])
    for fn in base_fn_list:
        if fn not in full_obs_list:
            full_obs_list.append(fn)
    for obs in (metadata['sci']+
                metadata['std']):
        for key in obs.keys():
            # Fred's update 2 (another consequence of it ...)
            if key != 'type' : 
                for fn in obs[key]:
                    if fn not in full_obs_list:
                        full_obs_list.append(fn)
    return full_obs_list

def get_sci_obs_list(metadata):
    sci_obs_list = []
    for obs in metadata['sci']:
        for fn in obs['sci']:
            if fn not in sci_obs_list:
                sci_obs_list.append(fn)
    return sci_obs_list

# ------------------- Fred's update 2 ----------------
def get_std_obs_list(metadata,type = 'all'):
    std_obs_list = []
    for obs in metadata['std']:
        for fn in obs['sci']:
            if fn not in std_obs_list and type == 'all' :
                std_obs_list.append(fn)
            if fn not in std_obs_list and (type in obs['type']) :
                std_obs_list.append(fn)
    return std_obs_list

def get_sky_obs_list(metadata):
    sky_obs_list = []
    # Fred's update 2 -> to also use the ones from std * !
    for obs in metadata['sci']+metadata['std']:
        if 'sky' not in obs.keys():
            continue
        for fn in obs['sky']:
            if fn not in sky_obs_list:
                sky_obs_list.append(fn)
    return sky_obs_list

# ---------------- Fred's update -------------------
def get_associated_calib(metadata, this_fn, type):
    for obs in (metadata['sci'] + metadata['std']) :
        if 'sky' in obs.keys():
            sky = obs['sky']
        else:
            sky = []
        for fn in (obs['sci']+sky) :
            if fn == this_fn:
                if type in obs.keys():
                    if obs[type] != '' :
                        return obs[type]
    return False
        
#------------------
# primary ones!
def get_primary_sci_obs_list(metadata):
    sci_obs_list = [obs['sci'][0]
                    for obs in metadata['sci']]
    return sci_obs_list

def get_primary_std_obs_list(metadata, type='all'):
    if type =='all' :
        std_obs_list = [obs['sci'][0]
                        for obs in metadata['std']]
    elif type == 'telluric' or type == 'flux' :
        std_obs_list = []
        for obs in metadata['std']:
            if obs['sci'][0] not in std_obs_list and (type in obs['type']) :
                std_obs_list.append(obs['sci'][0])
    else :
       print 'Standard star type not understood !'
       print 'I will crash now ...'
    return std_obs_list

#------------------------------------------------------------------------
# NAMES FOR MASTER CALIBRATION FILES!!!
superbias_fn     = '%s_superbias.fits' % calib_prefix
superbias_fit_fn = '%s_superbias_fit.fits' % calib_prefix
super_dflat_raw  = '%s_super_domeflat_raw.fits' % calib_prefix 
super_dflat_fn   = '%s_super_domeflat.fits' % calib_prefix
super_dflat_mef  = '%s_super_domeflat_mef.fits' % calib_prefix
super_tflat_raw  = '%s_super_twiflat_raw.fits' % calib_prefix 
super_tflat_fn   = '%s_super_twiflat.fits' % calib_prefix
super_tflat_mef  = '%s_super_twiflat_mef.fits' % calib_prefix
slitlet_def_fn   = '%s_slitlet_defs.pkl' % calib_prefix
wsol_out_fn      = '%s_wave_soln.fits' % calib_prefix
wire_out_fn      = '%s_wire_soln.fits' % calib_prefix
flat_resp_fn     = '%s_resp_mef.fits' % calib_prefix
calib_fn         = '%s_calib.pkl' % calib_prefix
tellcorr_fn      = '%s_tellcorr.pkl' % calib_prefix

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# DEFINE THE PROCESSING STEPS

#------------------------------------------------------
# Subtract overscan
def run_overscan_sub(metadata, prev_suffix, curr_suffix):
    full_obs_list = get_full_obs_list(metadata)
    # this is the only time I hard code that this step should happen first
    for fn in full_obs_list:
        in_fn = '%s%s.fits' % (data_dir, fn)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        if skip_done and os.path.isfile(out_fn):
            continue
        print 'Subtracting Overscan for %s' % in_fn.split('/')[-1]
        pywifes.subtract_overscan(in_fn, out_fn, data_hdu=my_data_hdu)
    return

#------------------------------------------------------
# repair bad pixels!
def run_bpm_repair(metadata, prev_suffix, curr_suffix):
    full_obs_list = get_full_obs_list(metadata)
    for fn in full_obs_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        if skip_done and os.path.isfile(out_fn):
            continue
        print 'Repairing red bad pixels for %s' % in_fn.split('/')[-1]
        pywifes.repair_red_bad_pix(in_fn, out_fn, data_hdu=my_data_hdu)
    return

#------------------------------------------------------
# Generate super-bias
def run_superbias(metadata, prev_suffix, curr_suffix,
                  method='row_med', **args):
    bias_list = [
        '%s%s.p%s.fits' % (out_dir, x, prev_suffix)
        for x in metadata['bias']]
    print 'Calculating Global Superbias'
    pywifes.imcombine(bias_list, superbias_fn, data_hdu=my_data_hdu)
    # decide what bias model you will actually subtract - could be just data
    if method == 'fit' or method == 'row_med':
        # Fit a smart surface to the bias or take the median
        # A bit experimental so far ... but you know what you are doing, right ?
        pywifes.generate_wifes_bias_fit(
            superbias_fn, superbias_fit_fn, 
            data_hdu=my_data_hdu, method=method, **args)
    else:
        pywifes.imcopy(superbias_fn, superbias_fit_fn)
    # generate local superbiases for any science frames
    sci_obs_list  = get_sci_obs_list(metadata)
    std_obs_list  = get_std_obs_list(metadata)
    for fn in (sci_obs_list + std_obs_list):
        # Find if there is an associated bias (or more ?)
        local_biases = get_associated_calib(metadata,fn, 'bias')
        if local_biases:
            local_bias_fn = get_associated_calib(metadata,fn,'bias')[0]
            print 'Calculating Local Superbias for %s' % local_bias_fn
            local_superbias = '%s%s.fits' % (out_dir, local_bias_fn+'.lsb')
            local_superbias_fit = '%s%s.fits' % (out_dir, local_bias_fn+'.lsb_fit')
            if os.path.isfile(local_superbias_fit):
                continue
            # step 1 - coadd biases
            local_biases_filename = [
                '%s%s.p%s.fits' % (out_dir, x, prev_suffix)
                for x in local_biases]
            pywifes.imcombine(local_biases_filename,local_superbias, data_hdu=my_data_hdu)
            # step 2 - generate fit!
            if method == 'fit' or method == 'row_med':
                pywifes.generate_wifes_bias_fit(local_superbias, 
                                                local_superbias_fit,
                                                data_hdu=my_data_hdu,
                                                method=method, **args)
            else:
                pywifes.imcopy(local_superbias, local_superbias_fit)
    return

#----------------------------------------------------
# Subtract bias
def run_bias_sub(metadata, prev_suffix, curr_suffix,
                 method='sub', **args):
    full_obs_list = get_full_obs_list(metadata)
    sci_obs_list  = get_sci_obs_list(metadata)
    std_obs_list  = get_std_obs_list(metadata)
    sky_obs_list  = get_sky_obs_list(metadata)
    for fn in full_obs_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        if skip_done and os.path.isfile(out_fn):
            continue
        # figure out which bias to subtract
        local_biases = get_associated_calib(metadata,fn, 'bias')
        if local_biases:
            local_bias_fn = get_associated_calib(metadata,fn,'bias')[0]
            local_superbias = '%s%s.fits' % (out_dir, local_bias_fn+'.lsb')
            bias_fit_fn = '%s%s.fits' % (out_dir, local_bias_fn+'.lsb_fit')
            bias_type = 'local'
        else:
            bias_fit_fn = superbias_fit_fn
            bias_type = 'global'
        # subtract it!
        print 'Subtracting %s superbias for %s' % (
            bias_type, in_fn.split('/')[-1])
        if method == 'copy':
            pywifes.imcopy(in_fn, out_fn)
        else:
            pywifes.imarith(in_fn, '-', bias_fit_fn, out_fn, 
                            data_hdu=my_data_hdu)
    return

#------------------------------------------------------
# Generate super-flat
def run_superflat(metadata, prev_suffix, curr_suffix,
                  source, scale=None, method='median'):
    if source == 'dome':
        flat_list = [
            '%s%s.p%s.fits' % (out_dir, x, prev_suffix)
            for x in metadata['domeflat']]
        out_fn = super_dflat_raw
    elif source == 'twi':
        flat_list = [
            '%s%s.p%s.fits' % (out_dir, x, prev_suffix)
            for x in metadata['twiflat']]
        out_fn = super_tflat_raw
    else:
        raise ValueError, 'Flatfield type not recognized'
    print 'Generating co-add %sflat' % source
    pywifes.imcombine(flat_list, out_fn,
                      data_hdu=my_data_hdu,
                      scale=scale,
                      method=method)
    return

# Fred's flat cleanup
def run_flat_cleanup(metadata, prev_suffix, curr_suffix,
                   type=['dome','twi'],offsets=[0.,0.], **args):
    # check the slitlet definition file
    if os.path.isfile(slitlet_def_fn):
        slitlet_fn = slitlet_def_fn
    else:
        slitlet_fn=None
    if 'dome' in type :
        print 'Correcting master domeflat',super_dflat_fn.split('/')[-1]
        pywifes.interslice_cleanup(super_dflat_raw,super_dflat_fn, slitlet_fn,
                                   offset=offsets[type.index('dome')],
                                   method='2D',**args)
    if 'twi' in type :
        print 'Correcting master twilight flat',super_tflat_fn.split('/')[-1]
        pywifes.interslice_cleanup(super_tflat_raw,super_tflat_fn, slitlet_fn,
                                   offset=offsets[type.index('twi')],
                                   method='2D',**args)
    return

#------------------------------------------------------
# Fit slitlet profiles
def run_slitlet_profile(metadata, prev_suffix, curr_suffix, **args):
    if os.path.isfile(super_dflat_fn):
        flatfield_fn  = super_dflat_fn
    else:
        flatfield_fn  = super_dflat_raw
    output_fn = slitlet_def_fn
    pywifes.derive_slitlet_profiles(flatfield_fn,
                                    output_fn,
                                    data_hdu=my_data_hdu,
                                    **args)
    return

#------------------------------------------------------
# Create MEF files
def run_superflat_mef(metadata, prev_suffix, curr_suffix, source):
    if source == 'dome':
        if os.path.isfile(super_dflat_fn):
            in_fn  = super_dflat_fn
        else:
            in_fn = super_dflat_raw
        out_fn = super_dflat_mef
    elif source == 'twi':
        if os.path.isfile(super_tflat_fn):
            in_fn  = super_tflat_fn
        else :
            in_fn = super_tflat_raw
        out_fn = super_tflat_mef
    else:
        raise ValueError, 'Flatfield type not recognized'
    # check the slitlet definition file
    if os.path.isfile(slitlet_def_fn):
        slitlet_fn = slitlet_def_fn
    else:
        slitlet_fn=None
    # run it!
    print 'Generating MEF %sflat' % source
    pywifes.wifes_slitlet_mef(in_fn, out_fn,
                              data_hdu=my_data_hdu,
                              slitlet_def_file=slitlet_fn)
    return

def run_slitlet_mef(metadata, prev_suffix, curr_suffix, ns=False):
    full_obs_list = get_full_obs_list(metadata)
    sci_obs_list  = get_sci_obs_list(metadata)
    std_obs_list  = get_std_obs_list(metadata)
    sky_obs_list  = get_sky_obs_list(metadata)
    ns_proc_list = sci_obs_list+std_obs_list
    # check the slitlet definition file
    if os.path.isfile(slitlet_def_fn):
        slitlet_fn = slitlet_def_fn
    else:
        slitlet_fn=None
    for fn in full_obs_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        if skip_done and os.path.isfile(out_fn):
            continue
        print 'Creating MEF file for %s' % in_fn.split('/')[-1]
        if ns and fn in ns_proc_list:
            sky_fn = '%s%s.s%s.fits' % (out_dir, fn, curr_suffix)
            pywifes.wifes_slitlet_mef_ns(in_fn, out_fn, sky_fn,
                                         data_hdu=my_data_hdu,
                                         slitlet_def_file=slitlet_fn)
        else:
            pywifes.wifes_slitlet_mef(in_fn, out_fn, data_hdu=my_data_hdu,
                                      slitlet_def_file=slitlet_fn)
        gc.collect()
    return

#----------------- Fred's update -------------------
# Wavelength solution
def run_wave_soln(metadata, prev_suffix, curr_suffix, **args):
    # First, generate the master arc solution, based on generic arcs
    wsol_in_fn  = '%s%s.p%s.fits' % (out_dir, metadata['arc'][0],
                                     prev_suffix)
    print 'Deriving master wavelength solution from %s' % wsol_in_fn.split('/')[-1]
    pywifes.derive_wifes_wave_solution(wsol_in_fn, wsol_out_fn,
                                       **args)
    # local wave solutions for science or standards
    sci_obs_list  = get_sci_obs_list(metadata)
    std_obs_list  = get_std_obs_list(metadata)
    for fn in sci_obs_list + std_obs_list:
        # Check if the file has a dedicated arc associated with it ...
        # Only for Science and Std stars for now (sky not required at this stage)
        # (less critical for the rest anyway ...)
        #Mike I addition - reduce all local arcs (rather than just the first)
        local_arcs = get_associated_calib(metadata,fn, 'arc')   
        for local_arc in local_arcs :
            local_arc_fn = '%s%s.p%s.fits' % (out_dir, local_arc, prev_suffix)
            local_wsol_out_fn = '%s%s.wsol.fits' % (out_dir, local_arc)
            if os.path.isfile(local_wsol_out_fn):
                continue
            print 'Deriving local wavelength solution for %s' % local_arc
            pywifes.derive_wifes_wave_solution(local_arc_fn, local_wsol_out_fn,**args)
    return

#------------------------------------------------------
# Wire solution
def run_wire_soln(metadata, prev_suffix, curr_suffix):
    # Global wire solution
    wire_in_fn  = '%s%s.p%s.fits' % (out_dir,
                                     metadata['wire'][0],
                                     prev_suffix)
    print 'Deriving global wire solution from %s' % wire_in_fn.split('/')[-1]
    pywifes.derive_wifes_wire_solution(wire_in_fn, wire_out_fn)
    # Wire solutions for any specific obsevations
    sci_obs_list  = get_sci_obs_list(metadata)
    std_obs_list  = get_std_obs_list(metadata)
    for fn in sci_obs_list + std_obs_list:
        # Check if the file has a dedicated wire associated with it ...
        # Only for Science and Std stars for now (sky not required at this stage)
        # (less critical for the rest anyway ...)
        local_wires = get_associated_calib(metadata,fn, 'wire')            
        if local_wires :
            local_wire_fn = '%s%s.p%s.fits' % (out_dir, local_wires[0], prev_suffix)
            local_wire_out_fn = '%s%s.wire.fits' % (out_dir, local_wires[0])
            if os.path.isfile(local_wire_out_fn):
                continue
            print 'Deriving local wire solution for %s' % local_wires[0]
            pywifes.derive_wifes_wire_solution(local_wire_fn, local_wire_out_fn)
    return

#------------------------------------------------------
# Cosmic Rays
def run_cosmic_rays(metadata, prev_suffix, curr_suffix,
                    ns=False, multithread=False):
    from lacosmic import lacos_wifes
    # now run ONLY ON SCIENCE TARGETS AND STANDARDS
    sci_obs_list  = get_sci_obs_list(metadata)
    sky_obs_list  = get_sky_obs_list(metadata)
    std_obs_list  = get_std_obs_list(metadata)
    for fn in sci_obs_list+sky_obs_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        print 'Cleaning cosmics in %s' % in_fn.split('/')[-1]
        # skip files which are already done
        #if os.path.isfile(out_fn):
        #    continue
        lacos_wifes(in_fn, out_fn, wsol_fn=wsol_out_fn, niter=3,
                    sig_clip=10.0, obj_lim=10.0, sig_frac=0.2,
                    multithread=multithread)
        if ns:
            in_fn  = '%s%s.s%s.fits' % (out_dir, fn, prev_suffix)
            out_fn = '%s%s.s%s.fits' % (out_dir, fn, curr_suffix)
            print 'Cleaning cosmics in %s' % in_fn.split('/')[-1]
            lacos_wifes(in_fn, out_fn, wsol_fn=wsol_out_fn, niter=3,
                        sig_clip=10.0, obj_lim=10.0, sig_frac=0.2,
                        multithread=multithread)
        gc.collect()
    for fn in std_obs_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        print 'Cleaning cosmics in %s' % in_fn.split('/')[-1]
        #lacos_wifes(in_fn, out_fn, niter=1, sig_frac=2.0)
        lacos_wifes(in_fn, out_fn, wsol_fn=wsol_out_fn, niter=3,
                    sig_clip=10.0, obj_lim=10.0, sig_frac=0.2,
                    multithread=multithread)
        if ns:
            in_fn  = '%s%s.s%s.fits' % (out_dir, fn, prev_suffix)
            out_fn = '%s%s.s%s.fits' % (out_dir, fn, curr_suffix)
            print 'Cleaning cosmics in %s' % in_fn.split('/')[-1]
            #lacos_wifes(in_fn, out_fn, niter=1, sig_frac=2.0)
            lacos_wifes(in_fn, out_fn, wsol_fn=wsol_out_fn, niter=3,
                        sig_clip=10.0, obj_lim=10.0, sig_frac=0.2,
                        multithread=multithread)
        gc.collect()
    return

#------------------------------------------------------
# Sky subtraction!
def run_sky_sub_ns(metadata, prev_suffix, curr_suffix):
    sci_obs_list  = get_sci_obs_list(metadata)
    std_obs_list  = get_std_obs_list(metadata)
    ns_proc_list = sci_obs_list+std_obs_list
    for fn in ns_proc_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        sky_fn = '%s%s.s%s.fits' % (out_dir, fn, prev_suffix)
        print 'Subtracting N+S sky frame for %s' % in_fn.split('/')[-1]
        pywifes.scaled_imarith_mef(
            in_fn, '-', sky_fn, out_fn,
            scale='exptime')
    return

def run_sky_sub(metadata, prev_suffix, curr_suffix, ns=False):
    if ns:
        run_sky_sub_ns(metadata, prev_suffix, curr_suffix)
    else:
        # subtract sky frames from science objects
        for obs in metadata['sci']:
            if len(obs['sky']) > 0:
                sky_fn = obs['sky'][0]
                sky_proc_fn = '%s%s.p%s.fits' % (out_dir, sky_fn, prev_suffix)
                for fn in obs['sci']:
                    in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
                    out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
                    print 'Subtracting sky frame for %s' % in_fn.split('/')[-1]
                    # subtract scaled sky frame!
                    pywifes.scaled_imarith_mef(
                        in_fn, '-', sky_proc_fn, out_fn,
                        scale='exptime')
            else:
                for fn in obs['sci']:
                    in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
                    out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
                    print 'Copying image %s' % in_fn.split('/')[-1]
                    # subtract scaled sky frame!
                    pywifes.imcopy(in_fn, out_fn)
        # copy stdstar frames
        std_obs_list = get_std_obs_list(metadata)
        for fn in std_obs_list:
            in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
            out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
            print 'Copying standard star image %s' % in_fn.split('/')[-1]
            pywifes.imcopy(in_fn, out_fn)
    return

#------------------------------------------------------
# Image coaddition for science and standards!
def run_obs_coadd(metadata, prev_suffix, curr_suffix,
                  method='sum', scale=None):
    for obs in (metadata['sci']+metadata['std']):
        # if just one, then copy it
        if len(obs['sci']) == 1:
            fn = obs['sci'][0]
            in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
            out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
            print 'Copying image %s' % in_fn.split('/')[-1]
            pywifes.imcopy(in_fn, out_fn)            
        # coadd sci frames!
        else:
            in_fn_list = ['%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
                          for fn in obs['sci']]
            out_fn = '%s%s.p%s.fits' % (out_dir, obs['sci'][0], curr_suffix)
            print 'Coadding images for %s' % in_fn_list[0].split('/')[-1]
            pywifes.imcombine_mef(in_fn_list, out_fn,
                                  scale=scale,
                                  method=method)
    return

#------------------------------------------------------
# Flatfield Response
def run_flat_response(metadata, prev_suffix, curr_suffix,
                      mode='all'):
    # now fit the desired style of response function
    print 'Generating flatfield response function'
    if mode == 'all':
        pywifes.wifes_2dim_response(super_dflat_mef,
                                    super_tflat_mef,
                                    flat_resp_fn,
                                    wsol_fn=wsol_out_fn)
    elif mode == 'dome':
        pywifes.wifes_response_poly(super_dflat_mef,
                                    flat_resp_fn,
                                    wsol_fn=wsol_out_fn)
    else:
        raise ValueError, 'Requested response mode not recognized'
    return

#------------------------------------------------------
# Flatfield Division
def run_flatfield(metadata, prev_suffix, curr_suffix):
    sci_obs_list = get_primary_sci_obs_list(metadata)
    std_obs_list = get_primary_std_obs_list(metadata)
    for fn in sci_obs_list+std_obs_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        print 'Flat-fielding image %s' % in_fn.split('/')[-1]
        pywifes.imarith_mef(in_fn, '/',
                            flat_resp_fn,
                            out_fn)
    return

#------------------------------------------------------
# Data Cube Generation
def run_cube_gen(metadata, prev_suffix, curr_suffix, **args):
    # now generate cubes
    sci_obs_list = get_primary_sci_obs_list(metadata)
    std_obs_list = get_primary_std_obs_list(metadata)
    for fn in sci_obs_list+std_obs_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        print 'Generating Data Cube for %s' % in_fn.split('/')[-1]
        # decide whether to use global or local wsol and wire files
        local_wires = get_associated_calib(metadata,fn, 'wire')
        if local_wires :
            wire_fn = '%s%s.wire.fits' % (out_dir, local_wires[0])
            print '(Note: using %s as wire file)' % wire_fn.split('/')[-1]
        else:
            wire_fn = wire_out_fn
        local_arcs = get_associated_calib(metadata,fn, 'arc')
        #Mike I addition: enable linear interpolation between wavelength
        #scales found from a before and after arc file.
        if local_arcs :
            if len(local_arcs)==1:
                wsol_fn = '%s%s.wsol.fits' % (out_dir, local_arcs[0])
                print '(Note: using %s as wsol file)' % wsol_fn.split('/')[-1]
            else:
                #Assume that the two arc files are before and after.
                #(WARNING: Probably not a great assumption - the code
                #behaviour could be different with e.g. two low SNR arcs
                #both before the science data)
                fn0 = '%s%s.wsol.fits' % (out_dir, local_arcs[0])
                fn1 = '%s%s.wsol.fits' % (out_dir, local_arcs[1])
                fits0 = pyfits.open(fn0)
                fits1 = pyfits.open(fn1)
                #Compute the weights w from the filenames, which have the times in them.
                #Unfortunately, there are no times in the headers. That would be a better
                #way to do this though!
                ix = fn0.rfind('-')
                t0 = float(fn0[ix-6:ix-4]) + \
                     float(fn0[ix-4:ix-2])/60.0 + \
                     float(fn0[ix-2:ix])/60.0/60.0
                ix = fn1.rfind('-')
                t1 = float(fn1[ix-6:ix-4]) + \
                     float(fn1[ix-4:ix-2])/60.0 + \
                     float(fn1[ix-2:ix])/60.0/60.0
                ix = fn.rfind('-')
                ts = float(fn[ix-6:ix-4]) + \
                     float(fn[ix-4:ix-2])/60.0 + \
                     float(fn[ix-2:ix])/60.0/60.0
                w0 = (t1 - ts)/(t1 - t0)
                w1 = (ts - t0)/(t1 - t0)
                #Now find a weighted average of the wavelength solution slices
                #and save to a file named according to the target name.
                n_slices = len(fits0) 
                for i in range(1,n_slices):
                    fits0[i].data = w0*fits0[i].data + w1*fits1[i].data
                wsol_fn = '%s%s.wsol.fits' % (out_dir, fn) 
                fits0.writeto(wsol_fn)
        else:
            wsol_fn = wsol_out_fn
        #Now that we've sorted the file names out, generate the cube!
        pywifes.generate_wifes_cube(
            in_fn, out_fn,
            wire_fn=wire_fn,
            wsol_fn=wsol_fn,
            ny_orig=76, offset_orig=2.0, **args)
        #print squirrel
    return

#------------------------------------------------------
# Standard star extraction
def run_extract_stars(metadata, prev_suffix, curr_suffix, type='all',**args):
    # for each std, extract spectrum as desired
    std_obs_list = get_primary_std_obs_list(metadata, type=type)
    #print std_obs_list
    for fn in std_obs_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.x%s.dat'  % (out_dir, fn, prev_suffix)
        print 'Extract %s standard star from %s' % (type, in_fn.split('/')[-1])
        pywifes.extract_wifes_stdstar(in_fn,
                                      save_fn=out_fn,
                                      save_mode='ascii',
                                      **args)
    return

# Sensitivity Function fit
def run_derive_calib(metadata, prev_suffix, curr_suffix, method = 'poly',**args):
    std_obs_list = get_primary_std_obs_list(metadata, type='flux')
    std_cube_list = ['%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
                     for fn in std_obs_list]
    extract_list = ['%s%s.x%s.dat' % (out_dir, fn, prev_suffix)
                    for fn in std_obs_list]
    print 'Deriving sensitivity function'
    best_calib = pywifes.derive_wifes_calibration(
        std_cube_list,
        calib_fn,
        extract_in_list=extract_list,method=method,
        **args)
    return

# Applying Calibration
def run_flux_calib(metadata, prev_suffix, curr_suffix,
                   mode='pywifes', **args):
    # calibrate all sci and std obs
    sci_obs_list = get_primary_sci_obs_list(metadata)
    std_obs_list = get_primary_std_obs_list(metadata)
    for fn in sci_obs_list+std_obs_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        print 'Flux-calibrating cube %s' % in_fn.split('/')[-1]
        pywifes.calibrate_wifes_cube(
            in_fn, out_fn, calib_fn, mode)
    return

#------------------------------------------------------
# Telluric - derive
def run_derive_telluric(metadata, prev_suffix, curr_suffix, **args):
    std_obs_list = get_primary_std_obs_list(metadata, 'telluric')
    std_cube_list = ['%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
                     for fn in std_obs_list]
    extract_list = ['%s%s.x%s.dat' % (out_dir, fn, prev_suffix)
                    for fn in std_obs_list]
    print 'Deriving telluric correction'
    pywifes.derive_wifes_telluric(std_cube_list,
                                  tellcorr_fn,
                                  extract_in_list=extract_list,
                                  **args)
    return

# Telluric - apply
def run_telluric_corr(metadata, prev_suffix, curr_suffix, **args):
    # calibrate all sci and std obs
    sci_obs_list = get_primary_sci_obs_list(metadata)
    std_obs_list = get_primary_std_obs_list(metadata)
    for fn in sci_obs_list+std_obs_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        print 'Correcting telluric in %s' % in_fn.split('/')[-1]
        pywifes.apply_wifes_telluric(
            in_fn, out_fn, tellcorr_fn)
    return

#------------------------- Fred's update -----------------
# Save final cube in suitable fits format
def run_save_3dcube(metadata, prev_suffix, curr_suffix, **args):
    # now generate cubes
    sci_obs_list = get_primary_sci_obs_list(metadata)
    for fn in sci_obs_list:
        in_fn  = '%s%s.p%s.fits' % (out_dir, fn, prev_suffix)
        out_fn = '%s%s.p%s.fits' % (out_dir, fn, curr_suffix)
        print 'Saving 3D Data Cube for %s' % in_fn.split('/')[-1]
        pywifes.generate_wifes_3dcube(
            in_fn, out_fn, **args)
    return

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# RUN THE PROCESSING STEPS
prev_suffix = None
for step in proc_steps:
    step_name   = step['step']
    step_run    = step['run']
    step_suffix = step['suffix']
    step_args   = step['args']
    func_name = 'run_'+step_name
    func = globals()[func_name]
    if step_run:
        func(obs_metadata,
             prev_suffix = prev_suffix,
             curr_suffix = step_suffix,
             **step_args)
    if step_suffix != None:
        prev_suffix = step_suffix

#------------------------------------------------------------------------
#------------------------------------------------------------------------


#------------------------- Fred's update --------------------------------
duration = datetime.datetime.now() - start_time
print 'All done in %.01f seconds.' % duration.total_seconds()
