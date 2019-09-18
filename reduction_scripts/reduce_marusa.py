#! /usr/bin/env python

"""
Example:
python reduce_red_data_python3.py config.py /priv/mulga1/marusa/2m3data/20190302
"""

import sys
import os
import pickle
from astropy.io import fits as pyfits
import pywifes
import gc
import datetime
import warnings
import numpy as np
import imp
#~ from importlib import reload

#~ import calibration_filenames_date as cal
#~ cal=cal.result

#------------------------------------------------------------------------
start_time = datetime.datetime.now()
#------------------------------------------------------------------------

# Config
config_name=sys.argv[1]
config = imp.load_source(config_name.replace('.py', ''), config_name)
reload(config)
#~ metadata=config.reduce_metadata

# This is not so simple as output folder is not the same
# What is this?
#absolute_paths=True

prefix=config.prefix

obsdate = sys.argv[2]

# Input folder with raw data
data_dir = os.path.join(config.input_root, obsdate) #sys.argv[2]


# Output folder (should already exist and metadata should be there)

###### OUTPUT FOLDER ################
# Get obsdate
#~ path = sys.argv[2]
#~ obsdate = path.split('/')[-1]
#~ if len(obsdate)<2: # in case path ends with /
    #~ obsdate = path.split('/')[-2] # Hope that works
print('OBSDATE', obsdate)

root_obsdate = os.path.join(config.output_root, '%s'%obsdate)

# Create folder with date
root_bool = os.path.isdir(root_obsdate) and os.path.exists(root_obsdate)
print 'ROOT_BOOL', root_obsdate, os.path.isdir(root_obsdate), os.path.exists(root_obsdate)
print 'TEST', os.path.isdir('/data/mash/marusa/reduction_wifes/pipeline/reduction_scripts/'), os.path.exists('/data/mash/marusa/reduction_wifes/pipeline/reduction_scripts/')
if not root_bool:
    os.mkdir(root_obsdate)
print('root_obsdate', root_obsdate)

# Add band (grating)
out_dir = os.path.join(root_obsdate, 'reduced_%s'%config.band)

if prefix is not None and len(prefix)>0:
    print ('prefix')
    out_dir= out_dir + '_%s'%prefix
out_dir_bool = os.path.isdir(out_dir) and os.path.exists(out_dir)
if not out_dir_bool:
    os.mkdir(out_dir)
print('out_dir', out_dir)
######################################

#~ path = sys.argv[2]
#~ obsdate = path.split('/')[-1]
#~ if len(obsdate)<1:
    #~ obsdate = path.split('/')[-2]
#~ print('OBSDATE', obsdate)
#~ if prefix is not None and len(prefix)>0:
    #~ out_dir = os.path.join(config.output_root, '%s/reduced_r_%s'%(obsdate, prefix))
#~ else:
    #~ out_dir = os.path.join(config.output_root, '%s/reduced_r'%obsdate)
#~ print config.output_root, obsdate, prefix
#~ out_dir_bool = os.path.isdir(out_dir) and os.path.exists(out_dir)
#~ print('Trying to create', out_dir)
#~ if not out_dir_bool:
    #~ os.mkdir(out_dir)
    #~ print('DIR %s CREATED.'%out_dir)
#~ print('out_dir', out_dir)

if config.band == 'r':
    calib_prefix = os.path.join(out_dir,'wifesR')
elif config.band == 'b':
    calib_prefix = os.path.join(out_dir,'wifesB')

# Load metadata
if prefix is not None and len(prefix)>0:
    if config.band == 'r':
        metadata_filename=os.path.join(out_dir, '%s_metadata_WiFeSRed.py'%prefix)
    elif config.band == 'b':
        metadata_filename=os.path.join(out_dir, '%s_metadata_WiFeSBlue.py'%prefix)
else:
    if config.band == 'r':
        metadata_filename=os.path.join(out_dir, 'metadata_WiFeSRed.py')
    elif config.band == 'b':
        metadata_filename=os.path.join(out_dir, 'metadata_WiFeSBlue.py')
print('metadata_filename', metadata_filename, config.metadata_filename)
obs_metadata = imp.load_source('obs_metadata', metadata_filename).night_data
#~ mode = imp.load_source('obs_metadata', metadata_filename).mode

# Some WiFeS specific things
my_data_hdu=0

# SET MULTITHREAD ?
multithread=config.multithread

# SET SKIP ALREADY DONE FILES ?
skip_done=config.skip_done

proc_steps = config.proc_steps


def find_the_mode(): # delete this


    #~ selected_cal_dates={'DARK':20190304, 'ZERO':20190302, 'FLAT': 20190304}
    selected_cal_dates={'DARK':20190304}

    # Calibrations
    try:
        c=cal[mode]    
    except:
        c=None

    for kk, vv in selected_cal_dates.iteritems(): # For imagetype kk, dict of dates vv
        t=c[kk]
        d=t[vv]
        print kk
        for x in d:
            print 'copy', x, os.path.join(root, x.split('/')[-1])#, x, root
            copyfile(x, os.path.join(root, x.split('/')[-1]))
        print
    print
    print
    print




    if c:
        print '#----------------------------------------------------'
        print 'Available calibrations:'
        for kk, vv in c.iteritems():
            print kk
            for kkk, vvv in vv.iteritems():
                print kkk, len(vvv)
            print

    print
    print
    print '#########################################################'
    #~ print
    #~ print






#------------------------------------------------------------------------
#------------------------------------------------------------------------
# METADATA WRANGLING FUNCTION
def get_full_obs_list(metadata):
    full_obs_list = []
    base_fn_list = (
        metadata['bias']+
        metadata['arc']+
        metadata['wire']+
        #metadata['dark']+ # WHY??? This has been commented out originally.
        metadata['domeflat']+
        metadata['twiflat'])
    for fn in base_fn_list:
        if fn not in full_obs_list:
            full_obs_list.append(fn)
    for obs in (metadata['sci']+
                metadata['std']):
        for key in obs.keys():
            # Fred's update 2 (another consequence of it ...)
            if (key != 'type') and (key != 'name'): 
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
       print('Standard star type not understood !')
       print('I will crash now ...')
    return std_obs_list

def check_if_obsdates_of_data_dir_and_filename_are_the_same(data_dir, fn):
    obsdate_fn = fn.split('-')[1].split('.')[0]
    obsdate_data_dir = data_dir.split('/')[-1]
    
    if len(obsdate_data_dir)<1:
        obsdate_data_dir = data_dir.split('/')[-2] # There must be a '/' at the very end of data_dir...
        
    if obsdate_fn==obsdate_data_dir:
        return False # No correction needed
    else:
        data_dir=data_dir.replace(obsdate_data_dir, obsdate_fn)
        return data_dir
     

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
        
        # Check if we have the right folder (calibrations might be somewhere else)
        change_dir = check_if_obsdates_of_data_dir_and_filename_are_the_same(data_dir, fn)
        if change_dir:
            in_fn = os.path.join(change_dir, '%s.fits' %fn)
        else:
            in_fn = os.path.join(data_dir, '%s.fits' %fn)
            
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        if skip_done and os.path.isfile(out_fn):
            continue
        print('Subtracting Overscan for %s' % in_fn.split('/')[-1])
        pywifes.subtract_overscan(in_fn, out_fn, data_hdu=my_data_hdu)
    return

#------------------------------------------------------
# repair bad pixels!
def run_bpm_repair(metadata, prev_suffix, curr_suffix):
    full_obs_list = get_full_obs_list(metadata)
    for fn in full_obs_list:
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        if skip_done and os.path.isfile(out_fn):
            continue
        
        if config.band == 'r':
            print('Repairing red bad pixels for %s' % in_fn.split('/')[-1])
            pywifes.repair_red_bad_pix(in_fn, out_fn, data_hdu=my_data_hdu)
        elif config.band == 'b':
            print('Repairing blue bad pixels for %s' % in_fn.split('/')[-1])
            pywifes.repair_blue_bad_pix(in_fn, out_fn, data_hdu=my_data_hdu)
    return

#------------------------------------------------------
# Generate super-bias
def run_superbias(metadata, prev_suffix, curr_suffix,
                  method='row_med', **args):
    bias_list = [
        os.path.join(out_dir, '%s.p%s.fits' % (x, prev_suffix))
        for x in metadata['bias']]
    print('Calculating Global Superbias')
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
            print('Calculating Local Superbias for %s' % local_bias_fn)
            local_superbias = '%s%s.fits' % (out_dir, local_bias_fn+'.lsb')
            local_superbias_fit = os.path.join(out_dir, '%s.fits' % (local_bias_fn+'.lsb_fit'))
            if os.path.isfile(local_superbias_fit):
                continue
            # step 1 - coadd biases
            local_biases_filename = [
                os.path.join(out_dir, '%s.p%s.fits' % (x, prev_suffix))
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
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        if skip_done and os.path.isfile(out_fn):
            continue
        # figure out which bias to subtract
        local_biases = get_associated_calib(metadata,fn, 'bias')
        if local_biases:
            local_bias_fn = get_associated_calib(metadata,fn,'bias')[0]
            local_superbias = os.path.join(out_dir, '%s.fits' % (local_bias_fn+'.lsb'))
            bias_fit_fn = os.path.join(out_dir, '%s.fits' % (local_bias_fn+'.lsb_fit'))
            bias_type = 'local'
        else:
            bias_fit_fn = superbias_fit_fn
            bias_type = 'global'
        # subtract it!
        print('Subtracting %s superbias for %s' % (
            bias_type, in_fn.split('/')[-1]))
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
            os.path.join(out_dir, '%s.p%s.fits' % (x, prev_suffix))
            for x in metadata['domeflat']]
        out_fn = super_dflat_raw
    elif source == 'twi':
        flat_list = [
            os.path.join(out_dir, '%s.p%s.fits' % (x, prev_suffix))
            for x in metadata['twiflat']]
        out_fn = super_tflat_raw
    else:
        raise ValueError('Flatfield type not recognized')
    print('Generating co-add %sflat' % source)
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
        print('Correcting master domeflat',super_dflat_fn.split('/')[-1])
        pywifes.interslice_cleanup(super_dflat_raw, super_dflat_fn, slitlet_fn,
                                   offset=offsets[type.index('dome')],
                                   method='2D',**args)
    if 'twi' in type :
        print('Correcting master twilight flat',super_tflat_fn.split('/')[-1])
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
        raise ValueError('Flatfield type not recognized')
    # check the slitlet definition file
    if os.path.isfile(slitlet_def_fn):
        slitlet_fn = slitlet_def_fn
    else:
        slitlet_fn=None
    # run it!
    print('Generating MEF %sflat' % source)
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
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        if skip_done and os.path.isfile(out_fn):
            continue
        print('Creating MEF file for %s' % in_fn.split('/')[-1])
        if ns and fn in ns_proc_list:
            sky_fn = os.path.join(out_dir, '%s.s%s.fits' % (fn, curr_suffix))
            pywifes.wifes_slitlet_mef_ns(in_fn, out_fn, sky_fn,
                                         data_hdu=my_data_hdu,
                                         slitlet_def_file=slitlet_fn)
        else:
            pywifes.wifes_slitlet_mef(in_fn, out_fn, data_hdu=my_data_hdu,
                                      slitlet_def_file=slitlet_fn)
        gc.collect()
    return

#------------------------------------------------------
# Wavelength solution
def run_wave_soln(metadata, prev_suffix, curr_suffix, **args):
    # First, generate the master arc solution, based on generic arcs
    wsol_in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (metadata['arc'][0],
                                     prev_suffix))
    print('Deriving master wavelength solution from %s' % wsol_in_fn.split('/')[-1])
    pywifes.derive_wifes_wave_solution(wsol_in_fn, wsol_out_fn,
                                       **args)
    # local wave solutions for science or standards
    sci_obs_list  = get_sci_obs_list(metadata)
    std_obs_list  = get_std_obs_list(metadata)
    for fn in sci_obs_list + std_obs_list:
        # Check if the file has a dedicated arc associated with it ...
        # Only for Science and Std stars for now (sky not required at this stage)
        # (less critical for the rest anyway ...)
        # As per Mike I. pull request: if two arcs are present, find a solution
        # for both to later interpolate between them.
        # Restrict it to the first two arcs in the list (in case the feature is
        # being unknowingly used).
        local_arcs = get_associated_calib(metadata,fn, 'arc')            
        if local_arcs :
            for i in range(np.min([2,np.size(local_arcs)])):
                local_arc_fn = os.path.join(out_dir, '%s.p%s.fits' % (local_arcs[i], prev_suffix))
                local_wsol_out_fn = os.path.join(out_dir, '%s.wsol.fits' % (local_arcs[i]))
                if os.path.isfile(local_wsol_out_fn):
                    continue
                print('Deriving local wavelength solution for %s' % local_arcs[i])
                pywifes.derive_wifes_wave_solution(local_arc_fn, local_wsol_out_fn,
                                                    **args)
    return

#------------------------------------------------------
# Wire solution
def run_wire_soln(metadata, prev_suffix, curr_suffix):
    # Global wire solution
    wire_in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (
                                     metadata['wire'][0],
                                     prev_suffix))
    print('Deriving global wire solution from %s' % wire_in_fn.split('/')[-1])
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
            local_wire_fn = os.path.join(out_dir, '%s.p%s.fits' % (local_wires[0], prev_suffix))
            local_wire_out_fn = os.path.join(out_dir, '%s.wire.fits' % (local_wires[0]))
            if os.path.isfile(local_wire_out_fn):
                continue
            print('Deriving local wire solution for %s' % local_wires[0])
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
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        print('Cleaning cosmics in %s' % in_fn.split('/')[-1])
        # skip files which are already done
        #if os.path.isfile(out_fn):
        #    continue
        lacos_wifes(in_fn, out_fn, wsol_fn=wsol_out_fn, niter=3,
                    sig_clip=10.0, obj_lim=10.0, sig_frac=0.2,
                    multithread=multithread)
        if ns:
            in_fn  = os.path.join(out_dir, '%s.s%s.fits' % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, '%s.s%s.fits' % (fn, curr_suffix))
            print('Cleaning cosmics in %s' % in_fn.split('/')[-1])
            lacos_wifes(in_fn, out_fn, wsol_fn=wsol_out_fn, niter=3,
                        sig_clip=10.0, obj_lim=10.0, sig_frac=0.2,
                        multithread=multithread)
        gc.collect()
    for fn in std_obs_list:
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        print('Cleaning cosmics in %s' % in_fn.split('/')[-1])
        #lacos_wifes(in_fn, out_fn, niter=1, sig_frac=2.0)
        lacos_wifes(in_fn, out_fn, wsol_fn=wsol_out_fn, niter=3,
                    sig_clip=10.0, obj_lim=10.0, sig_frac=0.2,
                    multithread=multithread)
        if ns:
            in_fn  = os.path.join(out_dir, '%s.s%s.fits' % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, '%s.s%s.fits' % (fn, curr_suffix))
            print('Cleaning cosmics in %s' % in_fn.split('/')[-1])
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
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        sky_fn = os.path.join(out_dir, '%s.s%s.fits' % (fn, prev_suffix))
        print('Subtracting N+S sky frame for %s' % in_fn.split('/')[-1])
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
                sky_proc_fn = os.path.join(out_dir, '%s.p%s.fits' % (sky_fn, prev_suffix))
                for fn in obs['sci']:
                    in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
                    out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
                    print('Subtracting sky frame for %s' % in_fn.split('/')[-1])
                    # subtract scaled sky frame!
                    pywifes.scaled_imarith_mef(
                        in_fn, '-', sky_proc_fn, out_fn,
                        scale='exptime')
            else:
                for fn in obs['sci']:
                    in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
                    out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
                    print('Copying image %s' % in_fn.split('/')[-1])
                    # subtract scaled sky frame!
                    pywifes.imcopy(in_fn, out_fn)
        # copy stdstar frames
        std_obs_list = get_std_obs_list(metadata)
        for fn in std_obs_list:
            in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
            print('Copying standard star image %s' % in_fn.split('/')[-1])
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
            in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
            print('Copying image %s' % in_fn.split('/')[-1])
            pywifes.imcopy(in_fn, out_fn)            
        # coadd sci frames!
        else:
            in_fn_list = [os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
                          for fn in obs['sci']]
            out_fn = os.path.join(out_dir, '%s.p%s.fits' % (obs['sci'][0], curr_suffix))
            print('Coadding images for %s' % in_fn_list[0].split('/')[-1])
            pywifes.imcombine_mef(in_fn_list, out_fn,
                                  scale=scale,
                                  method=method)
    return

#------------------------------------------------------
# Flatfield Response
def run_flat_response(metadata, prev_suffix, curr_suffix,
                      mode='all'):
    # now fit the desired style of response function
    print('Generating flatfield response function')
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
        raise ValueError('Requested response mode not recognized')
    return

#------------------------------------------------------
# Flatfield Division
def run_flatfield(metadata, prev_suffix, curr_suffix):
    sci_obs_list = get_primary_sci_obs_list(metadata)
    #~ sci_obs_list = get_sci_obs_list(metadata) # MZ: TESTING. Something still wrong here. Use this line if coadd=False
    std_obs_list = get_primary_std_obs_list(metadata)
    for fn in sci_obs_list+std_obs_list:
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        print('Flat-fielding image %s' % in_fn.split('/')[-1])
        pywifes.imarith_mef(in_fn, '/',
                            flat_resp_fn,
                            out_fn)
    return

#------------------------------------------------------
# Data Cube Generation
def run_cube_gen(metadata, prev_suffix, curr_suffix, **args):
    # now generate cubes
    sci_obs_list = get_primary_sci_obs_list(metadata)
    #~ sci_obs_list = get_sci_obs_list(metadata) # MZ: TESTING... Something still wrong here. Use this line if coadd=False
    std_obs_list = get_primary_std_obs_list(metadata)
    for fn in sci_obs_list+std_obs_list:
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        print('Generating Data Cube for %s' % in_fn.split('/')[-1])
        # decide whether to use global or local wsol and wire files
        local_wires = get_associated_calib(metadata,fn, 'wire')
        if local_wires :
            wire_fn = os.path.join(out_dir, '%s.wire.fits' % (local_wires[0]))
            print('(Note: using %s as wire file)' % wire_fn.split('/')[-1])
        else:
            wire_fn = wire_out_fn
        local_arcs = get_associated_calib(metadata,fn, 'arc')
        if local_arcs :
            # Do I have two arcs ? Do they surround the Science file ?
            # Implement linear interpolation as suggested by Mike I.
            if len(local_arcs) ==2:
                # First, get the Science time
                f = pyfits.open(in_fn)
                sci_header = f[0].header
                sci_time = sci_header['DATE-OBS']
                # Now get the arc times
                arc_times = ['','']
                for i in range(2):
                    # Fetch the arc time from the "extra" pkl file
                    local_wsol_out_fn_extra = os.path.join(out_dir, '%s.wsol.fits_extra.pkl' % (local_arcs[i]) )
                    f = open(local_wsol_out_fn_extra, 'rb')
                    try:
                        f_pickled = pickle.load(f, protocol=2) # TODO: see if this works
                    except:
                        f_pickled = pickle.load(f) # Python 3
                    f.close()
                    arc_times[i] = f_pickled[-1][0]
                 
                # Now, make sure the Science is between the arcs:
                t0 = datetime.datetime( np.int(arc_times[0].split('-')[0]),
                                        np.int(arc_times[0].split('-')[1]),
                                        np.int(arc_times[0].split('-')[2].split('T')[0]),
                                        np.int(arc_times[0].split('T')[1].split(':')[0]),
                                        np.int(arc_times[0].split(':')[1]),
                                        np.int(arc_times[0].split(':')[2].split('.')[0]),
                                        )
                t1 = datetime.datetime( np.int(sci_time.split('-')[0]),
                                        np.int(sci_time.split('-')[1]),
                                        np.int(sci_time.split('-')[2].split('T')[0]),
                                        np.int(sci_time.split('T')[1].split(':')[0]),
                                        np.int(sci_time.split(':')[1]),
                                        np.int(sci_time.split(':')[2].split('.')[0]),
                                        )
                t2 = datetime.datetime( np.int(arc_times[1].split('-')[0]),
                                        np.int(arc_times[1].split('-')[1]),
                                        np.int(arc_times[1].split('-')[2].split('T')[0]),
                                        np.int(arc_times[1].split('T')[1].split(':')[0]),
                                        np.int(arc_times[1].split(':')[1]),
                                        np.int(arc_times[1].split(':')[2].split('.')[0]),
                                        )
                ds1 = (t1 - t0).total_seconds()
                ds2 = (t2 - t1).total_seconds()
                if ds1>0 and ds2>0:
                    # Alright, I need to interpolate betweent the two arcs
                    # There is a difference between blue and red arm: # MZ
                    if config.band == 'r':
                        w1 = ds2/(ds1+ds2)
                        w2 = ds1/(ds1+ds2)
                    elif config.band == 'b':
                        w1 = ds1/(ds1+ds2)
                        w2 = ds2/(ds1+ds2)
                     
                    # Open the arc solution files 
                    fn0 = os.path.join(out_dir, '%s.wsol.fits' % (local_arcs[0]))
                    fn1 = os.path.join(out_dir, '%s.wsol.fits' % (local_arcs[1]))
                    fits0 = pyfits.open(fn0)
                    fits1 = pyfits.open(fn1)


                    for i in range(1,len(fits0)):
                         fits0[i].data = w1*fits0[i].data + w2*fits1[i].data

                    wsol_fn = os.path.join(out_dir, '%s.wsol.fits' % (fn) )
                    fits0.writeto(wsol_fn, overwrite=True)

                    print('(2 arcs found)')
                    print('(Note: using %sx%s.wsol.fits + %sx%s.wsol.fits as wsol file)' % (np.round(w1,2),local_arcs[0],np.round(w2,2),local_arcs[1]))
                           
                else:
                    # Arcs do not surround the Science frame
                    # Revert to using the first one instead
                    wsol_fn = os.path.join(out_dir, '%s.wsol.fits' % (local_arcs[0]))
                    print('(2 arcs found, but they do not bracket the Science frame!)')
                    print('(Note: using %s as wsol file)' % wsol_fn.split('/')[-1])
                    
            else:
                # Either 1 or more than two arcs present ... only use the first one !
                wsol_fn = os.path.join(out_dir, '%s.wsol.fits' % (local_arcs[0]))
                print('(Note: using %s as wsol file)' % wsol_fn.split('/')[-1] )

        else:
            wsol_fn = wsol_out_fn

        # All done, let's generate the cube
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
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.x%s.dat'  % (fn, prev_suffix))
        print('Extract %s standard star from %s' % (type, in_fn.split('/')[-1]))
        pywifes.extract_wifes_stdstar(in_fn,
                                      save_fn=out_fn,
                                      save_mode='ascii',
                                      **args)
    return

# Sensitivity Function fit
def run_derive_calib(metadata, prev_suffix, curr_suffix, method = 'poly',**args):
    std_obs_list = get_primary_std_obs_list(metadata, type='flux')
    std_cube_list = [os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
                     for fn in std_obs_list]
    extract_list = [os.path.join(out_dir, '%s.x%s.dat' % (fn, prev_suffix))
                    for fn in std_obs_list]
    print('Deriving sensitivity function')
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
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        print('Flux-calibrating cube %s' % in_fn.split('/')[-1])
        pywifes.calibrate_wifes_cube(
            in_fn, out_fn, calib_fn, mode)
    return

#------------------------------------------------------
# Telluric - derive
def run_derive_telluric(metadata, prev_suffix, curr_suffix, **args):
    std_obs_list = get_primary_std_obs_list(metadata, 'telluric')
    std_cube_list = [os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
                     for fn in std_obs_list]
    extract_list = [os.path.join(out_dir, '%s.x%s.dat' % (fn, prev_suffix))
                    for fn in std_obs_list]
    print('Deriving telluric correction')
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
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        print('Correcting telluric in %s' % in_fn.split('/')[-1])
        pywifes.apply_wifes_telluric(
            in_fn, out_fn, tellcorr_fn)
    return

#------------------------- Fred's update -----------------
# Save final cube in suitable fits format
def run_save_3dcube(metadata, prev_suffix, curr_suffix, **args):
    # now generate cubes
    sci_obs_list = get_primary_sci_obs_list(metadata)
    for fn in sci_obs_list:
        in_fn  = os.path.join(out_dir, '%s.p%s.fits' % (fn, prev_suffix))
        out_fn = os.path.join(out_dir, '%s.p%s.fits' % (fn, curr_suffix))
        print('Saving 3D Data Cube for %s' % in_fn.split('/')[-1])
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
        #~ try:
            #~ print('$$$ WHAT TO REDUCE ', step_name, get_primary_sci_obs_list(obs_metadata))
        #~ except:
            #~ print('$$$ WHAT TO REDUCE ', step_name)
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
print('All done in %.01f seconds.' % duration.total_seconds())
