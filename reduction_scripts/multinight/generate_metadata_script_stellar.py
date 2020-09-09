#! /usr/bin/env python

"""
Example:
python generate_metadata_script_marusa.py config.py /priv/mulga1/marusa/2m3data/20190114
"""

import numpy as np
import sys, getopt
import os
import glob
from astropy.io import fits as pyfits
import wifes_calib
sys.path.insert(1, '/data/mash/marusa/reduction_wifes/pipeline/utils/')
import pywifes_utils as pu
import imp

keywords = ['GRATINGB', 'GRATINGR', 'BEAMSPLT', 'CCDSUM', 'CCDSEC']
keywords_dark_zero = ['CCDSEC', 'CCDSUM']

obsdate = sys.argv[2]

print('########################################################')
print('########################################################')
print('################### %s ###########################'%obsdate)
print('########################################################')
print('########################################################')


# Config file (should be in the output folder)
config = imp.load_source('config', sys.argv[1])

# Calibration files in case some are missing for this night: this is a file containing a dictionary with filenames
sys.path.insert(1, '/data/mash/marusa/2m3data/wifes/') # TODO: THIS SHOULDN'T BE LIKE THIS!!
import calibration_all_filenames as cal
cal=cal.result

# Some nights have unusually high bias levels. Print a warning if reducing such a night.
# Also, be careful if including calibrations from such a night!
try:
    if config.badcalib_filename is not None:
        badcalib = np.loadtxt(config.badcalib_filename, comments='#', dtype=int)
        badcalib_obsdates = set(badcalib[:,0])
        if int(obsdate) in badcalib_obsdates:
            print('WARNING WARNING WARNING: Biases in this night have unusually high values!')
except:
    pass

# List of standard stars
stdstar_list = wifes_calib.ref_fname_lookup.keys()
stdstar_is_flux_cal = wifes_calib.ref_flux_lookup
stdstar_is_telluric = wifes_calib.ref_telluric_lookup


# Input folder with raw data
data_dir = os.path.join(config.input_root, obsdate)

print('#'+54*'-')
print('data_dir', data_dir)

###### OUTPUT FOLDER ################
root_obsdate = os.path.join(config.output_root, '%s'%obsdate)

# Create folder with date
root_bool = os.path.isdir(root_obsdate) and os.path.exists(root_obsdate)
if not root_bool:
    os.mkdir(root_obsdate)
print('root_obsdate', root_obsdate)

# Add band (grating) to the output folder name
out_dir = os.path.join(root_obsdate, 'reduced_%s'%config.band)

prefix=config.prefix

if prefix is not None and len(prefix)>0:
    print('prefix', prefix)
    out_dir += '_%s'%prefix
out_dir_bool = os.path.isdir(out_dir) and os.path.exists(out_dir)
if not out_dir_bool:
    os.mkdir(out_dir)
print('out_dir', out_dir)

print('#'+54*'-')
print('prefix: %s'%prefix)
print('#'+54*'-')

# Get list of all fits files in directory
# Take only WiFeS frames
all_files = os.listdir(data_dir)
all_files = [os.path.join(data_dir, x) for x in all_files if x.endswith('.fits') and 'T2m3ag' not in x and 'T2m3Ec'.lower() not in x.lower()]


def take_only_specific_objects(object_list, all_files):
    """
    If you have a folder with many different objects but want to reduce only a few specific
    stars, then specify them in the object_list
    
    Return:
    all_files without science frames that are not needed here
    
    """
    
    all_files2=[]
    for fn in all_files:
        try:
            header = pyfits.getheader(fn, 0)
        except:
            print('Cant open file', fn)
            continue
        try:
            objectname=header['OBJNAME']
            if objectname.replace(' ', '') not in object_list:
                continue
            else:
                all_files2.append(fn)
        except: # calibration files
            all_files2.append(fn)
    
    print('all_files, all_files2', len(all_files), len(all_files2))
    return all_files2



try:
    object_list = np.loadtxt(config.object_list_filename)
    print('Object list read from', config.object_list_filename)
    object_list = [x.replace(' ', '') for x in object_list]
    all_files = take_only_specific_objects(object_list, all_files)
    print('Taken only selected objectnames.')
except:
    try:
        object_list = config.object_list
        object_list = [x.replace(' ', '') for x in object_list]
        all_files = take_only_specific_objects(object_list, all_files)
        print('Taken only selected objectnames.')
    except:
        object_list=None

print('object_list', object_list)
print('ALLFILES ', len(all_files))

print('Excluding bad frames...')
all_files = pu.exclude_bad_frames(all_files, config.excluderun_filename)

print('Ready to start with %d files.'%len(all_files))

### IF CALIBRATION FILES ARE MISSING ##############
# In case calibration files are missing, specify dates (format 20190315) for each cal type
selected_cal_dates={}
try:
    opts, args = getopt.getopt(sys.argv[3:], "d:b:f:")
except getopt.GetoptError:
    print('test.py -i <inputfile> -o <outputfile>') # TODO
for opt, arg in opts:
    if opt in ("-d"):
        selected_cal_dates['DARK']=int(arg)
    elif opt in ("-b"):
        selected_cal_dates['BIAS']=int(arg)
    elif opt in ("-f"):
        selected_cal_dates['FLAT']=int(arg)
    
print 'SELECTED_CAL_DATES', selected_cal_dates
###################################################





def classify_files_into_modes():
    """
    Find all different modes (settings) taken during this night.
    Classify filenames into different modes
    
    Return:
    Dictionary: modes[mode] = [list of filenames]
    
    """

    modes=dict()
    
    for fn in all_files:
        try:
            header = pyfits.getheader(fn, 0)
        except:
            print('Cant open file', fn)
            continue

        try:
            k=[header[x] for x in keywords]
        except:
            print 'Cannot get full header for', fn, header['IMAGETYP']
            continue
        
        k=tuple(k) # Lists or sets cannot be dictionary keys
        try:
            modes[k].append(fn)
        except:
            modes[k]=[fn]

    if len(modes)<1:
        print('WARNING: 0 modes!')
        
    return modes

def classify_filenames(all_files):
    blue_obs = []
    red_obs = []
    obs_date = None
    for fn in all_files:
        obs = fn.replace('.fits', '').split('/')[-1]
        
        try:
            header = pyfits.getheader(fn, 0)
        except:
            continue
        
        # Obsdate
        if obs_date == None and header['IMAGETYP'].lower()=='object':
            obs_date = header['DATE-OBS'].split('T')[0].replace('-', '')             

        # Camera
        camera = header['CAMERA']
        
        # Classify into blue/red
        if camera == 'WiFeSBlue':
            if obs in blue_obs:
                continue
            else:
                blue_obs.append(obs)
        if camera == 'WiFeSRed':
            if obs in red_obs:
                continue
            else:
                red_obs.append(obs)
        
    return blue_obs, red_obs, obs_date

def classify_frames_into_imagetypes(frames=None):
    # classify each obs
    bias = []
    domeflat = []
    twiflat = []
    dark = []
    arc = []
    wire = []
    stdstar = {}
    science = {}

    objects={} # Needed to match science and arc
    arcs=[]

    for obs in sorted(frames):
        fn = os.path.join(data_dir, obs+'.fits')
        header = pyfits.getheader(fn, 0)
        
        
        try:
            obj_name = header['OBJNAME'].replace(' ', '')
        except:
            obj_name=None # zero (bias) exposure
        try:
            mjd = header['MJD-OBS']
        except:
            mjd=None
        run = header['RUN']
        imagetype = header['IMAGETYP'].upper()
        exptime = header['EXPTIME']
        
        print(fn, imagetype, obj_name)


        
        #~ if imagetype=='OBJECT':
            #~ if objectnames: # keep only selected objects
                #~ if obj_name not in objectnames:
                    #~ continue
            #~ if exclude_objectnames: # exclude selected objects
                #~ if obj_name in exclude_objectnames:
                    #~ continue
        
        
        #~ print ccdsumf, naxis2
        #---------------------------
        # check if it is within a close distance to a standard star
        # if so, fix the object name to be the good one from the list!
        
        # TODO: This can lead to errors!
        
        try:
            near_std, std_dist = wifes_calib.find_nearest_stdstar(fn)
            if std_dist < 100.0:
                obj_name = near_std
        except:
            pass
            
            
        #---------------------------
        # 1 - bias frames
        if imagetype.upper() == 'ZERO' or imagetype.upper() == 'BIAS':
            bias.append(obs)
        # 2 - quartz flats
        if imagetype == 'FLAT':
            domeflat.append(obs)
        # 3 - twilight flats
        if imagetype == 'SKYFLAT':
            twiflat.append(obs)
        # 4 - dark frames
        if imagetype == 'DARK':
            dark.append(obs)
        # 5 - arc frames
        if imagetype == 'ARC':
            arc.append(obs)
            
            # For arc-science matching
            arcs.append([run, mjd, imagetype, obs])
            
        # 6 - wire frames
        if imagetype == 'WIRE':
            wire.append(obs)
        # all else are science targets
        if imagetype == 'OBJECT':
            if obj_name in stdstar_list:
                # group standard obs together!
                if obj_name in stdstar.keys():
                    stdstar[obj_name].append(obs)
                else:
                    stdstar[obj_name] = [obs]
            else:
                # group science obs together!
                if obj_name in science.keys():
                    science[obj_name].append(obs)
                else:
                    science[obj_name] = [obs]

            # For arc-science matching
            if obj_name in objects.keys():
                objects[obj_name].append([run, mjd, exptime, imagetype, obs])
            else:
                objects[obj_name]=[[run, mjd, exptime, imagetype, obs]]

    arcs_per_star = match_object_and_arc(objects=objects, arcs=arcs)    
    
    return science, bias, domeflat, twiflat, dark, arc, wire, stdstar, arcs_per_star #, objects

# Marusa: Match science images and arcs. One arc before and one after the science exposure + all arcs in between.
# Disadvantage: if you e.g. observe 1 objects with arcs, then do other things and observe it later again, everything is combined together.
def match_object_and_arc(objects=None, arcs=None):
    """
    Find one arc before and one arc after science exposures + all arcs in between. Based on MJD comparison.
    
    objects: science exposures, format: object[objname] = [[run, MJD, exptime, imagetyp, filename]]]
    
    arcs: a list of arcs [run, MJD, imagetyp, filename]
    
    Note: This is not done with run numbers because they don't necessarily increase during the night - e.g. they reset to 0 at TAROS restart.
    """

    result={}
    
    if len(arcs)<1:
        print('ARCS:::', arcs)
    
    #~ arc_mjd=np.array([x[1] for x in arcs])

    for k, v in objects.items():
        # Sort all exposured for this object this night
        v=sorted(v)

        mjd_first_science = v[0][1] # first science exposure of this object this night
        mjd_last_science = v[-1][1] # last science exposure of this object this night
        exptime_last=v[-1][2]/3600.0/24.0 # convert exptime to days
        
        # Last arc must be taken after the last science exposure was finished
        mjd_last_science+=exptime_last


        # make sure the science and arc exposures were taken on the same night
        date_object = v[0][-1].split('.')[0].split('-')[1]
        arcs_tmp = [x for x in arcs if x[-1].split('.')[0].split('-')[1]==date_object] # Take only arcs on this night
        arc_mjd=np.array([x[1] for x in arcs_tmp])

        a=[]

        # Find arc taken just before the first science exposure
        diff_first=mjd_first_science-arc_mjd # Take first positive
        mask_first = diff_first>0
        if len(diff_first[mask_first])>0:
            value_first=sorted(diff_first[mask_first])[0]
            index_first=np.where(diff_first==value_first)[0][0]
            a=[arcs_tmp[index_first][-1]]
        else:
            pass # First image of the night. First arc was taken later.

        # Find arc taken just after the last science exposure.
        diff_last=arc_mjd-mjd_last_science # Take closest to 0 from the negative side
        mask_last = diff_last>0
        if len(diff_last[mask_last])>0:
            value_last=sorted(diff_last[mask_last])[0]
            index_last=np.where(diff_last==value_last)[0][0]
            a.append(arcs_tmp[index_last][-1])
        else:
            pass # No arc at the end of the night. This shouldn't happen.
        
        # Check if there were any other arcs between first and last science exposure.
        mask = (arc_mjd>mjd_first_science) & (arc_mjd<mjd_last_science)
        for m, ar in zip(mask, arcs_tmp):
            if m:
                a.append(ar[-1])

        if len(a)>0:
            result[k]=a
        else:
            print('WARNING: NO ARC FOUND FOR %s %s.'%(k, v[0][-1]))
    
    return result
    
def test_if_all_essential_files_are_available(camera=None, science=None, arcs=None, dark=None, bias=None, flat_dome=None, flat_twi=None, std_obs=None, wire=None):
    result={}
    
    if len(science)<1:
        print('**** WARNING (%s): No science frames found.'%camera)
        
    if len(arcs)<1:
        print('**** WARNING (%s): No arc frames found.'%camera)
        
    if len(bias)<config.calmin:
        print('**** WARNING (%s): No bias frames found.'%camera)
        bias_key=False
        result['BIAS']=False
    else:
        bias_key=True
        result['BIAS']=bias
        
    if len(dark)<config.calmin:
        print('WARNING (%s): No dark frames found.'%camera)
        dark_key=False
        result['DARK']=False
    else:
        dark_key=True
        result['DARK']=dark
        
    if len(flat_dome)<config.calmin:
        print('**** WARNING (%s): No dome flat frames found.'%camera)
        flat_dome_key=False
        result['FLAT']=False
    else:
        flat_dome_key=True
        result['FLAT']=flat_dome
        
    if len(flat_twi)<1:
        print('WARNING (%s): No twilight flat frames found.'%camera)
        
    if len(wire)<1:
        print('WARNING (%s): No wire frames found.'%camera)
        
    if len(std_obs)<1:
        print('**** WARNING (%s): No std_obs frames found.'%camera)
    
    # Standards
    if len(std_obs)>0:
        flux=False
        telluric=False
        for x in std_obs:
            if stdstar_is_flux_cal[x]:
                flux=True
            if stdstar_is_telluric[x]:
                telluric=True
        if not flux:
            print('**** WARNING (%s): No FLUX std_obs frames found.'%camera)
        if not telluric:
            print('**** WARNING (%s): No TELLURIC std_obs frames found.'%camera)
    
    
    
    return result
   
def propose_missing_calib_files(mode=None, calstat=None):
    """
    Find corresponding calibration files in other nights.
    
    calstat: stats on calibrations: calstat[imagetype]=False/len(images)
    """

    keyword_indices = [keywords.index(x) for x in keywords_dark_zero]
    
    # Mode (settings) for darks and zeros because the requirements are not that strict
    mode_dark_zero = tuple([mode[x] for x in keyword_indices])
    
    # What calib files are missing?
    for imagetype, status in calstat.iteritems():
        if status: # Set lower limit on number of calib files needed.
            if len(status)<config.calmin:
                missing=True
            else:
                missing=False
        else:
            missing=True
        
        if missing: # Missing. Find them. What if c==None?
            print('Missing', imagetype)
            if imagetype.upper() not in ['DARK', 'ZERO', 'BIAS']:
            #~ if imagetype.upper() not in ['ZERO', 'BIAS']:
                c=cal[mode]
            else:
                c=cal[mode_dark_zero]
                
            try:
                dates=c[imagetype] # dates available for this particular imagetype
                print('Missing %s calibration file. Available:'%imagetype)
                for k, v in dates.iteritems():
                    if k in badcalib_obsdates:
                        print(k, len(v), '(This night might have biases with unusually high values!)') # len is both blue and red!
                    else:
                        print(k, len(v)) # len is both blue and red!
                    #~ filenames=v[date]

                print()
            except:
                if imagetype=='BIAS': # zero instead of bias
                    try:
                        dates=c['ZERO'] # dates available for this particular imagetype
                        print('Missing %s calibration file. Available:'%imagetype)
                        for k, v in dates.iteritems():
                            if k in badcalib_obsdates:
                                print(k, len(v), '(This night might have biases with unusually high values!)') # len is both blue and red!
                            else:
                                print(k, len(v)) # len is both blue and red!
                            #~ filenames=v[date]

                        print()
                    except:
                        print('*** (case1) NO %s AVAILABLE FOR THIS MODE!!!'%imagetype)

                else:
                    print('*** (case2) NO %s AVAILABLE FOR THIS MODE!!!'%imagetype)

def include_missing_calib_files(mode=None, calstat=None, camera=None):
    """
    calstat: stats on calibrations: calstat[imagetype]=False/len(images)
    Include just the correct name. The path to the files is reconstructed in the reduction code from the date in the filename.
    """
    result={}

    keyword_indices = [keywords.index(x) for x in keywords_dark_zero]
    
    mode_dark_zero = tuple([mode[x] for x in keyword_indices])


    # What calib files are missing?
    for imagetype, status in calstat.iteritems():
        if imagetype.upper() not in ['DARK', 'ZERO', 'BIAS']:
            c=cal[mode]
        else:
            c=cal[mode_dark_zero]
                
        if imagetype not in selected_cal_dates:
            continue
        else:
            if imagetype not in ['BIAS', 'ZERO']:
                date_wanted = selected_cal_dates[imagetype] # check BIAS-ZERO
            else:
                try:
                    date_wanted = selected_cal_dates['ZERO'] # check BIAS-ZERO
                except:
                    date_wanted = selected_cal_dates['BIAS'] # check BIAS-ZERO
            
        if status: # Calibration files available, so don't do anything. BUT: Set lower limit on number of calib files needed.
            if len(status)<3:
                missing=True
                TODO=True
        
        else: # Missing. Find them. What if c==None?
            try:
                dates=c[imagetype] # dates available for this particular imagetype
                filenames = dates[date_wanted]
                
                # Take only a selected band (blue or red)
                if camera=='WiFeSRed':
                    filenames = [x for x in filenames if 'T2m3wr' in x]
                elif camera=='WiFeSBlue':
                    filenames = [x for x in filenames if 'T2m3wb' in x]
                
                # Delete path. It is added later in the reduction code
                filenames = [x.split('/')[-1].replace('.fits', '') for x in filenames]
                
                print('Adding %d %s images:'%(len(filenames), imagetype))
                result[imagetype]=filenames
                for x in filenames:
                    print(x)
                print()
            except:
                if imagetype=='BIAS': # zero instead of bias
                    try:
                        dates=c['ZERO'] # dates available for this particular imagetype
                        filenames = dates[date_wanted]

                        # Take only a selected band (blue or red)
                        if camera=='WiFeSRed':
                            filenames = [x for x in filenames if 'T2m3wr' in x]
                        elif camera=='WiFeSBlue':
                            filenames = [x for x in filenames if 'T2m3wb' in x]

                        filenames = [x.split('/')[-1].replace('.fits', '') for x in filenames]

                        print('Adding %d %s images:'%(len(filenames), imagetype))
                        result[imagetype]=filenames
                        for x in filenames:
                            print(x)
                        print()
                    except:
                        print('*** (include case1) NO %s AVAILABLE FOR THIS MODE!!!'%imagetype)

                else:
                    print('*** (include case2) NO %s AVAILABLE FOR THIS MODE!!!'%imagetype)    

    return result
 
def write_metadata(science=None, bias=None, domeflat=None, twiflat=None, dark=None, arc=None, arcs_per_star=None, wire=None, camera=None, std_obs=None, nmode=0, kmode=None, number_of_modes=0):
    if len(science)>0:
        pass
    else:
        print('No science images for the mode', kmode)
        return False
    #------------------------------------------------------
    # write to metadata save script!

    # filename
    if prefix is not None and len(prefix)>0:
        metadata_filename=os.path.join(out_dir, '%s_mode_%d_metadata_%s.py'%(prefix, nmode, camera)) # TODO
    else:
        if number_of_modes==1:
            metadata_filename=os.path.join(out_dir, 'metadata_%s.py'%camera)
        elif number_of_modes>1 and nmode==0:
            metadata_filename=os.path.join(out_dir, 'metadata_%s.py'%camera)
        elif number_of_modes>1 and nmode>0:
            metadata_filename=os.path.join(out_dir, 'mode_%d_metadata_%s.py'%(nmode, camera))

    f = open(metadata_filename, 'w')
    
    dsplit = '#' + 54*'-' + '\n'

    #------------------
    # mode
    f.write(dsplit)
    f.write('mode = '+str(kmode)+'\n\n')

    #------------------
    # calibrations
    f.write(dsplit)

    # 1 - bias
    f.write('bias_obs = [\n')
    for obs in bias:
        f.write('   \'%s\', \n' % obs)
    f.write('    ]\n')
    f.write('\n')

    # 2 - domeflat
    f.write('domeflat_obs = [\n')
    for obs in domeflat:
        f.write('   \'%s\', \n' % obs)
    f.write('    ]\n')
    f.write('\n')

    # 3 - twiflat
    f.write('twiflat_obs = [\n')
    for obs in twiflat:
        f.write('   \'%s\', \n' % obs)
    f.write('    ]\n')
    f.write('\n')

    # 4 - dark
    f.write('dark_obs = [\n')
    for obs in dark:
        f.write('   \'%s\', \n' % obs)
    f.write('    ]\n')
    f.write('\n')

    # 5 - arc
    f.write('arc_obs = [\n')
    for obs in arc:
        f.write('   \'%s\', \n' % obs)
    f.write('    ]\n')
    f.write('\n')

    # 6 - wire
    f.write('wire_obs = [\n')
    for obs in wire:
        f.write('   \'%s\', \n' % obs)
    f.write('    ]\n')
    f.write('\n')

    #------------------
    # science
    #~ f.write(dsplit)
    f.write('sci_obs = [\n')
    for obj_name in science.keys():
        obs_list = science[obj_name]
        obs_str = '\'%s\'' % obs_list[0]
        for i in range(1,len(obs_list)):
            obs = obs_list[i]
            obs_str += ',\n               \'%s\'' % obs
        obs_list_arc = arcs_per_star[obj_name]
        obs_arc_str = '\'%s\'' % obs_list_arc[0]
        for i in range(1,len(obs_list_arc)):
            obs = obs_list_arc[i]
            obs_arc_str += ',\n               \'%s\'' % obs
        f.write('    # %s\n' % obj_name)
        f.write('    {\'sci\'  : [%s],\n' % obs_str)
        f.write('     \'sky\'  : [],\n')
        f.write('     \'arc\'  : [%s]},\n' % obs_arc_str)


    f.write('    ]')
    f.write('\n\n')

    #------------------
    # stdstars
    #~ f.write(dsplit)
    f.write('std_obs = [')
    for obj_name in stdstar.keys():
        obs_list = stdstar[obj_name]
        obs_str = '\'%s\'' % obs_list[0]
        for i in range(1,len(obs_list)):
            obs = obs_list[i]
            obs_str += ',\n               \'%s\'' % obs
         # Also write arcs for tellurics 
        obs_list_arc = arcs_per_star[obj_name]
        obs_arc_str = '\'%s\'' % obs_list_arc[0]
        for i in range(1,len(obs_list_arc)):
            obs = obs_list_arc[i]
            obs_arc_str += ',\n               \'%s\'' % obs

        f.write('    # %s\n' % obj_name)
        f.write('    {\'sci\'  : [%s],\n' % obs_str)
        f.write('     \'name\' : [\'%s\'],\n' % obj_name)
        # Star is both telluric and flux standard
        if stdstar_is_flux_cal[obj_name] and stdstar_is_telluric[obj_name]:
            f.write('     \'type\' : [\'flux\', \'telluric\'],\n')
        # Star is flux standard but not telluric
        elif stdstar_is_flux_cal[obj_name] and not stdstar_is_telluric[obj_name]:
            f.write('     \'type\' : [\'flux\'],\n')
        # Star is telluric but not flux standard
        elif not stdstar_is_flux_cal[obj_name] and stdstar_is_telluric[obj_name]:
            f.write('     \'type\' : [\'telluric\'],\n')
        f.write('     \'arc\'  : [%s]},\n' % obs_arc_str)
    f.write('    ]\n')
    f.write('\n')

    #------------------
    # footers
    f.write(dsplit)
    #~ out_fn = 'wifesB_%s_metadata.pkl' % obs_date
    f.write('night_data = {\n')
    f.write('    \'bias\' : bias_obs,\n')
    f.write('    \'domeflat\' : domeflat_obs,\n')
    f.write('    \'twiflat\' : twiflat_obs,\n')
    f.write('    \'dark\' : dark_obs,\n')
    f.write('    \'wire\' : wire_obs,\n')
    f.write('    \'arc\'  : arc_obs,\n')
    f.write('    \'sci\'  : sci_obs,\n')
    f.write('    \'std\'  : std_obs}\n')
    f.write('\n')
    #~ f.write('f1 = open(\'%s\', \'w\')' % out_fn
    #~ f.write('pickle.dump(night_data, f1)'
    #~ f.write('f1.close()'

    f.close()
    
    print('########################################################')
    print('METADATA written in', metadata_filename)
    print('########################################################')
    return True

   
if __name__ == '__main__':
    # Data might be taken with different settings (modes)
    # 'modes' is a dictionary with modes as keys and filenames as values
    # Classify filenames into different modes
    modes = classify_files_into_modes()
    
    # Sort modes by the number of images
    keys = [k for k in modes.iterkeys()]
    if len(modes)>1:
        lens = [len(v) for v in modes.itervalues()]
        indices = list(np.argsort(lens))
        keys = [keys[index] for index in indices]
        keys=keys[::-1] # reverse order
    if len(modes)<1:
        print('NO MODES FOUND! check your data')

    # Prepare metadata file for each mode    
    M=len(modes)
    m=0
    #~ for mode, filenames in modes.iteritems():
    for mode in keys:
        filenames = modes[mode]
        # Classify filenames into blue/red
        print('CLASSIFY'), mode
        blue_obs, red_obs, obs_date = classify_filenames(filenames)
        print
        print
        print
        print
        print
    
        # Classify filenames into imagetypes
        if config.band == 'r':
            camera='WiFeSRed'
            science, bias, domeflat, twiflat, dark, arc, wire, stdstar, arcs_per_star = classify_frames_into_imagetypes(frames=red_obs)
        elif config.band == 'b':
            camera='WiFeSBlue'
            science, bias, domeflat, twiflat, dark, arc, wire, stdstar, arcs_per_star = classify_frames_into_imagetypes(frames=blue_obs)
        
                    
        
        # DATA CHECKS
        if len(science)>0:
            print
            print
            print '########################################################'
            print 'Mode', mode
            print '########################################################'


            print('#'+54*'-')
            #~ print('OBSDATE', obsdate)
            print('%s: %d objects found:'%(camera, len(science)))
            for k in science.keys():
                print(k)
            print('#'+54*'-')
        else:
            print('No science data found', mode)
            continue # Just ignore modes with no science data
        
        # Availability of calibration files
        calstat = test_if_all_essential_files_are_available(camera=camera, science=science, arcs=arc, dark=dark, bias=bias, flat_dome=domeflat, flat_twi=twiflat, std_obs=stdstar, wire=wire)

        try:
            # If calibration files are missing, find them in other nights          
            propose_missing_calib_files(mode=mode, calstat=calstat)
            
            # Update array with missing data. Only include correct filenames (data is not copied here!)
            missing_cal = include_missing_calib_files(mode=mode, calstat=calstat, camera=camera)
            for imagetype, filenames in missing_cal.iteritems():
                if imagetype=='BIAS':
                    bias=filenames
                    print('new biases'), bias
                elif imagetype=='DARK':
                    dark=filenames
                    print('new darks'), dark
                elif imagetype=='FLAT':
                    domeflat=filenames
                    print('new flats'), domeflat
                else:
                    print("Update missing data: We've got a problem here.", imagetype)
        except:
            print('There is no additional calibration data available for this mode.')
        
        
        # SORT filenames
        bias=sorted(bias)
        domeflat=sorted(domeflat)
        dark=sorted(dark)
        arc=sorted(arc)
        

        # WRITE METADATA        
        success = write_metadata(camera=camera, science=science, bias=bias, domeflat=domeflat, twiflat=twiflat, dark=dark, arc=arc, arcs_per_star=arcs_per_star, wire=wire, std_obs=stdstar, nmode=m, kmode=mode, number_of_modes=M)
        

        m+=1 # mode numbers
