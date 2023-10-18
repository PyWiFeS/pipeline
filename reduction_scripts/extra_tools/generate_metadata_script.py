#! /usr/bin/env python

import sys
import os
from astropy.io import fits as pyfits
import wifes_calib

stdstar_list = wifes_calib.ref_fname_lookup.keys()

# get directory to be parsed
try:
    data_dir = os.path.abspath(sys.argv[1])+'/'
except:
    data_dir = os.getcwd()+'/'
try:
    naxis2_use = int(sys.argv[2])
except:
    naxis2_use = 0

# get list of all fits files in directory
all_files = os.listdir(data_dir)
blue_obs = []
red_obs = []
obs_date = None
for fn in all_files:
    obs = fn.replace('.fits', '')
    if obs_date == None:
        try:
            f = pyfits.open(data_dir+fn)
            obs_date = f[0].header['DATE-OBS'].split('T')[0].replace('-', '')
            f.close()
        except:
            continue
        #obs_date = obs[7:15]
    try:
        f = pyfits.open(data_dir+fn)
        camera = f[0].header['CAMERA']
        naxis2 = f[0].header['NAXIS2']
        f.close()
    except:
        continue
    if naxis2_use!=0:
        if naxis2 != naxis2_use:
            continue
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
#print blue_obs
#print red_obs
#print obs_date

#---------------------------------------------
#  ***  BLUE CHANNEL  ***
#---------------------------------------------
# classify each obs
blue_bias = []
blue_domeflat = []
blue_twiflat = []
blue_dark = []
blue_arc = []
blue_wire = []
blue_stdstar = {}
blue_science = {}

for obs in blue_obs:
    fn = data_dir+obs+'.fits'
    f = pyfits.open(fn)
    imagetype = f[0].header['IMAGETYP'].upper()
    obj_name = f[0].header['OBJECT']
    f.close()
    #---------------------------
    # check if it is within a close distance to a standard star
    # if so, fix the object name to be the good one from the list!
    try:
        near_std, std_dist = wifes_calib.find_nearest_stdstar(fn)
        if std_dist < 100.0:
            obj_name = near_std
    except:
        pass
    #---------------------------
    # 1 - bias frames
    if (imagetype == 'BIAS') or (imagetype == 'ZERO'):
        blue_bias.append(obs)
    # 2 - quartz flats
    if imagetype == 'FLAT':
        blue_domeflat.append(obs)
    # 3 - twilight flats
    if imagetype == 'SKYFLAT':
        blue_twiflat.append(obs)
    # 4 - dark frames
    if imagetype == 'DARK':
        blue_dark.append(obs)
    # 5 - arc frames
    if imagetype == 'ARC':
        blue_arc.append(obs)
    # 6 - wire frames
    if imagetype == 'WIRE':
        blue_wire.append(obs)
    # all else are science targets
    if imagetype == 'OBJECT':
        if obj_name in stdstar_list:
            # group standard obs together!
            if obj_name in blue_stdstar.keys():
                blue_stdstar[obj_name].append(obs)
            else:
                blue_stdstar[obj_name] = [obs]
        else:
            # group science obs together!
            if obj_name in blue_science.keys():
                blue_science[obj_name].append(obs)
            else:
                blue_science[obj_name] = [obs]

#------------------------------------------------------
# write to metadata save script!
f = open('save_blue_metadata.py', 'w')

dsplit = '#' + 54*'-'

#------------------
# headers
f.write('import pickle' + '\n')
f.write('' + '\n')

#------------------
# calibrations
f.write(dsplit+ '\n')

# 1 - bias
f.write('bias_obs = [' + '\n')
for obs in blue_bias:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

# 2 - domeflat
f.write('domeflat_obs = [' + '\n')
for obs in blue_domeflat:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

# 3 - twiflat
f.write('twiflat_obs = [' + '\n')
for obs in blue_twiflat:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

# 4 - dark
f.write('dark_obs = [' + '\n')
for obs in blue_dark:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

# 5 - arc
f.write('arc_obs = [' + '\n')
for obs in blue_arc:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

# 6 - wire
f.write('wire_obs = [' + '\n')
for obs in blue_wire:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

#------------------
# science
f.write(dsplit+ '\n')
f.write('sci_obs = [' + '\n')
for obj_name in blue_science.keys():
    obs_list = blue_science[obj_name]
    obs_str = '\'%s\'' % obs_list[0]
    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    f.write('    # %s' % obj_name + '\n')
    f.write('    {\'sci\'  : [%s],' % obs_str + '\n')
    f.write('     \'sky\'  : []},' + '\n')


f.write('    ]' + '\n')
f.write('' + '\n')

#------------------
# stdstars
f.write(dsplit+ '\n')
f.write('std_obs = [' + '\n')
for obj_name in blue_stdstar.keys():
    obs_list = blue_stdstar[obj_name]
    obs_str = '\'%s\'' % obs_list[0]
    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    f.write('    # %s' % obj_name + '\n')
    f.write('    {\'sci\'  : [%s],' % obs_str + '\n')
    f.write('     \'name\' : [\'%s\'],' % obj_name + '\n')
    f.write('     \'type\' : [\'flux\', \'telluric\']},' + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

#------------------
# footers
f.write(dsplit+ '\n')
out_fn = 'wifesB_%s_metadata.pkl' % obs_date
f.write('night_data = {' + '\n')
f.write('    \'bias\' : bias_obs,' + '\n')
f.write('    \'domeflat\' : domeflat_obs,' + '\n')
f.write('    \'twiflat\' : twiflat_obs,' + '\n')
f.write('    \'dark\' : dark_obs,' + '\n')
f.write('    \'wire\' : wire_obs,' + '\n')
f.write('    \'arc\'  : arc_obs,' + '\n')
f.write('    \'sci\'  : sci_obs,' + '\n')
f.write('    \'std\'  : std_obs}' + '\n')
f.write('' + '\n')
f.write('f1 = open(\'%s\', \'wb\')' % out_fn + '\n')
f.write('pickle.dump(night_data, f1)' + '\n')
f.write('f1.close()' + '\n')

f.close()

#---------------------------------------------
#  ***  RED CHANNEL  ***
#---------------------------------------------
# classify each obs
red_bias = []
red_domeflat = []
red_twiflat = []
red_dark = []
red_arc = []
red_wire = []
red_stdstar = {}
red_science = {}

for obs in red_obs:
    fn = data_dir+obs+'.fits'
    f = pyfits.open(fn)
    imagetype = f[0].header['IMAGETYP'].upper()
    obj_name = f[0].header['OBJECT']
    f.close()
    #---------------------------
    # check if it is within a close distance to a standard star
    # if so, fix the object name to be the good one from the list!
    try:
        near_std, std_dist = wifes_calib.find_nearest_stdstar(fn)
        if std_dist < 100.0:
            obj_name = near_std
    except:
        pass
    #---------------------------
    # 1 - bias frames
    if (imagetype == 'BIAS') or (imagetype == 'ZERO'):
        red_bias.append(obs)
    # 2 - quartz flats
    if imagetype == 'FLAT':
        red_domeflat.append(obs)
    # 3 - twilight flats
    if imagetype == 'SKYFLAT':
        red_twiflat.append(obs)
    # 4 - dark frames
    if imagetype == 'DARK':
        red_dark.append(obs)
    # 5 - arc frames
    if imagetype == 'ARC':
        red_arc.append(obs)
    # 6 - wire frames
    if imagetype == 'WIRE':
        red_wire.append(obs)
    # all else are science targets
    if imagetype == 'OBJECT':
        if obj_name in stdstar_list:
            # group standard obs together!
            if obj_name in red_stdstar.keys():
                red_stdstar[obj_name].append(obs)
            else:
                red_stdstar[obj_name] = [obs]
        else:
            # group science obs together!
            if obj_name in red_science.keys():
                red_science[obj_name].append(obs)
            else:
                red_science[obj_name] = [obs]


#------------------------------------------------------
# write to metadata save script!
f = open('save_red_metadata.py', 'w')

dsplit = '#' + 54*'-'

#------------------
# headers
f.write('import pickle' + '\n')
f.write('' + '\n')

#------------------
# calibrations
f.write(dsplit + '\n')

# 1 - bias
f.write('bias_obs = [' + '\n')
for obs in red_bias:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

# 2 - domeflat
f.write('domeflat_obs = [' + '\n')
for obs in red_domeflat:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

# 3 - twiflat
f.write('twiflat_obs = [' + '\n')
for obs in red_twiflat:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

# 4 - dark
f.write('dark_obs = [' + '\n')
for obs in red_dark:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

# 5 - arc
f.write('arc_obs = [' + '\n')
for obs in red_arc:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

# 6 - wire
f.write('wire_obs = [' + '\n')
for obs in red_wire:
    f.write('   \'%s\',' % obs + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

#------------------
# science
f.write(dsplit+ '\n')
f.write('sci_obs = [' + '\n')
for obj_name in red_science.keys():
    obs_list = red_science[obj_name]
    obs_str = '\'%s\'' % obs_list[0]
    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    f.write('    # %s' % obj_name + '\n')
    f.write('    {\'sci\'  : [%s],' % obs_str + '\n')
    f.write('     \'sky\'  : []},' + '\n')


f.write('    ]' + '\n')
f.write('' + '\n')

#------------------
# stdstars
f.write(dsplit+ '\n')
f.write('std_obs = [' + '\n')
for obj_name in red_stdstar.keys():
    obs_list = red_stdstar[obj_name]
    obs_str = '\'%s\'' % obs_list[0]
    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    f.write('    # %s' % obj_name + '\n')
    f.write('    {\'sci\'  : [%s],' % obs_str + '\n')
    f.write('     \'name\' : [\'%s\'],' % obj_name + '\n')
    f.write('     \'type\' : [\'flux\', \'telluric\']},' + '\n')
f.write('    ]' + '\n')
f.write('' + '\n')

#------------------
# footers
f.write(dsplit+ '\n')
out_fn = 'wifesR_%s_metadata.pkl' % obs_date
f.write('night_data = {' + '\n')
f.write('    \'bias\' : bias_obs,' + '\n')
f.write('    \'domeflat\' : domeflat_obs,' + '\n')
f.write('    \'twiflat\' : twiflat_obs,' + '\n')
f.write('    \'dark\' : dark_obs,' + '\n')
f.write('    \'wire\' : wire_obs,' + '\n')
f.write('    \'arc\'  : arc_obs,' + '\n')
f.write('    \'sci\'  : sci_obs,' + '\n')
f.write('    \'std\'  : std_obs}' + '\n')
f.write('' + '\n')
f.write('f1 = open(\'%s\', \'wb\')' % out_fn + '\n')
f.write('pickle.dump(night_data, f1)' + '\n')
f.write('f1.close()' + '\n')

f.close()
