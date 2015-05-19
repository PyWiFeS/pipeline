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
        f.close()
    except:
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
    if imagetype == 'ZERO':
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
print >> f, 'import pickle'
print >> f, ''

#------------------
# calibrations
print >> f, dsplit

# 1 - bias
print >> f, 'bias_obs = ['
for obs in blue_bias:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

# 2 - domeflat
print >> f, 'domeflat_obs = ['
for obs in blue_domeflat:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

# 3 - twiflat
print >> f, 'twiflat_obs = ['
for obs in blue_twiflat:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

# 4 - dark
print >> f, 'dark_obs = ['
for obs in blue_dark:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

# 5 - arc
print >> f, 'arc_obs = ['
for obs in blue_arc:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

# 6 - wire
print >> f, 'wire_obs = ['
for obs in blue_wire:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

#------------------
# science
print >> f, dsplit
print >> f, 'sci_obs = ['
for obj_name in blue_science.keys():
    obs_list = blue_science[obj_name]
    obs_str = '\'%s\'' % obs_list[0]
    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    print >> f, '    # %s' % obj_name
    print >> f, '    {\'sci\'  : [%s],' % obs_str
    print >> f, '     \'sky\'  : []},'


print >> f, '    ]'
print >> f, ''

#------------------
# stdstars
print >> f, dsplit
print >> f, 'std_obs = ['
for obj_name in blue_stdstar.keys():
    obs_list = blue_stdstar[obj_name]
    obs_str = '\'%s\'' % obs_list[0]
    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    print >> f, '    # %s' % obj_name
    print >> f, '    {\'sci\'  : [%s],' % obs_str
    print >> f, '     \'name\' : [\'%s\'],' % obj_name
    print >> f, '     \'type\' : [\'flux\', \'telluric\']},'
print >> f, '    ]'
print >> f, ''

#------------------
# footers
print >> f, dsplit
out_fn = 'wifesB_%s_metadata.pkl' % obs_date
print >> f, 'night_data = {'
print >> f, '    \'bias\' : bias_obs,'
print >> f, '    \'domeflat\' : domeflat_obs,'
print >> f, '    \'twiflat\' : twiflat_obs,'
print >> f, '    \'dark\' : dark_obs,'
print >> f, '    \'wire\' : wire_obs,'
print >> f, '    \'arc\'  : arc_obs,'
print >> f, '    \'sci\'  : sci_obs,'
print >> f, '    \'std\'  : std_obs}'
print >> f, ''
print >> f, 'f1 = open(\'%s\', \'w\')' % out_fn
print >> f, 'pickle.dump(night_data, f1)'
print >> f, 'f1.close()'

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
    if imagetype == 'ZERO':
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
print >> f, 'import pickle'
print >> f, ''

#------------------
# calibrations
print >> f, dsplit

# 1 - bias
print >> f, 'bias_obs = ['
for obs in red_bias:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

# 2 - domeflat
print >> f, 'domeflat_obs = ['
for obs in red_domeflat:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

# 3 - twiflat
print >> f, 'twiflat_obs = ['
for obs in red_twiflat:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

# 4 - dark
print >> f, 'dark_obs = ['
for obs in red_dark:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

# 5 - arc
print >> f, 'arc_obs = ['
for obs in red_arc:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

# 6 - wire
print >> f, 'wire_obs = ['
for obs in red_wire:
    print >> f, '   \'%s\',' % obs
print >> f, '    ]'
print >> f, ''

#------------------
# science
print >> f, dsplit
print >> f, 'sci_obs = ['
for obj_name in red_science.keys():
    obs_list = red_science[obj_name]
    obs_str = '\'%s\'' % obs_list[0]
    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    print >> f, '    # %s' % obj_name
    print >> f, '    {\'sci\'  : [%s],' % obs_str
    print >> f, '     \'sky\'  : []},'


print >> f, '    ]'
print >> f, ''

#------------------
# stdstars
print >> f, dsplit
print >> f, 'std_obs = ['
for obj_name in red_stdstar.keys():
    obs_list = red_stdstar[obj_name]
    obs_str = '\'%s\'' % obs_list[0]
    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    print >> f, '    # %s' % obj_name
    print >> f, '    {\'sci\'  : [%s],' % obs_str
    print >> f, '     \'name\' : [\'%s\'],' % obj_name
    print >> f, '     \'type\' : [\'flux\', \'telluric\']},'
print >> f, '    ]'
print >> f, ''

#------------------
# footers
print >> f, dsplit
out_fn = 'wifesR_%s_metadata.pkl' % obs_date
print >> f, 'night_data = {'
print >> f, '    \'bias\' : bias_obs,'
print >> f, '    \'domeflat\' : domeflat_obs,'
print >> f, '    \'twiflat\' : twiflat_obs,'
print >> f, '    \'dark\' : dark_obs,'
print >> f, '    \'wire\' : wire_obs,'
print >> f, '    \'arc\'  : arc_obs,'
print >> f, '    \'sci\'  : sci_obs,'
print >> f, '    \'std\'  : std_obs}'
print >> f, ''
print >> f, 'f1 = open(\'%s\', \'w\')' % out_fn
print >> f, 'pickle.dump(night_data, f1)'
print >> f, 'f1.close()'

f.close()
