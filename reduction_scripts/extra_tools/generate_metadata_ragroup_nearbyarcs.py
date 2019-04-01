#! /usr/bin/env python

import sys
import os
try:
	import pyfits
except:
	import astropy.io.fits as pyfits
import wifes_calib
import numpy

class MyList(object):
    def __init__(self):
        self.list = ["dummy"]
    def __iter__(self):
        return iter(self.list)
  #  def __append__(self):
  #      return append(self.list)

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
        obs_date = obs[7:15]
    if obs[:6] == 'T2m3wb':
        if obs in blue_obs:
            continue
        else:
            blue_obs.append(obs)
    if obs[:6] == 'T2m3wr':
        if obs in red_obs:
            continue
        else:
            red_obs.append(obs)


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
blue_stdstar = []
blue_science = {}
#ra = MyList()
#dec = MyList()
ra = []
dec = []
count=0
for obs in blue_obs:
    fn = data_dir+obs+'.fits'
    f = pyfits.open(fn)
    imagetype = f[0].header['IMAGETYP']
    object = f[0].header['OBJECT']
    f.close()
    # 1 - bias frames
    if imagetype == 'zero':
        blue_bias.append(obs)
    # 2 - quartz flats
    if imagetype == 'flat':
        blue_domeflat.append(obs)
    # 3 - twilight flats
    if imagetype == 'skyflat':
        blue_twiflat.append(obs)
    # 4 - dark frames
    if imagetype == 'dark':
        blue_dark.append(obs)
    # 5 - arc frames
    if imagetype == 'arc':
        blue_arc.append(obs)
    # 6 - wire frames
    if imagetype == 'wire':
        blue_wire.append(obs)
    # all else are science targets
    if imagetype == 'object':
        if object in stdstar_list:
            blue_stdstar.append(obs)
        else:
            thisdec = f[0].header['DEC']
            thisra  = f[0].header['RA']
            thisdec = thisdec[0:9]
            thisra = thisra[0:8]
            object = thisra+' '+thisdec
            if thisra in ra and thisdec in dec: ##if ra and dec verbatim already seen
                print(object)
                # TODO: fix this
                try:
                    blue_science[object].append(obs)
                except:
                    blue_science[object]= [obs]
            else:
                blue_science[object] = [obs]
            ra.append(thisra)
            dec. append(thisdec)

#------------------------------------------------------
# write to metadata save script!
f = open('save_blue_metadata.py', 'w')

dsplit = '#' + 54*'-'

#------------------
# headers
print('import pickle', file=f)
print('', file=f)

#------------------
# calibrations
print(dsplit, file=f)

# 1 - bias
print('bias_obs = [', file=f)
for obs in blue_bias:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

# 2 - domeflat
print('domeflat_obs = [', file=f)
for obs in blue_domeflat:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

# 3 - twiflat
print('twiflat_obs = [', file=f)
for obs in blue_twiflat:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

# 4 - dark
print('dark_obs = [', file=f)
for obs in blue_dark:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

# 5 - arc
print('arc_obs = [', file=f)
for obs in blue_arc:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

# 6 - wire
print('wire_obs = [', file=f)
for obs in blue_wire:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

#------------------
# science
print(dsplit, file=f)
print('sci_obs = [', file=f)
for object in blue_science.keys():
    curbarc = 1000
    curaarc = 1000
    
    obs_list = blue_science[object]
    obs_str = '\'%s\'' % obs_list[0]
    obsnumber = int(obs_str[24:28])
    for runarc in blue_arc:
        thisarc_str = '\'%s\'' % runarc
        thisarcnum = int(thisarc_str[24:28])
        if thisarcnum < obsnumber:
            thisarcdiff=abs(thisarcnum-obsnumber)
            curarcdiff = abs(curbarc-obsnumber)
            if thisarcdiff < curarcdiff:
                curbarc=thisarcnum
                curbarcstr=thisarc_str
        if thisarcnum > obsnumber:
            thisarcdiff=abs(thisarcnum-obsnumber)
            curarcdiff = abs(curaarc-obsnumber)
            if thisarcdiff < curarcdiff:
                curaarc=thisarcnum
                curaarcstr=thisarc_str
    
    sciarc_str = curbarcstr
    sciarc_str += ',\n               '+curaarcstr

    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    print('    # %s' % object, file=f)
    print('    {\'sci\'  : [%s],' % obs_str, file=f)
    print('     \'arc\'  : ['+sciarc_str+'],', file=f)
    print('     \'flat\' : [],', file=f)
    print('     \'sky\'  : [],', file=f)
    print('     \'wire\' : []},', file=f)


print('    ]', file=f)
print('', file=f)

#------------------
# stdstars
print(dsplit, file=f)
print('std_obs = [', file=f)
for obs in blue_stdstar:
    f2 = pyfits.open(data_dir+obs+'.fits')
    object = f2[0].header['OBJECT']
    f2.close()
    print('    # %s' % object, file=f)
    print('    {\'sci\'  : [\'%s\'],' % obs, file=f)
    print('     \'arc\'  : [],', file=f)
    print('     \'flat\' : [],', file=f)
    print('     \'type\' : [\'flux\', \'telluric\']},', file=f)
print('    ]', file=f)
print('', file=f)

#------------------
# footers
print(dsplit, file=f)
out_fn = 'wifesB_%s_metadata.pkl' % obs_date
print('night_data = {', file=f)
print('    \'bias\' : bias_obs,', file=f)
print('    \'domeflat\' : domeflat_obs,', file=f)
print('    \'twiflat\' : twiflat_obs,', file=f)
print('    \'dark\' : dark_obs,', file=f)
print('    \'wire\' : wire_obs,', file=f)
print('    \'arc\'  : arc_obs,', file=f)
print('    \'sci\'  : sci_obs,', file=f)
print('    \'std\'  : std_obs}', file=f)
print('', file=f)
print('f1 = open(\'%s\', \'w\')' % out_fn, file=f)
print('pickle.dump(night_data, f1)', file=f)
print('f1.close()', file=f)

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
red_stdstar = []
red_science = {}
#ra = MyList()
#dec= MyList()
ra = []
dec = []
#print ra
##ra = ra.append('stuff')
##dec = dec.append('stuff')
for obs in red_obs:
    fn = data_dir+obs+'.fits'
    f = pyfits.open(fn)
    ##print f[0].header
    imagetype = f[0].header['IMAGETYP']
    object = f[0].header['OBJECT']
    f.close()
    # 1 - bias frames
    if imagetype == 'zero':
        red_bias.append(obs)
    # 2 - quartz flats
    if imagetype == 'flat':
        red_domeflat.append(obs)
    # 3 - twilight flats
    if imagetype == 'skyflat':
        red_twiflat.append(obs)
    # 4 - dark frames
    if imagetype == 'dark':
        red_dark.append(obs)
    # 5 - arc frames
    if imagetype == 'arc':
        red_arc.append(obs)
    # 6 - wire frames
    if imagetype == 'wire':
        red_wire.append(obs)
    # all else are science targets
    if imagetype == 'object':
        if object in stdstar_list:
            red_stdstar.append(obs)
        else:
            thisdec = f[0].header['DEC']
            thisra  = f[0].header['RA']
            thisdec = thisdec[0:9]
            thisra = thisra[0:8]
            
            #if object == '':
            ##re-label objects
            object = thisra+' '+thisdec
            if thisra in ra and thisdec in dec: ##if ra and dec verbatim already seen
                # TODO: FIX THIS
                try:
                    red_science[object].append(obs)
                except:
                    red_science[object] = [obs]
            else:
                red_science[object] = [obs]
            
            ##append the current ra and dec to the lists for the next stars
            ra.append(thisra)
            dec.append(thisdec)    
#------------------------------------------------------
# write to metadata save script!
f = open('save_red_metadata.py', 'w')

dsplit = '#' + 54*'-'

#------------------
# headers
print('import pickle', file=f)
print('', file=f)

#------------------
# calibrations
print(dsplit, file=f)

# 1 - bias
print('bias_obs = [', file=f)
for obs in red_bias:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

# 2 - domeflat
print('domeflat_obs = [', file=f)
for obs in red_domeflat:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

# 3 - twiflat
print('twiflat_obs = [', file=f)
for obs in red_twiflat:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

# 4 - dark
print('dark_obs = [', file=f)
for obs in red_dark:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

# 5 - arc
print('arc_obs = [', file=f)
for obs in red_arc:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

# 6 - wire
print('wire_obs = [', file=f)
for obs in red_wire:
    print('   \'%s\',' % obs, file=f)
print('    ]', file=f)
print('', file=f)

#------------------
# science
print(dsplit, file=f)
print('sci_obs = [', file=f)
for object in red_science.keys():
    curbarc = 1000
    curaarc = 1000
    
    obs_list = red_science[object]
    obs_str = '\'%s\'' % obs_list[0]
    obsnumber = int(obs_str[24:28])
    for runarc in red_arc:
        thisarc_str = '\'%s\'' % runarc
        thisarcnum = int(thisarc_str[24:28])
        if thisarcnum < obsnumber:
            thisarcdiff=abs(thisarcnum-obsnumber)
            curarcdiff = abs(curbarc-obsnumber)
            if thisarcdiff < curarcdiff:
                curbarc=thisarcnum
                curbarcstr=thisarc_str
        if thisarcnum > obsnumber:
            thisarcdiff=abs(thisarcnum-obsnumber)
            curarcdiff = abs(curaarc-obsnumber)
            if thisarcdiff < curarcdiff:
                curaarc=thisarcnum
                curaarcstr=thisarc_str
    
    sciarc_str = curbarcstr
    sciarc_str += ',\n               '+curaarcstr
    obs_list = red_science[object]
    obs_str = '\'%s\'' % obs_list[0]
    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    print('    # %s' % object, file=f)
    print('    {\'sci\'  : [%s],' % obs_str, file=f)
    print('     \'arc\'  : ['+sciarc_str+'],', file=f)
    print('     \'flat\' : [],', file=f)
    print('     \'sky\'  : [],', file=f)
    print('     \'wire\' : []},', file=f)


print('    ]', file=f)
print('', file=f)

#------------------
# stdstars
print(dsplit, file=f)
print('std_obs = [', file=f)
for obs in red_stdstar:
    f2 = pyfits.open(data_dir+obs+'.fits')
    object = f2[0].header['OBJECT']
    f2.close()
    print('    # %s' % object, file=f)
    print('    {\'sci\'  : [\'%s\'],' % obs, file=f)
    print('     \'arc\'  : [],', file=f)
    print('     \'flat\' : [],', file=f)
    print('     \'type\' : [\'flux\', \'telluric\']},', file=f)
print('    ]', file=f)
print('', file=f)

#------------------
# footers
print(dsplit, file=f)
out_fn = 'wifesR_%s_metadata.pkl' % obs_date
print('night_data = {', file=f)
print('    \'bias\' : bias_obs,', file=f)
print('    \'domeflat\' : domeflat_obs,', file=f)
print('    \'twiflat\' : twiflat_obs,', file=f)
print('    \'dark\' : dark_obs,', file=f)
print('    \'wire\' : wire_obs,', file=f)
print('    \'arc\'  : arc_obs,', file=f)
print('    \'sci\'  : sci_obs,', file=f)
print('    \'std\'  : std_obs}', file=f)
print('', file=f)
print('f1 = open(\'%s\', \'w\')' % out_fn, file=f)
print('pickle.dump(night_data, f1)', file=f)
print('f1.close()', file=f)

f.close()
