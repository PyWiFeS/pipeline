#! /usr/bin/env python

import sys
import os
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

dsplit = '#' + 54*'-' + '\n'

#------------------
# headers
f.write('import pickle\n')
f.write('\n')

#------------------
# calibrations
f.write(dsplit)

# 1 - bias
f.write('bias_obs = [\n')
for obs in blue_bias:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

# 2 - domeflat
f.write('domeflat_obs = [\n')
for obs in blue_domeflat:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

# 3 - twiflat
f.write('twiflat_obs = [\n')
for obs in blue_twiflat:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

# 4 - dark
f.write('dark_obs = [\n')
for obs in blue_dark:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

# 5 - arc
f.write('arc_obs = [\n')
for obs in blue_arc:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

# 6 - wire
f.write('wire_obs = [\n')
for obs in blue_wire:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

#------------------
# science
f.write(dsplit)
f.write('sci_obs = [\n')
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
    f.write('    # %s \n' % object)
    f.write('    {\'sci\'  : [%s], \n' % obs_str)
    f.write('     \'arc\'  : ['+sciarc_str+'],\n')
    f.write('     \'flat\' : [],\n')
    f.write('     \'sky\'  : [],\n')
    f.write('     \'wire\' : []},\n')


f.write('    ]\n')
f.write('\n')

#------------------
# stdstars
f.write(dsplit)
f.write('std_obs = [\n')
for obs in blue_stdstar:
    f2 = pyfits.open(data_dir+obs+'.fits')
    object = f2[0].header['OBJECT']
    f2.close()
    f.write('    # %s \n' % object)
    f.write('    {\'sci\'  : [\'%s\'], \n' % obs)
    f.write('     \'arc\'  : [],\n')
    f.write('     \'flat\' : [],\n')
    f.write('     \'type\' : [\'flux\', \'telluric\']},\n')
f.write('    ]\n')
f.write('\n')

#------------------
# footers
f.write(dsplit)
out_fn = 'wifesB_%s_metadata.pkl' % obs_date
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
f.write('f1 = open(\'%s\', \'wb\') \n' % out_fn)
f.write('pickle.dump(night_data, f1)\n')
f.write('f1.close()\n')

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

dsplit = '#' + 54*'-' + '\n'

#------------------
# headers
f.write('import pickle\n')
f.write('\n')

#------------------
# calibrations
f.write(dsplit)

# 1 - bias
f.write('bias_obs = [\n')
for obs in red_bias:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

# 2 - domeflat
f.write('domeflat_obs = [\n')
for obs in red_domeflat:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

# 3 - twiflat
f.write('twiflat_obs = [\n')
for obs in red_twiflat:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

# 4 - dark
f.write('dark_obs = [\n')
for obs in red_dark:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

# 5 - arc
f.write('arc_obs = [\n')
for obs in red_arc:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

# 6 - wire
f.write('wire_obs = [\n')
for obs in red_wire:
    f.write('   \'%s\', \n' % obs)
f.write('    ]\n')
f.write('\n')

#------------------
# science
f.write(dsplit)
f.write('sci_obs = [\n')
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
    f.write('    # %s \n' % object)
    f.write('    {\'sci\'  : [%s], \n' % obs_str)
    f.write('     \'arc\'  : ['+sciarc_str+'],\n')
    f.write('     \'flat\' : [],\n')
    f.write('     \'sky\'  : [],\n')
    f.write('     \'wire\' : []},\n')


f.write('    ]\n')
f.write('\n')

#------------------
# stdstars
f.write(dsplit)
f.write('std_obs = [\n')
for obs in red_stdstar:
    f2 = pyfits.open(data_dir+obs+'.fits')
    object = f2[0].header['OBJECT']
    f2.close()
    f.write('    # %s' % object)
    f.write('    {\'sci\'  : [\'%s\'],' % obs)
    f.write('     \'arc\'  : [],\n')
    f.write('     \'flat\' : [],\n')
    f.write('     \'type\' : [\'flux\', \'telluric\']},\n')
f.write('    ]\n')
f.write('\n')

#------------------
# footers
f.write(dsplit)
out_fn = 'wifesR_%s_metadata.pkl' % obs_date
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
f.write('f1 = open(\'%s\', \'wb\') \n' % out_fn)
f.write('pickle.dump(night_data, f1)\n')
f.write('f1.close()\n')

f.close()
