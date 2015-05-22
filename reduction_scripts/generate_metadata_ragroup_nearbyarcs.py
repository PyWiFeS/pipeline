#! /usr/bin/env python

import sys
import os
import pyfits
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
print blue_obs
print red_obs
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
blue_stdstar = []
blue_science = {}
#ra = MyList()
#dec = MyList()
ra = []
dec = []
count=0
for obs in blue_obs:
    print obs
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
           # print thisra,thisdec
           # stop
            thisdec = thisdec[0:9]
            thisra = thisra[0:8]
            object = thisra+' '+thisdec
            #print object
            if thisra in ra and thisdec in dec: ##if ra and dec verbatim already seen
                blue_science[object].append(obs)
               # print blue_science[object]
               # stop
            else:
               #stop
                blue_science[object] = [obs]
            ##append the current ra and dec to the lists for the next stars
            #ra  = [ra,thisra]
            ra.append(thisra)
            dec. append(thisdec)

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
for object in blue_science.keys():
    curbarc = 1000
    curaarc = 1000
    
    obs_list = blue_science[object]
    obs_str = '\'%s\'' % obs_list[0]
    print 'dude'
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
    
    #print curaarcstr
    #print curbarcstr
    sciarc_str = curbarcstr
    sciarc_str += ',\n               '+curaarcstr
    #stop

    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    print >> f, '    # %s' % object
    print >> f, '    {\'sci\'  : [%s],' % obs_str
    print >> f, '     \'arc\'  : ['+sciarc_str+'],'
    print >> f, '     \'flat\' : [],'
    print >> f, '     \'sky\'  : [],'
    print >> f, '     \'wire\' : []},'


print >> f, '    ]'
print >> f, ''

#------------------
# stdstars
print >> f, dsplit
print >> f, 'std_obs = ['
for obs in blue_stdstar:
    f2 = pyfits.open(data_dir+obs+'.fits')
    object = f2[0].header['OBJECT']
    f2.close()
    print >> f, '    # %s' % object
    print >> f, '    {\'sci\'  : [\'%s\'],' % obs
    print >> f, '     \'arc\'  : [],'
    print >> f, '     \'flat\' : [],'
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
            print object
            if thisra in ra and thisdec in dec: ##if ra and dec verbatim already seen
                red_science[object].append(obs)
            else:
                red_science[object] = [obs]
            
           ## print red_science[object]
            ##append the current ra and dec to the lists for the next stars
            #ra  = [ra,thisra]
            #dec = [dec,thisdec]
            ra.append(thisra)
            dec.append(thisdec)    
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
for object in red_science.keys():
    curbarc = 1000
    curaarc = 1000
    
    obs_list = red_science[object]
    obs_str = '\'%s\'' % obs_list[0]
    ##print 'dude'
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
    
    #print curaarcstr
    #print curbarcstr
    sciarc_str = curbarcstr
    sciarc_str += ',\n               '+curaarcstr#print object
    obs_list = red_science[object]
    #print obs_list
    obs_str = '\'%s\'' % obs_list[0]
    for i in range(1,len(obs_list)):
        obs = obs_list[i]
        obs_str += ',\n               \'%s\'' % obs
    print >> f, '    # %s' % object
    print >> f, '    {\'sci\'  : [%s],' % obs_str
    print >> f, '     \'arc\'  : ['+sciarc_str+'],'
    print >> f, '     \'flat\' : [],'
    print >> f, '     \'sky\'  : [],'
    print >> f, '     \'wire\' : []},'


print >> f, '    ]'
print >> f, ''

#------------------
# stdstars
print >> f, dsplit
print >> f, 'std_obs = ['
for obs in red_stdstar:
    f2 = pyfits.open(data_dir+obs+'.fits')
    object = f2[0].header['OBJECT']
    f2.close()
    print >> f, '    # %s' % object
    print >> f, '    {\'sci\'  : [\'%s\'],' % obs
    print >> f, '     \'arc\'  : [],'
    print >> f, '     \'flat\' : [],'
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
