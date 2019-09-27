#! /usr/bin/env python

import numpy as np
import sys
import os
import glob
from astropy.io import fits as pyfits
from collections import OrderedDict

"""
Often you don't take all calib frames in the same night. Make a folder just with all calib frames from the run. Then select ones with the closest data (preferrably the same date as science frames).
"""

root = '/priv/mulga1/marusa/2m3data/wifes/'

# Prepare list of files
all_files=[]
for path, subdirs, files in os.walk(root):
    if path!=root:
        if '2019' in path: # Only echelle data in 2018
            for f in files:
                if 'ag' in f:
                    continue
                else:
                    filename=os.path.join(path, f)
                    all_files.append(filename)

keywords=['IMAGETYP', 'NAXIS1', 'NAXIS2', 'WINDOW', 'GRATINGB', 'GRATINGR', 'BEAMSPLT', 'CCDSIZE', 'CCDSEC', 'CCDSUM', 'TRIMSEC', 'DATASEC', 'DETSEC']
# ('DARK', 4202, 1028, 'REG_1x1_4202x2056+0+2056', 'B3000', 'R7000', 'RT480', '[1:4202,1:4112]', '[1:4202,2057:4112]', '1 2', '[1:4202,1:1028]', '[1:4202,1:1028]', '[1:4202,2057:4112]')


#~ keywords=['IMAGETYP', 'NAXIS1', 'NAXIS2', 'WINDOW', 'GRATINGB', 'GRATINGR', 'CCDSEC', 'CCDSUM', 'TRIMSEC', 'DATASEC', 'DETSEC']



def testing():
    modes=dict()
    modeset=dict()

    modesobjects=dict()
    modesetobjects=dict()

    for fn in all_files:
        try:
            f = pyfits.open(fn)
            header = f[0].header
            f.close()
        except:
            continue
        
        if header['IMAGETYP'].upper()=='OBJECT' or header['IMAGETYP'].upper()=='ARC': # Assuming that all the rest is what we seek for here.
            if header['IMAGETYP'].upper()=='OBJECT':
                try:
                    k=[header[x] for x in keywords]
                    k[0]=k[0].upper()
                except:
                    print fn, header['IMAGETYP'].upper()
                    continue
                
                k=tuple(k) # Lists or sets cannot be dictionary keys
                try:
                    modesobjects[k].append(fn)
                    #~ modeset[k[1:]].append(fn)
                    modesetobjects[k[1:]].append(k[0])
                except:
                    modesobjects[k]=[fn]
                    #~ modeset[k[1:]]=[fn]
                    modesetobjects[k[1:]]=[k[0]]

            continue

        try:
            k=[header[x] for x in keywords]
            k[0]=k[0].upper()
        except:
            print fn, header['IMAGETYP'].upper()
            continue
        
        k=tuple(k) # Lists or sets cannot be dictionary keys
        try:
            modes[k].append(fn)
            #~ modeset[k[1:]].append(fn)
            modeset[k[1:]].append(k[0])
        except:
            modes[k]=[fn]
            #~ modeset[k[1:]]=[fn]
            modeset[k[1:]]=[k[0]]

    for k, v in modes.iteritems():
        print k
        print v
        print


    print    
    print    
    print    
    print    
    for k, v in modeset.iteritems():
        print k, len(v)
        print set(v)
        print
            
    print len(modeset)

    print    
    print    
    print    
    print    
    for k, v in modesetobjects.iteritems():
        print k, len(v)
        print set(v)
        try:
            cal=modeset[k]
            l=len(cal)
            t=set(cal)
            
            # frequency
            
            print len(cal), set(cal)
        except:
            pass

        print
            
    print len(modesetobjects)


def prepare_pickle():   
    result=dict()

    for fn in all_files:
        try:
            f = pyfits.open(fn)
            header = f[0].header
            f.close()
        except:
            continue
        
        # Objects and arcs: continue
        if header['IMAGETYP'].upper()=='OBJECT' or header['IMAGETYP'].upper()=='ARC': # Assuming that all the rest is what we seek for here.
            continue

        # Get header info
        try:
            k=[header[x] for x in keywords]
            k[0]=k[0].upper()
        except:
            print fn, header['IMAGETYP'].upper()
            continue
        
        k=tuple(k) # Lists or sets cannot be dictionary keys
        
        try:
            d=result[k[1:]]
            
            try:
                d[k[0]].append(fn)
            except:
                d[k[0]]=[fn]
            
            result[k[1:]]=d

        except:
            result[k[1:]]={k[0]: [fn]}


    # Split by date
    for k, v in result.iteritems(): # for each mode
        
        c2=dict()
        for kk, vv in c.iteritems(): # for each imagetype
            for vvv in vv:
                date=vvv.split('/')[-1].split('.')[0].split('-')[1]
                
                try:
                    tmp=c2[kk]
                    try:
                        tmp[date].append(vvv)
                    except:
                        tmp[date]=[vvv]
                    c2[kk]=tmp
                except:
                    c2[kk]={date: [vvv]}

        for x, y in c2.iteritems():
            print x
            for yy, yyy in y.iteritems():
                print yy, len(yyy)





    # Print results
    result = OrderedDict(sorted(result.viewitems(), key=lambda x: len(x[1]), reverse=True))

    f=open('calibration_filenames.py', 'wb')
    dsplit = '#' + 54*'-' + '\n'

    f.write('result = {\n')
    for k, v in result.iteritems():
        line='\n     '+str(k)+': {\n'
        f.write(line)
        for t, names in v.iteritems():
            names=sorted(names)
            f.write('         "%s": ["%s",\n'%(t, names[0]))
            for n in names[1:-1]:
                f.write('             "%s",\n'%n)
            f.write('             "%s"],\n'%names[-1])
            #~ f.write('           ],\n')
        f.write('           },\n')
    f.write('}')
        

            
        
        
    f.close()
        
        
        
def prepare_result():   
    result=dict()

    for fn in all_files:
        try:
            f = pyfits.open(fn)
            header = f[0].header
            f.close()
        except:
            continue
        
        # Objects and arcs: continue
        if header['IMAGETYP'].upper()=='OBJECT' or header['IMAGETYP'].upper()=='ARC': # Assuming that all the rest is what we seek for here.
            continue

        # Get header info
        try:
            k=[header[x] for x in keywords]
            k[0]=k[0].upper()
        except:
            print fn, header['IMAGETYP'].upper()
            continue
        
        k=tuple(k) # Lists or sets cannot be dictionary keys
        
        try:
            d=result[k[1:]]
            
            try:
                d[k[0]].append(fn)
            except:
                d[k[0]]=[fn]
            
            result[k[1:]]=d

        except:
            result[k[1:]]={k[0]: [fn]}


    result = OrderedDict(sorted(result.viewitems(), key=lambda x: len(x[1]), reverse=True))

    f=open('calibration_filenames_date.py', 'wb')
    #~ f=open('calibration_filenames_date_less_keywords.py', 'wb')
    dsplit = '#' + 54*'-' + '\n'
    
    f.write('result = {\n')

    # Split by date
    for k, v in result.iteritems(): # for each mode
        line='\n     '+str(k)+': {\n'
        f.write(line) # mode
                
        c2=dict()
        for kk, vv in v.iteritems(): # for each imagetype
            for vvv in vv: # For each filename
                folder=vvv.split('/')[-2]
                date=vvv.split('/')[-1].split('.')[0].split('-')[1]
                if folder==date:
                    pass
                else:
                    continue # do not repeat files that have been copied there for easier calibration
                
                try:
                    tmp=c2[kk]
                    try:
                        tmp[date].append(vvv)
                    except:
                        tmp[date]=[vvv]
                    c2[kk]=tmp
                except:
                    c2[kk]={date: [vvv]}

        
        # Print in the file
        
        for t, y in c2.iteritems(): # c2[imagetype] = {date: [filenames]}
            f.write('         "%s": {\n'%(t)) # imagetype
            for yy, yyy in y.iteritems(): # date: [filenames]
                yyy=sorted(yyy)
                #~ f.write('           %s: [\n'%yy) # date
                f.write('           %s: ["%s",\n'%(yy, yyy[0])) # date
                for name in yyy[1:-1]:
                    f.write('               "%s",\n'%name)
                f.write('               "%s"],\n'%yyy[-1])
                #~ f.write('           ],\n')
            f.write('         },\n') # imagetype
        f.write('         },\n') # mode
        f.write('##################################################################\n')
    f.write('}')


 
        
        
    f.close()
        
        
        


prepare_result()
