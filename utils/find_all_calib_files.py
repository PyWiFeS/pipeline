#! /usr/bin/env python

import numpy as np
import sys
import os
import glob
from astropy.io import fits as pyfits
from collections import OrderedDict
from astropy.table import Table
import pywifes_utils as pu

"""
Often you don't take all calib frames in the same night. Make a folder just with all calib frames from the run. Then select ones with the closest data (preferrably the same date as science frames).
"""

# Input: wifes folder with all the nights (but no specific night)
root = '/data/mash/marusa/2m3data/wifes/'

# Output
output_root = '/data/mash/marusa/2m3reduced/wifes/'
output_filename = 'calibration_all_filenames.py'
output_filename = os.path.join(output_root, output_filename)
# Prepare list of files
all_files=[]

i=0
for path, subdirs, files in os.walk(root):
    if path!=root:
        for f in files:
            print path, f
            if 'T2m3ag' in f or 'T2M3Ec' in f or not f.endswith('.fits') or 'wsol' in f or 'wire_soln' in f or 'wave_soln' in f:
                continue
            else:
                filename=os.path.join(path, f)
                all_files.append(filename)
                i=1
                
                #~ if i>100:
                    #~ break
            #~ if i>100:
                #~ break
        #~ if i>100:
            #~ break

try:
    print('Excluding bad frames...')
    all_files = pu.exclude_bad_frames(all_files, '/data/mash/marusa/2m3data/wifes/list_of_bad_exposures_that_we_shouldn_use.dat') # TODO: HARDCODED
except:
    'No BAD files filename provided.'

print('Ready to start with %d frames.'%len(all_files))

keywords=['IMAGETYP', 'GRATINGB', 'GRATINGR', 'BEAMSPLT', 'CCDSUM', 'CCDSEC']
keywords_dark_zero = ['IMAGETYP', 'CCDSEC', 'CCDSUM']


def darks_and_zeros():
    """
    The only thing that matters for darks are CCDSUM and WINDOW. (? I'm guessing here.)
    """        

    keywords=keywords_dark_zero

    result=dict()

    for fn in all_files:
        try:
            f = pyfits.open(fn)
            header = f[0].header
            f.close()
            
            date = fn.split('/')
            date=date[-2]

        except:
            continue
        
        # Take only darks
        try:
            imagetyp=header['IMAGETYP'].upper()
        except:
            print('%s no imagetyp'%fn)
            continue
        if header['IMAGETYP'].upper() not in ['DARK', 'ZERO'] :
            continue

        # Get header info
        try:
            k=[header[x] for x in keywords]
            k[0]=k[0].upper()
        except:
            #~ print fn, header['IMAGETYP'].upper()
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
    
    
    return result
    
        
def prepare_result():   
    result=dict()

    for fn in all_files:
        try:
            f = pyfits.open(fn)
            header = f[0].header
            f.close()
        except:
            continue
        
        try:
            imagetyp=header['IMAGETYP'].upper()
        except:
            print('%s no imagetyp'%fn)
            continue
        
        # Objects and arcs: continue
        #~ if header['IMAGETYP'].upper()=='OBJECT' or header['IMAGETYP'].upper()=='ARC': # Assuming that all the rest is what we seek for here.
        if header['IMAGETYP'].upper()=='OBJECT' or header['IMAGETYP'].upper()=='ARC' or header['IMAGETYP'].upper()=='DARK' or header['IMAGETYP'].upper()=='ZERO': # Assuming that all the rest is what we seek for here.
            continue

        # Get header info
        try:
            k=[header[x] for x in keywords]
            k[0]=k[0].upper()
        except:
            print('Couldnt read header', fn, header['IMAGETYP'].upper())
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

    #~ print('result')
    #~ print(result)
    darks_and_zeros_dict = darks_and_zeros()
    
    for k, v in darks_and_zeros_dict.iteritems():
        result[k]=v
    
    #~ result = result.update(darks_and_zeros_dict)
    
    #~ print('result')
    #~ print(result)

    # Sort by the numb... TODO it doesnt work yet
    kys = result.keys()
    print kys
    #~ print(result[kys[0]])
    #~ print(result[kys[0]].itervalues())
    #~ print('printing')
    #~ for k in kys:
        #~ for imgtdict in result[k].itervalues():
            #~ print(imgtdict, len(imgtdict))
            #~ print()
    
    
    
    lngths = [np.sum([np.sum([len(vv) for vv in imgtlist]) for imgtlist in result[k]]) for k in kys]
    #~ print('lengths', lngths)
    sortidx = np.argsort(lngths)
    sortidx=sortidx[::-1] # reverse the order
    #~ print sortidx
    kys2 = [kys[x] for x in sortidx]
    kys=kys2
    #~ result = OrderedDict(sorted(result.viewitems(), key=lambda x: len(x[1]), reverse=True))

    # WRITE
    f=open(output_filename, 'wb')
    dsplit = '#' + 54*'-' + '\n'
    
    f.write('result = {\n')

    # Split by date
    #~ for k, v in result.iteritems(): # for each mode
    for k in kys: # for each mode
        v=result[k]
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
        
    print(output_filename, 'written.')
        


prepare_result()
#~ darks_and_zeros()
