#! /usr/bin/env python

import numpy as np
import os
import glob
from astropy.io import fits as pyfits


# get list of all fits files in directory

data_dir='../20190303/raw'
all_files = os.listdir(data_dir)
all_files = [os.path.join(data_dir, x) for x in all_files]


def find_all_modes():
    """
    Find all different modes taken during this night.
    By modes I mean different combinations of gratings, binning etc.
    Find images that go together.
    
    NOTE that some images were taken with echelle. --> Actually not. ECHELLE DATA ARE MISSING!?!
    """

    modes=dict()
    
    keywords=['NAXIS1', 'NAXIS2', 'WINDOW', 'GRATINGB', 'GRATINGR', 'BEAMSPLT', 'CCDSIZE', 'CCDSEC', 'CCDSUM', 'TRIMSEC', 'DATASEC', 'DETSEC']
    
    for fn in all_files:
        try:
            f = pyfits.open(fn)
            header = f[0].header
            f.close()
        except:
            continue

        #~ print fn
        #~ print header
        k=[header[x] for x in keywords]
        #~ print k
        
        k=tuple(k) # Lists or sets cannot be dictionary keys
        try:
            modes[k].append(fn)
        except:
            modes[k]=[fn]
            
        #~ print
        #~ print
        
        
        
    for k, v in modes.iteritems():
        print k
        print v
        print
        
    for k, v in modes.iteritems():
        print k

find_all_modes()
