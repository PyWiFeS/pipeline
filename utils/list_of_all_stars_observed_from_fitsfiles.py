"""
Search all subfolders of a folder and make a list of all object ids from headers.
"""

import numpy as np
from astropy.io import fits
import os
from astropy.table import Table

root = '/priv/mulga1/marusa/2m3data/'
young_stars_input_catalog_filename = 'young_stars_input_catalog.dat'

keywords = ['FILENAME', 'EXPTIME', 'INSTRUME', 'DATE-OBS', 'PROPID', 'OBJNAME']

t=[]


filenames = []
exptimes = []
instrumes = []
dateobss = []
propids = []
objnames = []
dates = []
detectors = []
airmass = []
ras = []
decs = []

for r, d, f in os.walk(root):
    for filename in f:
        if filename.endswith('.fits') and 'T2m3ag' not in filename:
            filename=os.path.join(r, filename)
            try:
                imagetyp = fits.getheader(filename, 0)['IMAGETYP']
                exptime = fits.getheader(filename, 0)['EXPTIME']

                RA = fits.getheader(filename, 0)['RA']
                if imagetyp=='object' and exptime>0:
                    
                    v = fits.getheader(filename, 0)['FILENAME']
                    filenames.append(v)
                    
                    v = fits.getheader(filename, 0)['INSTRUME']
                    instrumes.append(v)
                    
                    v = fits.getheader(filename, 0)['DATE-OBS']
                    dateobss.append(v)
                    
                    date = v.split('T')[0].replace('-', '')
                    dates.append(date)
                    
                    v = fits.getheader(filename, 0)['PROPID']
                    propids.append(int(v))
                    
                    v = fits.getheader(filename, 0)['OBJNAME']
                    objnames.append(v)
                    
                    v = fits.getheader(filename, 0)['EXPTIME']
                    exptimes.append(v)
                    
                    
                    v = fits.getheader(filename, 0)['DETECTOR']
                    detectors.append(v)
                    
                    
                    v = fits.getheader(filename, 0)['AIRMASS']
                    airmass.append(v)
                    
                    ras.append(RA)

                    v = fits.getheader(filename, 0)['DEC']
                    decs.append(v)
                    
                    
                    tmp=[]
                    for k in keywords:
                        v = fits.getheader(filename, 0)[k]
                        print v
                        tmp.append(v)
                    t.append(tmp)
                    print
            except:
                pass             




tab = Table([filenames], names=['filename'])
tab['OBJNAME'] = objnames
tab['RA'] = ras
tab['DEC'] = decs
tab['INSTRUME'] = instrumes
tab['DETECTOR'] = detectors
tab['DATE-OBS'] = dateobss
tab['DATE'] = dates
tab['AIRMASS'] = airmass
tab['EXPTIME'] = exptime
tab['PROPID'] = propids



### PROGRAM ##########
ys = np.loadtxt(young_stars_input_catalog_filename, comments='#', dtype=str)
gaia_id = [int(x) for x in ys[:,0]]
tmass_id = ys[:,1]
mask = np.logical_or(np.in1d(tab['OBJNAME'], gaia_id), np.in1d(tab['OBJNAME'], tmass_id))




print tab
tab.write('files.fits', format='fits', overwrite=True)


print t
print len(t)

objnames = [x[1] for x in t]
print set(objnames)
print len(set(objnames))
