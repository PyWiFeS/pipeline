"""
Convert all .p08.fits files to ascii spectra.
"""

import numpy as np
from astropy.io import fits
import sys
import os

# PyWiFeS
import process_stellar as ps # WHEN running this without ssh -Y, I get RuntimeError: Invalid DISPLAY variable

#~ root = '/priv/mulga1/marusa/2m3reduced/wifes/'
root = '/data/mash/marusa/2m3reduced/wifes/'
#root = "/priv/mulga2/arains/ys/wifes/reduced/"

try:
    obsdate = sys.argv[1]
    data_dir = os.path.join(root, obsdate)
except:
    data_dir = root
    obsdate=None

steps = ["08", "09", "10"]

print 'Converting to ascii:', data_dir

out_dir = os.path.join(data_dir, 'ascii')
if not os.path.isdir(out_dir) and not os.path.exists(out_dir):
    os.mkdir(out_dir)

for path, subdirs, files in os.walk(data_dir):
    for name in files:
        fl=os.path.join(path, name)
        step = fl.split(".")[-2][1:]
        # Only run on specified data reduction outputs
        if step in steps and fl.endswith('%s.fits' % step):
            print fl
            f=fits.open(fl)
            header = f[0].header
            objectid=header['OBJNAME']
            run=header['RUN']
            f.close()

            flux, wave = ps.read_and_find_star_p08(fl)
            spectrum, sig = ps.weighted_extract_spectrum(flux)

            obsdate=fl.split('.')[0].split('-')[-1]

            if 'T2m3wr' in name:
                filename = '%s_%s_%s_r.dat'%(obsdate, objectid, step)
            elif 'T2m3wb' in name:
                filename = '%s_%s_%s_b.dat'%(obsdate, objectid, step)
            
            filename = filename.replace(' ', '_')
            
            fln=os.path.join(out_dir, filename)
            print fln
            print
            np.savetxt(fln, np.transpose([wave, spectrum]))
