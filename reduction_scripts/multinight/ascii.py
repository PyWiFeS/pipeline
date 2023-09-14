"""
Extract ascii spectra.
"""
from __future__ import print_function, division
import numpy as np
from astropy.io import fits
import os

import process_stellar as ps

root = '/data/mash/marusa/2m3reduced/wifes/'

steps = ["08", "09", "10"]

"""Generate ascii spectra for each of the data reduction steps in steps.
Saves ascii spectra in a new folder in the night directory called 'ascii'.

Parameters:
-----------
root: string
    The base path to the reduced data (i.e. where the nightly folders are
    stored.)
night: string
    The night to extract spectra from (in the form 201XYYZZ).
steps: list of strings
    The PyWiFeS data reduction steps to extract and convert to 1D spectra.
"""


#~ print('Converting to ascii:', data_dir)


# Extract all files
for path, subdirs, files in os.walk(root):
    if 'reduced_b' in path or 'reduced_r' in path:
        pass
    else:
        continue

    # Sort out directory structures
    #~ print('path', path)
    #~ data_dir = os.path.join(root, night)
    out_dir = os.path.join(path, 'ascii')

    if not os.path.isdir(out_dir) and not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    for name in files:
        fl=os.path.join(path, name)
        
        if not fl.endswith('fits'):
            continue
        
        step = fl.split(".")[-2][1:]

        # Only run on specified data reduction outputs
        if step in steps and fl.endswith('%s.fits' % step):
            print(fl)
            f = fits.open(fl)
            header = f[0].header
            objectid = header['OBJNAME']
            run = header['RUN']
            obsdate = header['DATE-OBS'] # 2019-07-12T11:47:35.5
            night = obsdate.split('T')[0].replace('-', '')
            f.close()

            # Extract spectrum
            flux, wave = ps.read_and_find_star_p08(fl)
            spectrum, sig = ps.weighted_extract_spectrum(flux)

            # Determine output format depending on spectral arm
            if 'T2m3wr' in name:
                filename = '%s_%s_%s_r.dat'%(night, objectid, step)
            elif 'T2m3wb' in name:
                filename = '%s_%s_%s_b.dat'%(night, objectid, step)
            
            filename = filename.replace(' ', '_')
            
            # Save output
            fln = os.path.join(out_dir, filename)
            print(fln)
            np.savetxt(fln, np.transpose([wave, spectrum]))

