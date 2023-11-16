"""
Extract ascii spectra.
"""
from __future__ import print_function, division
import glob
import numpy as np
from astropy.io import fits
import os

import process_stellar as ps

#~ root = '/data/mash/marusa/2m3reduced/wifes/'
root = '/data/mash/marusa/2m3reduced/wifes/20190722/'

#~ steps = ["08", "09", "10"]
steps = ["08", "09"]

"""Make a fits data cube of the reduced and extracted stellar spectra, with
a different HDU for each data reduction step. Store this new cube in a 
folder called 'extracted_1D_cubes' in the night folder. By default:
 - 8: Final (non-fluxed or telluric corrected)
 - 9: Final (fluxed)
 -10: Final (fluxed and telluric corrected)

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

out_dir = os.path.join(root, 'extracted_fits_1d')

if not os.path.isdir(out_dir) and not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Get a list of the of the reduced files. Assume files have the naming
# convention T2m3w[r/b]-20190828.104313-0031.p08.fits' where [r/b] will be
# a single character indicating the arm of WiFeS used.
unique_obs = glob.glob(os.path.join(root, "*", 'reduced_*', "*%s.fits" % steps[0]))

# Generalise by removing the suffixes
unique_obs = [ob.replace(".p%s.fits" % steps[0], "") for ob in unique_obs]

# Now go through one observation at a time, and create a data cube 
# containing the original header information, and a HDU for each step of
# the data reduction listed in steps.
for ob in unique_obs:
    # Get the list of fits files to extract 1D spectra from
    fits_files = [ob + ".p%s.fits" % step for step in steps]

    # Get the header information from the first
    header = fits.getheader(fits_files[0])

    obj_id = header['OBJNAME'].replace(" ","")
    night = header['DATE-OBS'].split('T')[0].replace('-', '')

    # Construct a new fits file
    hdus = []
    hdus.append(fits.PrimaryHDU(header=header))

    print("Making cube for %s%s" % (ob, ".pXX.fits"))

    for fits_file, step in zip(fits_files, steps):
        if not os.path.exists(fits_file):
            print('%s does not exist!'%fits_file)
            continue
        # Extract the 1D spectrum
        print("\tExtracting 1D spectrum for step %s" % step)
        flux, wave = ps.read_and_find_star_p08(fits_file)
        spectrum, sig = ps.weighted_extract_spectrum(flux)

        # Make fits table from numpy record array
        data = np.array([wave, spectrum, sig]).T.tolist()
        rec = np.rec.array(data, names=["wave", "spectrum", "sigma"])

        hdus.append(fits.BinTableHDU.from_columns(rec))
        hdus[-1].header["REDSTEP"] = step
    
    # Determine output format depending on spectral arm
    if 'T2m3wr' in ob:
        output_filename = '%s_%s_r.fits' % (night, obj_id)
    elif 'T2m3wb' in ob:
        output_filename = '%s_%s_b.fits' % (night, obj_id)

    output_filepath = os.path.join(out_dir, output_filename)

    # Write the fits file
    print("Writing cube to %s \n" % output_filename)
    hdu_list = fits.HDUList(hdus)
    hdu_list.writeto(output_filepath, overwrite=True)



