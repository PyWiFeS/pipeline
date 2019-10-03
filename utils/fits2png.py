"""
Diagnostics: Convert all fits files in the folder into png/jpg and save them into separate folders based on their imagetype. This is to have a quick look if all imagetypes are correct.
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import sys
import os

# Folder with fits files
root = sys.argv[1]

for path, subdirs, files in os.walk(data_dir):
    for name in files:
        fl=os.path.join(path, name)
        if fl.endswith('.fits') or fl.endswith('.fit') or fl.endswith('.FITS') or fl.endswith('.FIT'):
            print fl
            
            # Read data
            f=fits.open(fl)
            header = f[0].header
            objectid=header['OBJNAME']
            run=header['RUN']
            imagetype=header['IMAGETYP']
            f.close()
            image_data = fits.getdata(image_file, ext=0)


            # Make an output folder
            out_dir = os.path.join(root, 'png', imagetype.lower())
            print 'out_dir', out_dir
            #~ if not os.path.isdir(out_dir) and not os.path.exists(out_dir):
                #~ os.mkdir(out_dir)

            # Plot and save figure
            #~ fig = plt.figure()
            #~ fig.imshow(image_data, cmap='gray')
            #~ plt.colorbar()
            fl_out = os.path.join(out_dir, name.split('.')[:-1]+'.png')
            print 'fl_out', fl_out
            #~ plt.savefig(fl_out)
