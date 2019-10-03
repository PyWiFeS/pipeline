"""
Diagnostics: Convert all fits files in the folder into png/jpg and save them into separate folders based on their imagetype. This is to have a quick look if all imagetypes are correct.
"""


import matplotlib.pyplot as plt
import sys

# Folder with fits files
root = sys.argv[1]


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




plt.figure()
plt.imshow(image_data, cmap='gray')
plt.colorbar()
