"""
Delete auxiliary files. Keep only .p08.fits
"""
import os
import sys

folder = sys.argv[1]

filenames = os.listdir(folder)

# What to delete: everything except p08-p11 and metadata
filenames = [os.path.join(folder, x) for x in filenames if not x.endswith('.p08.fits') and not x.endswith('.p09.fits') and not x.endswith('.p10.fits') and not x.endswith('.p11.fits') and 'metadata' not in x]


for x in filenames:
    print('Delete %s'%x)
    os.remove(x)
print('DONE.')
