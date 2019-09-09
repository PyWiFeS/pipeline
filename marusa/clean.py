"""
Delete auxiliary files. Keep only .p08.fits
"""
import os
import sys

folder = sys.argv[1]

filenames = os.listdir(folder)
filenames = [os.path.join(folder, x) for x in filenames if not x.endswith('.p08.fits')]

for x in filenames:
    os.remove(x)
