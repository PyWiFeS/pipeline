"""
Some fits files are from the Imager (?) from daytime testing??

Delete them in all subfolders.
"""

import os
import sys


path = sys.argv[1]

for r, d, f in os.walk(path):
    for file in f:
        if 'T2m3ag' in file:
            filename = os.path.join(r, file)
            print filename
            os.remove(filename)
