"""
Some fits files are from the Imager (?) from daytime testing??

Delete them in all subfolders.
"""

import os
import sys


path = sys.argv[1]

#~ files = []
# r=root, d=directories, f = files
for r, d, f in os.walk(path):
    for file in f:
        if 'T2m3ag' in file:
            filename = os.path.join(r, file)
            #~ files.append()
            print filename
            
            
            #~ os.remove("ChangedFile.csv")
