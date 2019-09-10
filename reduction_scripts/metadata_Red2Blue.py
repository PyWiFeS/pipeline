"""
Create Blue metadata from Red.
The thing is that blue and red fits files might have slightly different timestamps, so it is not enough just to change b and r in filenames.
"""

import sys
import os

root = '/priv/mulga1/marusa/2m3reduced/wifes/'

obsdate = sys.argv[1]





try:
    prefix = sys.argv[2]
    print prefix
    #~ filename_red = os .path.join(filename_red, prefix)
    
    red_folder = 'reduced_r_%s'%prefix
    print red_folder
except:
    red_folder = 'reduced_r'
    

red_folder = os.path.join(root, obsdate, red_folder)
print red_folder

files = []
for r, d, f in os.walk(red_folder):
    for file in f:
        if 'metadata' in file and '.pyc' not in file:
            files.append(os.path.join(r, file))

print files
# Read red metadata dict
