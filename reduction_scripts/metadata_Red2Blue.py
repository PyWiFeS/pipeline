"""
Create Blue metadata from Red.
The thing is that blue and red fits files might have slightly different timestamps, so it is not enough just to change b and r in filenames.

Usage: obsdate and prefix!!!

"""

import sys
import os
import imp




root = '/priv/mulga1/marusa/2m3reduced/wifes/'

obsdate = sys.argv[1]

def get_metadata(obsdate, band='r'):
    try:
        prefix = sys.argv[2]
        folder = 'reduced_%s_%s'%(band, prefix)
    except:
        folder = 'reduced_%s'%band
        

    folder = os.path.join(root, obsdate, folder)
    print folder

    files = []
    for r, d, f in os.walk(folder):
        for file in f:
            if 'metadata' in file and '.pyc' not in file and 'mode' not in file:
                files.append(os.path.join(r, file))

    #~ print files
    # Read red metadata dict
    for filename in files:
        print filename
        config = imp.load_source(filename.replace('.py', ''), filename)
    print
    
    return config

metadata_red = get_metadata(obsdate, band='r')
metadata_blue = get_metadata(obsdate, band='b')

red = metadata_red.night_data
blue = metadata_blue.night_data

for k, r in red.iteritems():
    b=blue[k]
    
    print k
    print b
    print r
    print
    print


print 
print metadata_blue.night_data
