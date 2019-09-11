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
    
    
    if k!='sci':
        run_r = set([x.split('-')[-1] for x in r])
        run_b = set([x.split('-')[-1] for x in b])
        
        diff = run_b.difference(run_r)
        if len(diff)>0:
            print k
            print diff
    else:
        rd=[]
        for x in r:
            print x
            for kk, rr in x.iteritems():
                for x in rr:
                    rd.append(x.split('-')[-1])
        print rd
        bl=[]
        for kk, bb in b.iteritems():
            for x in bb:
                bl.append(x.split('-')[-1])
            
        print k
        print rd
        print bl
    
    #~ r=sorted(r)
    #~ b=sorted(b)
    
    #~ print k
    #~ print b
    #~ print r
    #~ print
    #~ print
