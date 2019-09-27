"""
Print fits header

"""

from astropy.io import fits
import sys

filename = sys.argv[1]

try:
    keyword = sys.argv[2]
except:
    keyword=None

if keyword is None:
    print fits.getheader(filename, 0)
else:
    print fits.getheader(filename, 0)[keyword]
