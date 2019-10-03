"""
Edit header

Usage:
python edit_header.py T2m3wb-20190826.143357-0072.fits KEYWORD value

Imagetyp: arc, object

"""

from astropy.io import fits

import sys

filename = sys.argv[1]

keyword = sys.argv[2]
value = sys.argv[3]

fits.setval(filename, keyword, value=value)
