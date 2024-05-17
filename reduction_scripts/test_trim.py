import glob
import os
from pywifes.pywifes import trim_fits_file



lista  = glob.glob("*.fits")
# lista.remove('OBK-621856-WiFeS-Red--UT20240212T181235-0.fits')

for inimg_path in lista:
    trim_fits_file(inimg_path, outimg_prefix="trim_")


os.system('mv trim_*.fits ../.')