import os

__version__ = '0.7.0'

# If you set the environment variable 'PYWIFES_DIR' then it will be found
mdir = os.getenv('PYWIFES_DIR')
# OTHERWISE YOU MUST SET IT HERE BY HAND
if mdir == None:
    metadata_dir = '/Users/fvogt/Tools/Python/pywifes/pywifes_v0.5.6/reference_data/'
else:
    metadata_dir = mdir
