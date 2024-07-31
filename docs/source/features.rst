.. _features:

What's New in PyWiFeS
=====================

This pipeline automates many of the routine tasks involved in data reduction. It is written in Python and is open-source, allowing users to modify and extend the code to suit their needs.
PyWiFeS is an automated and upgraded version of the original PyWiFeS pipeline (see documentation `here <https://www.mso.anu.edu.au/pywifes/doku.php?id=documentation>`_), developed by `Childress et al. (2014) <http://adsabs.harvard.edu/abs/2014Ap%26SS.349..617C>`_


Upgrades from the Previous Version
----------------------------------

Main upgrades from the previous PyWiFeS version *as of July 2024*:

- Updated to Python 3, with bug fixes and headers for the new telescope setup.
- Cleaner repository (defunct scripts and directories were removed).
- Pip installable.
- Supports various observing configurations automatically:
    - `classic` mode.
    - `nod-and-shuffle` mode.
    - `sub-nod-and-shuffle` mode.
    - Stellar frame configuration (that is, only the middle part of the detector is used).
    - Binning along the spatial axis.
- New `.JSON` config files for each observing mode and grating.
  - The pipeline chooses the template automatically.
  - Users don't need to set anything in generate metadata or reduction scripts anymore.
  - Users can create their own `.JSON` file following the same structure as their preferred setup.
  - Another set of `.JSON` config files is provided for when the pipeline aims to generate the master calibration files only.
- Logger file to track data usage and pipeline performance.
- Implemented astrometry in the data cubes. The accuracy could be low (< 2 arcsec) in some cases. 
- Extraction and splice of the spectra and splice of the 3D astrometrised cubes are now implemented.
- Added multiprocessing for faster execution (Linux recommended).
- Multiple quality plots are automatically generated and saved.
- Organized output directory (`/data_products`).


Key Features
------------

- Fully automated data reduction processes.
- Provides tools for quality control and analysis.
- Supports batch processing of multiple objects in a single night.
- Extensible and customizable.
