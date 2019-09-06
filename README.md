# pipeline
The Python data reduction pipeline for WiFeS

# Modified pipeline to reduce young stars
## Prepare metadata
`python generate_metadata_script_marusa.py configRed.py /priv/mulga1/marusa/2m3data/wifes/20190314`
A new folder is created:
`/priv/mulga1/marusa/2m3reduced/wifes/20190114/reduced_r_ys` and `ys_metadata_WiFeSRed.py` is stored there. Review this file to check if everything is OK.
## Rename metadata file
## Reduce each night
`marusa@mash:~/reduce_young_stars/logs>python reduce_red_data_python3.py configRed.py /priv/mulga1/marusa/2m3data/20190314`
How does this script know what metadata file to read?
## Convert to ascii
`convert_fits_to_ascii.py`
## Wavelength calibration
TODO

# TODO
- Blue
- Simplify:
  - Copy all necessary files into `pipeline`
  - Make it more simple and robust.
- What happens if there are 2 or more stars in the slit?
- Add `convert_to_ascii` file and add obsdate_objectid.
- Create a script for: fits file with all essential info: RA, DEC, exact time (MJD?), airmass, ...
