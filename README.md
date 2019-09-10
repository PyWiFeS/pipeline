# pipeline
The Python data reduction pipeline for WiFeS

# Modified pipeline to reduce young stars
Scripts are in `pipeline/reduction_scripts/`

# What's changed in this version:
- This version is adjusted for reductions of many nights. In principle you only need to provide a path to the folder for each night, the rest of params in the files stay constant.
- generate metadata: a separate script that finds adjacent arcs for each science exposure. Output is not a pickle but python file that is imported into reduction script.
- config file: everything you need to set. You don't have to set anything in generate metadata or reduction script anymore.
## Prepare metadata
`python generate_metadata_script_marusa.py configRed.py /priv/mulga1/marusa/2m3data/wifes/20190314`
A new folder is created:
`/priv/mulga1/marusa/2m3reduced/wifes/20190114/reduced_r_ys` and `ys_metadata_WiFeSRed.py` is stored there. Review this file to check if everything is OK.
## Rename metadata file
## Reduce each night
`python reduce_red_data_python3.py configRed.py /priv/mulga1/marusa/2m3data/wifes/20190314`
How does this script know what metadata file to read?
## Convert to ascii
`convert_fits_to_ascii.py`
## Radial velocity determination and correction
TODO
## List of all files with objnames, airmass, instrument etc.
`pipeline/marusa/files.fits`, generated with `pipeline/marusa/list_of_all_stars_observed_from_fitsfiles.py` - it would be great if we could add 'Program' column to mark objects, tellurics, specphot etc. but would require more work.

# TODO
- Blue
- Simplify:
  - Copy all necessary files into `pipeline`
  - Make it more simple and robust.
- What happens if there are 2 or more stars in the slit?
