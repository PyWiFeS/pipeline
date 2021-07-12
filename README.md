# pipeline
The Python data reduction pipeline for WiFeS

### July 2021: How is this code different from the main pyWiFeS repository?
- Update to python3 and bug fixes (Adam fixed lots of them)
- We should keep the old structure, and only add a set of scripts to tackle multinight observations:
-- `reduction_scripts` still contains `reduce_blue_data.py`, `reduce_red_data.py`, `save_blue_metadata.py` and `save_red_metadata.py` just as before. We added `multinight` folder with our multinight scripts.



- There is one config.py file now for both blue and red. Filter (band: blue/red) is set in the command line.
- How to run: `python reduce_stellar.py config.py obsdate band`, e.g. `python reduce_stellar.py config.py 20181128 r` for red band

- reduction steps for the original pywifes are on http://www.mso.anu.edu.au/pywifes/doku.php?id=documentation

### Known problems
- Is there still a problem with pickle files?

### For Marusa's scripts to work, don't forget to add these lines (or modifies) in the `.bashrc` file:
`export PYTHONPATH="/priv/mulga1/marusa/reduction_wifes/pipeline/src/:$PYTHONPATH"
export PYWIFES_DIR=/priv/mulga1/marusa/reduction_wifes/pipeline/reference_data/
export PYTHONPATH="/priv/mulga1/marusa/reduction_wifes/pipeline/tools-master/:$PYTHONPATH"
export PYTHONPATH="/priv/mulga1/marusa/reduction_wifes/pipeline/reduction_scripts/:$PYTHONPATH"
export PYTHONPATH="/priv/mulga1/marusa/reduction_wifes/pipeline/reduction_scripts/extra_tools/:$PYTHONPATH"`

### How to run PyWiFeS reduction (stellar)
- Update list of calibration files if using: `utils/find_all_calib_files.py`.
- Edit `reduction_scripts/configBlue.py` and `reduction_scripts/configRed.py`. It is best to rename these files so they are unique and don't get overwritten when you do `git pull`.
- Generate a metadata file in `reduction_scripts`: `python generate_metadata_script_stellar.py configRed.py 20190321`
If you need to specify the nights with calibration files, use `python generate_metadata_script_stellar.py configRed.py 20190321 -b 2019018 -f 20190323`
- Check metadata file.
- Run `python reduce_marusa.py configBlue.py 20190321` in `reduction_scripts`.
- Convert to ascii: `utils/pywifes_utils.extract_stellar_spectra_ascii('/data/mash/marusa/2m3reduced/wifes/', '20190321')`
- Make a pdf with spectra: `plumage/plotting.plot_nightly_spectra('/data/mash/marusa/2m3reduced/wifes/', '20190321', plot_output_path='/data/mash/marusa/2m3reduced/wifes/')`
- Delete auxiliary reduction files (because there is a lot of data): `python utils/clean.py /data/mash/marusa/2m3reduced/wifes/20190321/reduced_b/` and same for `reduced_r`.

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
- `files.fits`: Marusa: Add a column `Program` (young stars), Adam (tellurics, spphot, standard, TESS), the rest is 'Other/...'
- `files.fits`: Marusa: A column `Useful`: crossmatch against a separate text file with bad data
- Adam: generate metadata: Automatically crossmatch against the list of tellurics. And add Bohlin 2013 HST flux standards to pywifes list of standards.
- Blue
- Simplify:
  - Copy all necessary files into `pipeline`
  - Make it more simple and robust.
- What happens if there are 2 or more stars in the slit?
- Diagnostics: It does happen that imagetype is wrong. E.g. Arcs are sometimes not arcs and then image reduction keeps crashing. Make a script to make a simple convert of all fits files in the folder into jpgs sorted into different folders based on their type. This would allow a quick look to check if all images are what they say they are.
# TODO2
- find_all_calib_files: things are hardcoded there
- make sure all paths are in an untracked file, so when you git pull it doesn't overwrite your settings
- put back the option to specify object name for the objects you want to reduce
- enable individual reduction, meaning that if there are many frames of the same object that night and you don't want to co-add them but reduce separately. Manage arcs in metadata: they should go together with only one science frame.
- Don't hardcode things. Put everything in config files!!!
- object_list: doesn't work right now!

- generate_metadata: It only assigns one arc to RZ Mic and many arcs to other science targets: /priv/mulga2/nstanden/20181018/...
- generate metadata: import pu: remove this
- include an option to chose a way of how arcs are matched with science: add coordinates
- Propose missing calibration files: it's BROKEN!!
