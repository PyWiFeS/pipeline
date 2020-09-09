## Data reduction scripts for easier reduction of multiple nights of stellar data
Data reduction steps:
- (1) Edit `config.py` file
- (2) [optional] Prepare a file with all available calibration files
- (3) Prepare metadata files
- (4) Run reduction scripts


### Data Calibration files (bias, flat, dark)
If you would like to use calibration files from other nights, 
- Update list of calibration files if using: `utils/find_all_calib_files.py`.

### Prepare `config.py` file
- Edit `reduction_scripts/configBlue.py` and `reduction_scripts/configRed.py`. It is best to rename these files so they are unique and don't get overwritten when you do `git pull`.

- `input_root`: folder where nightly folders are stored with raw data.
- `output_root`: where to save reductions with nightly data.
- `prefix`: [default: None][optional] Save to folders with this prefix.
- `calmin`: [default: 3]. Number of each of the calibration frames.
- `coadd_images`: [default: True]. Co-add images of the same object. # If false, then separate them in the metadata file. Take care with arcs.
- `multithread`: [default: False]. Does this work?
- `skip_done`: [default: True]. Skip reduction steps that were already done previously.
- `excluderun_filename`: [default: None] = '/data/mash/marusa/2m3data/wifes/list_of_bad_exposures_that_we_shouldn_use.dat'. # Run numbers of bad files that you don't want to include in the reduction
- `badcalib_filename`: [default: None] = '/data/mash/marusa/2m3data/wifes/list_of_high_biases_pay_attention.dat' # None. # List of bad calibration files
- `object_list_filename`: [default: None] = '/data/mash/marusa/2m3data/wifes/reduce_these_objects.dat'. Object list: reduce only these objects
- `object_list`: [default: None] = ['RZ Mic']. reduce only these objects

- [Is this true?] config file: everything you need to set. You don't have to set anything in generate metadata or reduction script anymore.

### Prepare metadata
- Generate a metadata file in `reduction_scripts`: `python generate_metadata_script_marusa.py configRed.py 20190321`
If you need to specify the nights with calibration files, use `python generate_metadata_script_marusa.py configRed.py 20190321 -b 2019018 -f 20190323`

- A new folder is created:
`/priv/mulga1/marusa/2m3reduced/wifes/20190114/reduced_r_ys` and `ys_metadata_WiFeSRed.py` is stored there. Review this file to check if everything is OK.

- [Is this true?] generate metadata: a separate script that finds adjacent arcs for each science exposure. Output is not a pickle but python file that is imported into reduction script.

- Check metadata file.
- Rename metadata file. (?)

# Flux and telluric calibration files
Check if they are on the list. If not, add them, and provide Hubble files.

### Reduce data
- Run `python reduce_marusa.py configBlue.py 20190321` in `reduction_scripts`.

## Optional 
### Radial velocity determination and correction
TODO

### Convert to ascii
`utils/pywifes_utils.extract_stellar_spectra_ascii('/data/mash/marusa/2m3reduced/wifes/', '20190321')`

??? `convert_fits_to_ascii.py`

### Plot a pdf with all the spectra of the night
`plumage/plotting.plot_nightly_spectra('/data/mash/marusa/2m3reduced/wifes/', '20190321', plot_output_path='/data/mash/marusa/2m3reduced/wifes/')`


### Delete auxiliary reduction files
There is a lot of data. `python utils/clean.py /data/mash/marusa/2m3reduced/wifes/20190321/reduced_b/` and same for `reduced_r`.

### List of all files with objnames, airmass, instrument etc.
`pipeline/marusa/files.fits`, generated with `pipeline/marusa/list_of_all_stars_observed_from_fitsfiles.py` - it would be great if we could add 'Program' column to mark objects, tellurics, specphot etc. but would require more work.

### Notes
`config.py`, `generate_metadata_script_marusa.py` and `reduce_red_data_python3.py` are customized files used to reduce young stars on `mash`.

### Advanced
`step: cube_gen` (`suffix: '08'`): `dw_set`, `wmin_set` and `wmax_set` are now set in the `reduce_stellar.py`, based on the band (Blue/Red).



# TODO [review this]
- `files.fits`: Marusa: Add a column `Program` (young stars), Adam (tellurics, spphot, standard, TESS), the rest is 'Other/...'
- `files.fits`: Marusa: A column `Useful`: crossmatch against a separate text file with bad data
- Adam: generate metadata: Automatically crossmatch against the list of tellurics. And add Bohlin 2013 HST flux standards to pywifes list of standards.
- Blue
- Simplify:
  - Copy all necessary files into `pipeline`
  - Make it more simple and robust.
- What happens if there are 2 or more stars in the slit?
- Diagnostics: It does happen that imagetype is wrong. E.g. Arcs are sometimes not arcs and then image reduction keeps crashing. Make a script to make a simple convert of all fits files in the folder into jpgs sorted into different folders based on their type. This would allow a quick look to check if all images are what they say they are.
# TODO2 [review this]
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
