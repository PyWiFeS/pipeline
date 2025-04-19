# PyWiFeS 

The automated Python data reduction pipeline for WiFeS. 

## Updates

What's been done [20250419]:

- allow user to propagate applied extinction correction into datacubes and extracted spectra

- allow user to separate nod-and-shuffle observations and process independently rather than subtracting

- improved edge-handling for atmospheric differential refraction remapping

- use enhanced flatfield algorithm when lacking twilight flats

- improve spectral extraction for binning modes other than 1x2


### Previously ###

What's been done [20250403]:

- update standard stars to use Nov 2024 edition of CALSPEC, where available

- refine standard star flux calibration procedure

- improve handling of coadded images, especially when exposure times differ between arms

- bugfixes in spliced spectrum wavelength definition, airmass-dependence of first telluric standard

[20250324]:

- allow user to propagate telluric spectrum into datacubes and extracted spectra

- add option to use 'native' wavelength scale of spectrograph and report wavelengths in new FITS extension

- add option to only extract spectra from datacubes without also splicing the cubes and spectra

- allow user specification of sky annulus size and output wavelength step in extraction

- allow user to splice cubes and spectra that do not overlap in wavelength

- restore WCS info for detection plots when PA near 90/270


From branch fji_deps_updates-and-spec_extract [20241205]:

- Improved source detection routine by implementing DAOStarFinder (`photutils`), which identifies sources based on local density maxima and fits a 2D Gaussian model.

- Updated the extraction parameters JSON file (`./pipeline_params/params_extract.json5`), making it more structured and including the parameters required for the new routine.

- Introduced the `--extract-params` command-line argument, allowing users to pass a tailored extraction configuration file. This is particularly useful for local users who need to adjust detection and extraction parameters (e.g., for detecting faint sources, changing the default aperture size for extraction).

- Update dependencies to photutils >= 2.0.0 (requiring scipy requirement update) and to allow any numpy version.

[20250115]:

- new NeAr and CuAr line list from NIST (Ritz wavelengths in air)

- option to output datacubes in vacuum wavelengths instead of air

- propagate various calibration keywords into science data

- clean up logging

- update method of standard star association with observations

[20241114]:

- all command-line arguments now use double-dash

- new `--reduce-both` option to reduce arms simultaneously

[20241023]:

- check for mixtures of gratings and beamsplitters (which can alter slitlet locations on CCD) and warn users

- smoother handling of processing from existing master_calib folder when also using night-specific flux/telluric standards in the processing

- cut darks to half-frame, if present

[20241014]:

- update flatfield algorithm when lacking twilight flats

[20241011]:

- improve handling of flux vs. telluric standards

- when running from existing master_calib files, do not contaminate the input folder with newly generated sensitivity or telluric corrections, but move to a master_calib folder in the active data_products directory

[20241002]:

- save smooth fit to flatfield (reflecting CCD, grating, beamsplitter throughput) and temporarily remove when fitting flux standard to reduce amplitude of fluctuations

- add one pass of filtering outlier pixels from Savitzky-Golay fitting to flatfield

- have sensitivity curve better avoid absorption at the reddest wavelengths of R3000 grating

- save sky spectrum from telluric standard to allow wavelength refinement when correcting science spectra

- allow sub-aperture Nod & Shuffle extraction to find positive and negative peaks, and extract both (inverting the flux for the negative peaks)

- more reliably associate red and blue spectra amongst multiple apertures by forcing lists to be sorted and allowing for differences in DATE-OBS timestamp between arms

[20240919]:

- streamlined fitting of flatfield data, now using double Savitzky-Golay filter

- does not search for raw files in subfolders (making it easier to hide away data that shouldn't be processed)

- add interactive diagnostic plot to fitting of flux standard star

[20240830]:

- merge in documentation and other updates from official repository.

- modify methodology for source extraction.

- add user option to include DQ extension when splicing cubes and extracting sources.

- add extension with subtracted sky to extracted spectral outputs.

- constrain datatypes for spliced and extracted cubes and spectra.

- handle saturated pixels in raw frames (flag in VAR and DQ frames for SCIENCE, STANDARD, SKY images; interpolate over for other image types).

[20240823]:

- update JSON to JSON5 to allow (C-style) comments to aid user configuration setup.

- reinstitute R7000 Littrow ghost masking offsets, but remove default usage.

- revise overscan region for B5/R5 epoch to avoid region of decay from high counts, and revise overscan estimate to use mean of central 50% of distribution per row.

- apply extinction correction to flux standard star (previously calculated but not applied).

- associate any sky exposures with science exposures if IMAGETYP = 'SKY' (currently a manual modification of the image header by the user) and OBJECT keywords match. Use of `--coadd-mode prompt` allows manual association of SKY and OBJECT exposures. Use of `--coadd-mode none` blindly associates one SKY and one OBJECT using the order they appear in the image list. Multiple sky exposures are scaled by their exposure times and median-combined. Subtraction from science frame scales effective sky exposure time to that of science frame.

- avoid recopying existing, older files into 'intermediate' folder.

[20240809]:

- simplify JSON files by removing separate nod & shuffle variants

- add option to avoid coadding all science frames with the same OBJECT name (useful for time-series, position shifts, mosaics, etc.)

- make greedy definition of STANDARD stars rely on command-line argument to enable (rough) sky proximity to known standard, otherwise (new default) require IMAGETYP keyword to be 'STANDARD'

- adjust Littrow ghost positions for 7000 gratings

- partially revert obs_coadd summation to better handle position shifts between (unguided) standard star exposures, but retain 'nansafesum' method for when observations are well aligned

- simplify code by merging single-thread and multi-threaded cube generation modules

- correct atmospheric differential refraction calculation to be valid for y-binning options x2, and update standard weather parameters based on SkyMapper records

- streamline STANDARD extraction code

[20240802]:

- add very cold pixels (typically having counts at least ~4x lower than neighbors) to bad pixel mask

- optional masking of Littrow ghosts in dome and sky flats

- dome flats can vary in counts by up to ~1000 within a triplet, so scale them to effectively use median-stacking

- twilight flats have scattered light up to several hundred, so modify the cleanup threshold

- twilight flats can vary in color from night to night, but stacking near the middle is hopefully sufficient to constrain the illumination

- wires likely vary in the same way as dome flats, so scale them before stacking

- increase buffer and nsig_lim in flat_cleanup to avoid clipping real scattered light

- chunk the imcombine process into x-sections if equivalent of > 5 unbinned images

- add creation of superarc to improve S/N

- add keywords for number of images imcombined

- add option to set flux pixels to NaN if masked in DQ extenstion

- avoid making unnecessary MEF files (bias, flat, dark)

- enable diagnostic plots for flat combination; useful because:
  - there are distinct amounts of x-axis shifts in fringe pattern of red dome flats
  - there are distinct families of red-end behavior of blue dome flats

[20240724]:

- updated bad pixel mask and simplified code to add more in the future; wait to interpolate over bad pixels in STD and OBJECT frames until VAR and DQ extensions are created and can be flagged with bad pixels

- when coadding exposures, do reasonable things with the VAR and DQ extensions

- minimise file bloat by limiting BITPIX of output frames to float32 for SCI and VAR, and int16 for DQ

- handle compressed raw data inputs (fpacked [.fz] or gzipped [.gz], just using astropy.io.fits capabilities)

- avoid elevated overscan levels for rows with high counts, and estimate the overcan for those rows with Chebyshev polynomial determined from interslice regions

- add creation of superwire to coadd the multiple wires typically available and improve fitting of spatial distortions

- fix call to LACosmic to use correct readnoise and to retain previously flagged pixels in DQ extension

- refresh code to work with stellar readout mode and old Taros data

- remove just-calib JSON files; instead just skip the steps not related to calibrations from the normal JSON files (to avoid potential differences between just-calib and normal/from-master usage)

- make use of twilight flats the default for both arms in the JSON files

- make skip-done option more comprehensive in avoiding repetition, and check if input is newer than output to also trigger running of step

- make automatic extraction and splicing rely on the presence of a command-line argument

- pull some utilities into a new wifes_utils module

- for early Automated Observing, checks and corrects for shifts in red arm readout

- switch to NaN-safe statistics in most places

- add a variety of PyWiFeS settings to image headers

- corrected various plots (axis ratios, axis label)

- edits for code legibility and (most) PEP8 style conventions

Tested with 1x2 Full Frame (B/R3000, UB/RI7000), 1x1 Full Frame (B/R3000), 1x2 Full Frame Nod & Shuffle (B/R3000), 1x2 Stellar (B/R3000), 1x2 TAROS Full Frame (B/R3000), 1x2 TAROS Stellar (B/R3000)


### [May 2024 updates - use only AUTOMATION BRANCH]


#### Main upgrades from the previous PyWiFeS version ([master](https://github.com/PyWiFeS/pipeline/tree/master) branch)
- Updated to Python 3, with bug fixes and headers for the new telescope setup.
- Cleaner repository (defunct scripts and directories were removed).
- Pip installable.
- Adjusted to read and reduce data automatically all at once in both arms and for all the observing configurations, including:
    - `classic` mode.
    - `nod-and-shuffle` mode.
    - `sub-nod-and-shuffle` mode.
    - Stellar frame configuration (that is, only the middle part of the detector is used).
    - Binning along the spatial axis.
- New `.JSON` config files provided for each grating. 
  - The pipeline chooses the template automatically.
  - Users don't need to set anything in generate metadata or reduction scripts anymore.
  - Users can create their own `.JSON` file following the same structure as their preferred setup.
- Logger file to track the data usage and pipeline performance.
- Astrometry is now implemented in the data cubes although the accuracy is uncertain.
- Extraction and splice of the spectra and splice of the 3D astrometrised cubes are now implemented.
- Multiprocessing can be enabled so the pipeline can run faster (see caveats below). 
- Multiple quality plots are automatically generated and saved.
- Outputs organised in a new `/data_products` directory as per below.


## User Manual
For more information, we refer the users to the [**PyWiFeS User Manual**](https://www.mso.anu.edu.au/pywifes/doku.php?id=documentation). This manual explains the general structure of the pipeline, the steps of the data reduction, or technical details about the Python modules and functions.

## Installation
1. Download or clone this fork of the `pipeline` repository in the `automation` branch:
    ```sh
   git clone -b automation https://github.com/PyWiFeS/pipeline.git
   ```
2. Set up a python environment (via conda, for example) with:

    -python 3.10

    -scipy 1.9.1

    -numpy < 2.0

    -pip

3. From the pipeline root directory, run:
    ```sh
   pip install .
   ```
4. Point the environment variable `PYWIFES_DIR` to your reference data directory. There are a few possible ways to do that:
    1. In your conda env, run the following:
    ```conda env config vars set PYWIFES_DIR=/your/path/to/pipeline/reference_data```
    ```conda env config vars unset PYTHONPATH```
    Then deactivate and reactivate the conda env.
    2. Add the following line to your `~/.bashrc` so it will run on login:
    ```sh
    export PYWIFES_DIR=/.../.../pipeline/reference_data
    ```
    3. Or run the command manually before 'Running the Pipeline'.
    4. Alternatively, if `PYWIFES_DIR` is not set, the pipeline searches the program's *install* directory.
    For this approach to work, you would instead need to install with `pip install -e .`
5.  If desired, set up an alias for the main reduction routine `reduce_data.py`: 
    ```sh
    alias pywifes-reduce='/.../pipeline/reduction_scripts/reduce_data.py'
    ```    

## Running the Pipeline
1. Put all raw data and calibration files in a dedicated directory, e.g. `/.../working_directory/my_raw_data`
2. Run the main reduction routine, giving the raw data directory path as an input parameter. The pipeline will run both arms automatically and choose the observing mode by checking the headers.
    ```sh
   pywifes-reduce my_raw_data
   ```

### User-Defined Reduction Parameters

**Set reduction steps**

To specify the reduction steps for blue and red data other than the defaults, users can provide the paths to the respective JSON or JSON5 files using the `--red-params` and `--blue-params` flags as follows:
   

    pywifes-reduce my_raw_data --red-params /.../user_red_param_file.json5 --blue-params /.../user_blue_param_file.json5


**Reduce Both Arms in Parallel**

Processing may be sped up by processing both arms simultaneously. Obviously, this entails utilising more of the machine's resources. To enable, use the `--run-both` flag as follows:

    
    pywifes-reduce --run-both my_raw_data

    
**Reduce Data Using Master Calibration Files**

To perform data reduction using master calibration files from previous reductions, use the `--from-master` flag along with the path to the directory containing all calibration files. Both blue and red data should be stored together in the same directory for coherence. If no path is specified after the `--from-master` flag, the default directory `./data_products/master_calib` is assumed.

    pywifes-reduce my_raw_data --from-master /.../.../master_calibrations_directory


**Controlling How OBJECT Frames Are Coadded**

The default pipeline behaviour is to coadd exposures that share the same OBJECT header keyword, which is unhelpful when there are position shifts of the source in the IFU (due to tracking offsets, dithering, mosaicking, PA changes, etc.) or the observations form a time series. With the `--coadd-mode` command line argument, one may choose:

    `all`: (default) coadd all observations with the same OBJECT name;
    `none`: treat all IMAGETYP='OBJECT' exposures separately;
    `prompt`: prompt the user to select which exposures to coadd together.

With the `prompt` option, the user makes the choice independently for each arm. The choices are saved for the next time the pipeline is run with that dataset (as long as the `--coadd-mode prompt` option is used each time).

Example:

    pywifes-reduce my_raw_data --coadd-mode prompt

**Association of SKY Frames with OBJECT Frames**

While Nod & Shuffle observations automatically subtract the sky in the `run_sky_sub` step, other OBJECT frames need to have one or more SKY frames associated with it. At present, this relies on a manual modification of the FITS headers of such SKY images: set IMAGETYP = 'SKY'. 

Under the default `--coadd-mode all` option, associations will be made if IMAGETYP='SKY' and the OBJECT keyword matches that of the science (IMAGETYP='OBJECT') frame. 

If the `--coadd-mode prompt` option is used, the user can associate one or more SKY images with the OBJECT images of their choice. As above, this is done independently for each arm, and the choices are saved in case the pipeline is interrupted.

When `--coadd-mode none` is used, a single sky frame is associated with each science frame having the same OBJECT keyword -- the first sky with the first science, the second sky with the second science, etc.

When multiple sky frames are associated with a science image, they are scaled by their exposure time (to an equivalent EXPTIME of 1 second) and median-combined. The sky frame is scaled to the science exposure time and subtracted.

**Other parameters**

`--skip-done`: Skip already completed steps. This will check for the existence of the output files from each step, as well as whether the output files are newer than the input files -- if either part fails, the step will be executed. Note that this will not skip the extraction or splicing steps, but will overwrite existing files if `--extract` or `--extract-and-splice` are used.

`--just-calib`: Triggers the data reduction in the absence of on-sky data (both science and calibration). It only produces basic calibration files.

`--greedy-stds`: Treat observations vaguely near known standard stars as STANDARD frames even if IMAGETYP = 'OBJECT'. If this option is not set, only IMAGETYP = 'STANDARD' frames are used as standards.

`--extract`: Automatically locate sources in the output datacubes, extract sources with parameters defined in JSON5 file:

    /.../pipeline/pipeline_params/params_extract.json5

The subtracted sky spectrum is included in the output. If the inputs are Nod & Shuffle frames, no sky is subtracted. Users may choose whether to propagate the data quality (DQ), applied telluric correction (TELLURICMODEL), and corrected atmospheric extinction (EXTINCTION) to the output spectrum.

`--extract-and-splice`: In addition to the extraction described above, splice together the datacubes and extracted spectra. By default, the pipeline uses 2nd-order Lanczos (sinc) interpolation to map the red arm onto the finer wavelength spacing of the blue arm (the red arm wavelength spacing is 60% coarser in the default pipeline setup). The user may specify an alternative wavelength spacing in the extraction JSON5 file, and both arms will be interpolated to that scale.


### Extra usabilities
#### Multiprocessing
When multiprocessing is enabled, the pipeline *may* do the job faster. This will depend on the operative system used to run the pipeline. The multiprocessing setup is recommended for **Linux** users, as they will see a significant improvement in the computation time. On the other side, Mac OS users might get a similar running time (or just slightly faster) than in one-process mode. 
To enable the multithreading option, please follow these steps:

1. Open the `.json5` file that corresponds to your grating. That is, `reduction_scripts/pipeline_parms/params_<grating>.json5`.
2. Set `"multithread": true` in all the cases. There should be a total of 6 `"multithread"`, 3 for each of the blue and red arms in the following steps: `"step": "wave_soln"`, `"step": "cosmic_rays"`, and `"step": "cube_gen"`.
3. [Optional] Set `max_processes` to the *maximum* number of sub-processes you would like to launch for `"step": "cosmic_rays"`, and `"step": "cube_gen"`. If `-1`, the pipeline will use as many processes as there are hardware & logical cores on your device, which may be larger than the number of *available* cores, e.g. for Slurm users. Limiting the number of sub-processes can improve the efficiency and availability of your device.
4. Run the pipeline following the instructions above.

#### Skip steps 

Some steps in the data reduction process can be skipped by setting `"run": false` in the corresponding step in the `.json` files. However, in some cases, the step cannot be skipped as it is required for the pipeline to continue reducing the data. For example, the wavelength solution is always required for a successful data reduction. Other steps such as the flux calibration, the extraction of the standard star, or the telluric correction can be skipped in case of, for example, missing calibration files. 


### DATA REDUCED
The pipeline will generate the `data_products` directory within the working directory 
`/.../.../working_directory` containing the reduced data, a logger file to track the information from the reduction process, and the following structure (the Splice files only appearing if the `--extract-and-splice` flag is used): 

- data_products
    - `pywifes_logger.log` 
    - `xxx-Blue-UTxxx.cube.fits`
    - `xxx-Red--UTxxx.cube.fits`
    - `xxx-Splice-UTxxx.cube.fits`
    - ... 
    - `xxx-Blue-UTxxx.spec.apx.fits`
    - `xxx-Red--UTxxx.spec.apx.fits`
    - `xxx-Splice-UTxxx.spec.apx.fits`
    
    - plots
        - `UTxxx_detection_plot.png`
        - `xxx-Blue-UTxxx_.spec.ap1.png`
        - `xxx-Red-UTxxx_.spec.ap1.png`
        - ...
        - `xxx-Splice-UTxxx_.spec.ap1.png`
        - ... 
        - blue
            - `bias.png`
            - `flat_response.png`
            - ...
        - red
            - `bias.png`
            - `flat_response.png`
            - ...

    - intermediate
        - blue
            - `xxx-Blue-UTxxx.p00.fits`
            - `xxx-Blue-UTxxx.p01.fits`
            - ...
            - `xxx-Blue-UTxxx.p10.fits`
        - red
            - `xxx-Red--UTxxx.p00.fits`
            - `xxx-Red--UTxxx.p01.fits`
            - ...
            - `xxx-Red--UTxxx.p10.fits`
        - raw_data_temp
            - `xxx-Blue-UTxxx.fits`
            - ...
            - `xxx-Red--UTxxx.fits`
            - ...

    - master_calib 
        - `wifes_blue_<master_calibration>_files.fits`
        - ...
        - `wifes_red_<master_calibration>_files.fits`
        - ...

`data_products` contains the `plots` directory with the final figures of the data reduction: the 2D extracted spectra, the spliced spectra, and the sources detection plots. The figures generated during the calibration steps are saved in a different directory for each arm (`data_products/plots/arm`). 
Then, the `data_products/intermediate` directory with the calibration files generated during the data reduction process and saves the data separately for each red and blue arm. Also in `intermediate` there is the temporary directory `raw_data_temp` aimed to store the raw data and any pre-treated images (e.g. cut down calibration frames to stellar mode size, when needed) during the data reduction process. `raw_data_temp` is automatically removed when the pipeline is successfully completed. 

Finally, we find `data_products/master_calib`, which is a directory with all master calibration files produced in the data reduction. They are stored to be used in further reductions if required.


## Reporting Issues or Suggestions
If you encounter any issues or have suggestions for improving the pipeline, please [**open a new issue**](https://github.com/PyWiFeS/pipeline) in the `issues` tab and fill out the provided template. Your feedback is very valuable!
