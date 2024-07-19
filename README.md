# PyWiFeS 

The automated Python data reduction pipeline for WiFeS. 

## [July 2024 updates for this fork of the automation branch]

Forked from PyWiFeS/pipeline [commit 45c69d8] in July 2024


What's been done:

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

Tested with 1x2 Full Frame (B/R3000), 1x1 Full Frame (B/R3000), 1x2 Full Frame Nod & Shuffle (B/R3000), 1x2 Stellar (B/R3000), 1x2 TAROS Full Frame (B/R3000)


***Future Development Underway***

- improve fitting of flatfield and standard star shapes

- additional backwards compatibility with TAROS data

- add option to use [Astro-SCRAPPY](https://astroscrappy.readthedocs.io/en/latest/), a fast Cython/C implementation of LACosmic (requires installation via pip)

- allow for position shifts (and strange PSF) between coadded frames when extracting the STD (it also means the N*nanmean estimate of the sum will go wrong where there are NaNs)

- modify wavelength treatment in cube generation to better reflect native resolution of the VPH gratings

- test additional gratings, beam-splitters, observing modes


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
- New `.JSON` config files provided for each observing mode (`classic` and `nod-and-shuffle`) and grating. 
  - The pipeline chooses the template automatically.
  - Users don't need to set anything in generate metadata or reduction scripts anymore.
  - Users can create their own `.JSON` file following the same structure as their preferred setup.
  - Another set of `.JSON` config files is provided for when the pipeline aims to generate the master calibration files only.
- Logger file to track the data usage and pipeline performance.
- Astrometry is now implemented in the data cubes although the accuracy is 
- Extraction and splice of the spectra and splice of the 3D astrometrised cubes are now implemented.
- Multiprocessing can be enabled so the pipeline can run faster (see caveats below). 
- Multiple quality plots are automatically generated and saved.
- Outputs organised in a new `/data_products` directory as per below.


## User Manual
For more information, we refer the users to the [**PyWiFeS User Manual**](https://www.mso.anu.edu.au/pywifes/doku.php?id=documentation). This manual explains the general structure of the pipeline, the steps of the data reduction, or technical details about the Python modules and functions.

## Installation
1. Download or clone this fork of the `pipeline` repository in the `automation` branch:
    ```sh
   git clone -b automation https://github.com/conken/pipeline.git
   ```
2. Set up a python environment (via conda, for example) with:
    -python 3.10
    -scipy 1.9.1
    -pip

3. From the pipeline root directory, run:
    ```sh
   pip install .
   ```
4. Point the environment variable `PYWIFES_DIR` to your reference data directory. There are a few possible ways to do that:
    1. In your conda env, run the following:
    '''conda env config vars set PYWIFES_DIR=/your/path/to/pipeline/reference_data'''
    '''conda env config vars unset PYTHONPATH'''
    Then deactivate and reactivate the conda env.
    2. Add the following line to your `~/.bashrc` so it will run on login:
    ```sh
    export PYWIFES_DIR=/.../.../pipeline/reference_data
    ```
    3. Or run the command manually before 'Running the Pipeline'.
    4. Alternatively, if `PYWIFES_DIR` is not set, the pipeline searches the program's *install* directory.
    For this approach to work, you would instead need to install with `pip install -e .`
5.  Set up an alias for the main reduction routine `reduce_data.py`: 
    ```sh
    alias pywifes-reduce='/.../.../pipeline/reduction_scripts/reduce_data.py'
    ```    

## Running the Pipeline
1. Put all raw data and calibration files in a dedicated directory, e.g. `/.../.../working_directory/my_raw_data`
2. Run the main reduction routine, giving the raw data directory path as an input parameter. The pipeline will run both arms automatically and choose the observing mode by checking the headers.
    ```sh
   pywifes-reduce my_raw_data
   ```

### User-Defined Reduction Parameters

**Set reduction steps**

To specify the reduction steps for blue and red data other than the defaults, users can provide the paths to the respective JSON files using the `--red-params` and `--blue-params` flags as follows:
   

    pywifes-reduce my_raw_data --red-params /.../.../user_red_param_file.json --blue-params /.../.../user_blue_param_file.json


**Reduce Data Using Master Calibration Files**

To perform data reduction using master calibration files from previous reductions, use the `--from-master` flag along with the path to the directory containing all calibration files. Both blue and red data should be stored together in the same directory for coherence. If no path is specified after the `--from-master` flag, the default directory `./data_products/master_calib` is assumed.

    pywifes-reduce my_raw_data --from-master /.../.../master_calibrations_directory


**Other parameters**

`-skip-done`: Skip already completed steps. This will check for the existence of the output files from each step, as well as whether the output files are newer than the input files -- if either part fails, the step will be executed.

`-just-calib`: Triggers the data reduction in the absence of on-sky data (both science and calibration). It only produces basic calibration files.

'-extract-and-splice': Automatically locate sources in the output datacubes, extract sources with parameters defined in JSON files:

    /../../pipeline/pipeline_params/params_extract_{obs_mode}.json

The sky annulus is hard-coded to extend from 3 to 4 times the JSON-specified source radius. The pipeline uses 2nd-order Lanczos (sinc) interpolation to map the red arm onto the finer wavelength spacing of the blue arm (the red arm wavelength spacing is 60% coarser in the default JSON setup).

### Extra usabilities
#### Multiprocessing
When multiprocessing is enabled, the pipeline *may* do the job faster. This will depend on the operative system used to run the pipeline. The multiprocessing setup is recommended for **Linux** users, as they will see a significant improvement in the computation time. On the other side, Mac OS users might get a similar running time (or just slightly faster) than in one-process mode. 
To enable the multithreading option, please follow these steps:

1. Open the `.json` file that corresponds to your data observing mode and grating. That is, `reduction_scripts/pipeline_parms/params_class_<grating>.json` for classic mode, or `reduction_scripts/params_ns_<grating>.json` for nod-and-shuffle.
2. Set `"multithread": true` in all the cases. There should be a total of 6 `"multithread"`, 3 for each of the blue and red arms in the following steps: `"step": "wave_soln"`, `"step": "cosmic_rays"`, and `"step": "cube_gen"`.
3. [Optional] Set `max_processes` to the *maximum* number of sub-processes you would like to launch for `"step": "cosmic_rays"`, and `"step": "cube_gen"`. If `-1`, the pipeline will use as many processes as there are hardware & logical cores on your device, which may be larger than the number of *available* cores, e.g. for Slurm users. Limiting the number of sub-processes can improve the efficiency and availability of your device.
4. Run the pipeline following the instructions above.

#### Skip steps 

Some steps in the data reduction process can be skipped by setting `"run": false` in the corresponding step in the `.json` files. However, in some cases, the step cannot be skipped as it is required for the pipeline to continue reducing the data. For example, the wavelength solution is always required for a successful data reduction. Other steps such as the flux calibration, the extraction of the standard star, or the telluric correction can be skipped in case of, for example, missing calibration files. 


### DATA REDUCED
The pipeline will generate the `data_products` directory within the working directory 
`/.../.../working_directory` containing the reduced data, a logger file to track the information from the reduction process, and the following structure (the Splice files only appearing if the `-extract-and-splice` flag is used): 

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
If you encounter any issues or have suggestions for improving the pipeline, please [**open a new issue**](https://github.com/conken/pipeline/issues) in the `issues` tab and fill out the provided template. Your feedback is very valuable!






