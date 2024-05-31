# PyWiFeS 

The automated Python data reduction pipeline for WiFeS. 

### [May 2024 updates - use only AUTOMATION BRANCH]

***Note: This version is in ACTIVE development.***


## Main upgrades from the previous PyWiFeS version ([master](https://github.com/PyWiFeS/pipeline/tree/master) branch)
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
For more information, we refer the users to the **PyWiFeS User Manual** [LINK]. This manual explains the general structure of the pipeline, the steps of the data reduction, or technical details about the Python modules and functions.

## Installation
1. Download or clone the `pipeline` repository in the `automation` branch:
    ```sh
   git clone -b automation https://github.com/PyWiFeS/pipeline.git
   ```
2. In a Python environment and from the pipeline root directory, run:
    ```sh
   pip install .
   ```
3. Point the environment variable `PYWIFES_DIR` to your reference data directory. There are a few possible ways to do that:
    1. Add the following line to your `~/.bashrc` so it will run on login:
    ```sh
    export PYWIFES_DIR=/Users/.../pipeline/reference_data
    ```
    2. Or run the command manually before 'Running the Pipeline'.
    3. Alternatively, if `PYWIFES_DIR` is not set, the pipeline searches the program's *install* directory.
    For this approach to work, you would instead need to install with `pip install -e .`
4.  Set up an alias for the main reduction routine `reduce_data.py`: 
    ```sh
    alias pywifes-reduce='/Users/.../pipeline/reduction_scripts/reduce_data.py'
    ```    

## Running the Pipeline
1. Put all raw data and calibration files in a dedicated directory, e.g. `/Users/.../working_directory/my_raw_data`
2. Run the main reduction routine, giving the raw data directory path as an input parameter. The pipeline will run both arms automatically and choose the observing mode by checking the headers.
    ```sh
   pywifes-reduce my_raw_data
   ```

### User-Defined Reduction Parameters

**Set reduction steps**

To specify the reduction steps for blue and red data, users can provide the paths to the respective JSON files using the `--red-params` and `--blue-params` flags as follows:
   

    pywifes-reduce my_raw_data --red-params /Users/.../user_red_param_file.json --blue-params /Users/.../user_blue_param_file.json


**Reduce Data Using Master Calibration Files**

To perform data reduction using master calibration files from previous reductions, use the `--from-master` flag along with the path to the directory containing all calibration files. Both blue and red data should be stored together in the same directory for coherence. If no path is specified after the `--from-master` flag, the default directory `./data_products/master_calib` is assumed.

    pywifes-reduce my_raw_data --from-master /Users/.../master_calibrations_directory


**Other parameters**

`-skip-done`: Skip already completed steps. 

`-just-calib`: Triggers the data reduction in the absence of on-sky data (both science and calibration). It only produces basic calibration files.

### Extra usabilities
#### Multiprocessing
When multiprocessing is enabled, the pipeline *may* do the job faster. This will depend on the operative system used to run the pipeline. The multiprocessing setup is recommended for **Linux** users, as they will see a significant improvement in the computation time. On the other side, Mac OS users might get a similar running time (or just slightly faster) than in one-process mode. 
To enable the multithreading option, please follow these steps:

1. Open the `.json` file that corresponds to your data observing mode and grating. That is, `reduction_scripts/pipeline_parms/params_class_<grating>.json` for classic mode, or `reduction_scripts/params_ns_<grating>.json` for nod-and-shuffle.
2. Set `"multithread": true` in all the cases. There should be a total of 6 `"multithread"`, 3 for each of the blue and red arms in the following steps: `"step": "wave_soln"`, `"step": "cosmic_rays"`, and `"step": "cube_gen"`.
3. [Optional] Set `max_processes` to the *maximum* number of sub-processes you would like to launch for `"step": "cosmic_rays"`, and `"step": "cube_gen"`. If `-1`, the pipeline will use as many processes as there are hardware & logical cores on your device, which may be larger than the number of *available* cores, e.g. for Slurm users. Limiting the number of sub-processes can improve the efficiency and availability of your device.
4. Run the pipeline following the instructions below.

#### Skip steps 

Some steps in the data reduction process can be skipped by setting `"run": false` in the corresponding step in the `.json` files. However, in some cases, the step cannot be skipped as it is required for the pipeline to continue reducing the data. For example, the wavelength solution is always required for a successful data reduction. Other steps such as the flux calibration, the extraction of the standard star, or the telluric correction can be skipped in case of, for example, missing calibration files. 


### DATA REDUCED
The pipeline will generate the `data_products` directory within the working directory 
`/Users/.../working_directory` containing the reduced data, a logger file to track the information from the reduction process, and the following structure: 

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
        - `UTxxx_.spec.ap1.png`
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

`data_products` contains the `plots` directory with the output figures and quality plots obtained during the data reduction, saved in a separated for in each arm when required (`data_products/plots/arm`). 
Then, the `data_products/intermediate` directory with the calibration files generated during the data reduction process and saves the data separately for each red and blue arm. Also in `intermediate` there is the temporary directory `raw_data_temp` aimed to store the raw data and any pre-treated images (e.g. cut down calibration frames to stellar mode size, when needed) during the data reduction process. `raw_data_temp` is automatically removed when the pipeline is successfully completed. 

Finally, we find `data_products/master_calib`, which is a directory with all master calibration files produced in the data reduction. They are stored to be used in further reductions if required.


### TO DO
- Include unit tests. 

## Reporting Issues or Suggestions
If you encounter any issues or have suggestions for improving the pipeline, please [**open a new issue**](https://github.com/PyWiFeS/pipeline/issues) in the `issues` tab and fill out the provided template. Your feedback is very valuable!






