# Pipeline
The Python data reduction pipeline for WiFeS

## Development Version - use only AUTOMATION BRANCH

**Note:** This version is in active development.

### December 2023 Update
#### How is this code different from the [master](https://github.com/PyWiFeS/pipeline/tree/master) PyWiFeS repository?
- Updated to Python 3, with bug fixes and headers for the new telescope setup.
- Cleaner repository (some defunct scripts to be removed).
- Adjusted to read and reduce data automatically all at once.
- New `.JSON` config files provided for each observing mode (`classic` and `nod-and-shuffle`).
  - The pipeline chooses the template automatically.
  - Users don't need to set anything in generate metadata or reduction scripts anymore.
  - Users can create their own `.JSON` file following the same structure with their preferred setup.
- Astrometry is now implemented in the data cubes (still some small offsets in RA and Dec that will be fixed)


#### Known Problems
- **Multiprocessing**
    - **Update Mar 2024**: Multiprocessing can be enabled by so the pipeline *may* do the job faster. To enable the multithreading option, please follow these steps:
        1. Open the `.json` file that corresponds to the observing mode of your data. That is, `reduction_scripts/params_class.json` for classic mode, or `reduction_scripts/params_ns.json` for nod-and-shuffle.
        2. Set `"multithread": true` in all the cases. There should be a total of 6 `"multithread"`, 3 for each of the blue and red arms in the following steps: `"step": "wave_soln"`, `"step": "cosmic_rays"`, and `"step": "cube_gen"`.
        3. [Optional] Set `max_processes` to the *maximum* number of sub-processes you would like to launch for `"step": "cosmic_rays"`, and `"step": "cube_gen"`. If `-1`, the pipeline will use as many processes as there are hardware & logical cores on your device, which maybe be larger than the number of *available* cores, e.g. for Slurm users. Limiting the number of sub-processes can improve efficiency and availability of your device.
        4. Run the pipeline following the instructions below.

- **Skip steps**: not all the steps in the data reduction process can be skipped (i.e., `"run": false` in the `.json` files). For example, the wavelength solution is always required for a succesfull data reduction. Other steps such as the flux calibration, the extraction of the standar star, or the telliric correction can be skipped in case of missing calibration files. 


---

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
4.  Set up an alias for the main reduction routine: 
    ```sh
    alias pywifes-reduce='/Users/.../pipeline/reduction_scripts/reduce_data.py'
    ```    

## Running the Pipeline
1. Put all raw data and calibration files in a dedicated directory, e.g. `/Users/.../working_directory/my_raw_data`
2. Run the main reduction routine, giving the raw data directory path as an input parameter. The pipeline will run both arms automatically and choose the observing mode by checking the headers.
    ```sh
   pywifes-reduce my_raw_data
   ```


**DATA REDUCED:**
The pipeline will generate the `data_products` directory within the working directory 
`/Users/.../working_directory` containing the reduced data and following the structure above: 

- data_products
    - `xxx-Blue-UTxxx.cube.fits`
    - `xxx-Red--UTxxx.cube.fits`
    - `xxx-Splice-UTxxx.cube.fits`
    - `UTxxx_detected_apertures_plot.pdf` 
    - ... 
    - `xxx-Blue-UTxxx.spec.apx.fits`
    - `xxx-Red--UTxxx.spec.apx.fits`
    - `xxx-Splice-UTxxx.spec.apx.fits`
    -
    - intermediate
        - blue
            - `wifes_blue_<master_calibration>_files.fits`
            - `xxx-Blue-UTxxx.p00.fits`
            - `xxx-Blue-UTxxx.p01.fits`
            - ...
            - `blue_file_name.p10.fits`
        - red
            - `wifes_red_<master_calibration>_files.fits`
            - `xxx-Red--UTxxx.p00.fits`
            - `xxx-Red--UTxxx.p01.fits`
            - ...
            - `xxx-Red--UTxxx.p10.fits`

`data_products` also contains the `intermediate` directory with the files generated during the data reduction process that are saved separately for each red and blue arm. Those intermediate files are master calibration files (e.g., master bias, master flats, ...) and other calibration files generated during the data reduction named as `xxx-Blue/Red-UTxxx.p00.fits, xxx-Blue/Red-UTxxx.p01.fits, ..., xxx-Blue/Red-UTxxx.p10.fits`. 


### TO DO

## Reporting Issues or Suggestions
If you encounter any issues or have suggestions for improving the pipeline, please [**open a new issue**](https://github.com/PyWiFeS/pipeline/issues) in the `issues` tab and fill out the provided template. Your feedback is very valuable!






