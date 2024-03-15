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
        1. Ensure you are using Python 3.10 or later. Below Python 3.10, you may encounter warnings or errors.
        2. Open the `.json` file that corresponds to the observing mode of your data. That is, `reduction_scripts/params_class.json` for classic mode, or `reduction_scripts/params_ns.json` for nod-and-shuffle.
        3. Set `"multithread": true` in all the cases. There should be a total of 6 `"multithread"`, 3 for each of the blue and red arms in the following steps: `"step": "wave_soln"`, `"step": "cosmic_rays"`, and `"step": "cube_gen"`.
        4. [Optional] Set `max_processes` to the *maximum* number of sub-processes you would like to launch for `"step": "cosmic_rays"`, and `"step": "cube_gen"`. If `-1`, the pipeline will use as many processes as there are hardware & logical cores on your device, which maybe be larger than the number of *available* cores, e.g. for Slurm users. Limiting the number of sub-processes can improve efficiency and availability of your device.
        5. Run the pipeline following the instructions below.

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
1. Put all raw data and calibration files in the same directory: `/Users/.../my_directory/my_raw_data`
2. Run the main reduction routine, giving the raw data directory path as an input parameter. The pipeline will run both arms automatically and choose the observing mode by checking the headers.
    ```sh
   pywifes-reduce my_raw_data
   ```


**DATA REDUCED:**
The pipeline will automatically generate the following directory structure containing the reduced data within our working directory `/Users/.../my_directory`: 

- data_products
    - `red_file_name.cube.fits`
    - `blue_file_name.cube.fits`
    - intermediate
        - blue
            - `blue_master_calibration_files.fits`
            - `blue_file_name.p00.fits`
            - `blue_file_name.p01.fits`
            - ...
            - `blue_file_name.p10.fits`
        - red
            - `master_calibration_files.fits`
            - `red_file_name.p00.fits`
            - `red_file_name.p01.fits`
            - ...
            - `red_file_name.p10.fits`

The main diretory contains a subfolder for intermediate files generated during the data reduction process that are saved separately for each red and blue arm. Those intermediate files are master calibration files (e.g., master bias, master flats) and other calibration files generated during the data reduction named as `file_name.p00.fits, file_name.p01.fits, ..., file_name.p10.fits`. Initially, the final data cubes (`file_name.cube.fits`) are saved in these intermediate directories. Once the data reduction finishses its job succesfully, the data cubes are moved to `\data_products`.


### TO DO
- Extract spectra from the data cubes
- Splice blue and red spectra

## Reporting Issues or Suggestions
If you encounter any issues or have suggestions for improving the pipeline, please [**open a new issue**](https://github.com/PyWiFeS/pipeline/issues) in the `issues` tab and fill out the provided template. Your feedback is very valuable!






