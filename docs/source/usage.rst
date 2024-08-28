.. _usage:

Running the Pipeline
====================

Quick Start
-----------

The pipeline is designed to be user-friendly and easy to use. The main routine is called `pywifes-reduce` and is used to reduce the data. The pipeline will automatically determine the observing mode and the arm used (blue or red) by checking the headers of the raw data. The pipeline will then run the corresponding reduction script for the arm and observing mode.

1. Put all raw data and calibration files in a dedicated directory, e.g., `/Users/.../working_directory/my_raw_data`.
2. Run the main reduction routine, giving the raw data directory path as an input parameter. The pipeline will run both arms automatically and choose the observing mode by checking the headers.
    .. code-block:: bash

        pywifes-reduce /Users/.../working_directory/my_raw_data


User-Defined Reduction Parameters
---------------------------------

**Set reduction steps**

To specify the reduction steps for blue and red data, users can provide the paths to the respective JSON files using the `--red-params` and `--blue-params` flags as follows:
   .. code-block:: bash

        pywifes-reduce /Users/.../working_directory/my_raw_data --red-params /Users/.../user_red_param_file.json --blue-params /Users/.../user_blue_param_file.json

**Reduce Data Using Master Calibration Files**

To perform data reduction using master calibration files from previous reductions, use the `--from-master` flag along with the path to the directory containing all calibration files. Both blue and red data should be stored together in the same directory for coherence. If no path is specified after the `--from-master` flag, the default directory `./data_products/master_calib` is assumed.
   .. code-block:: bash
      
       pywifes-reduce /Users/.../working_directory/my_raw_data --from-master /Users/.../master_calibrations_directory

**Other parameters**

- `-skip-done`: Skip already completed steps.
- `-just-calib`: Triggers the data reduction in the absence of on-sky data (both science and calibration). It only produces basic calibration files.

Extra Usabilities
-----------------

**Multiprocessing**

When multiprocessing is enabled, the pipeline *may* do the job faster. This will depend on the operating system used to run the pipeline. The multiprocessing setup is recommended for **Linux** users, as they will see a significant improvement in computation time. On the other hand, Mac OS users might get a similar running time (or just slightly faster) than in one-process mode. 

To enable the multithreading option, please follow these steps:

1. Open the `.json` file that corresponds to your data observing mode and grating. That is, `reduction_scripts/pipeline_parms/params_class_<grating>.json` for classic mode, or `reduction_scripts/params_ns_<grating>.json` for nod-and-shuffle.
2. Set `"multithread": true` in all the cases. There should be a total of 6 `"multithread"`, 3 for each of the blue and red arms in the following steps: `"step": "wave_soln"`, `"step": "cosmic_rays"`, and `"step": "cube_gen"`.
3. [Optional] Set `max_processes` to the *maximum* number of sub-processes you would like to launch for `"step": "cosmic_rays"`, and `"step": "cube_gen"`. If `-1`, the pipeline will use as many processes as there are hardware & logical cores on your device, which may be larger than the number of *available* cores, e.g., for Slurm users. Limiting the number of sub-processes can improve the efficiency and availability of your device.
4. Run the pipeline following the instructions below.

**Skip steps** 

Some steps in the data reduction process can be skipped by setting `"run": false` in the corresponding step in the `.json` files. However, in some cases, the step cannot be skipped as it is required for the pipeline to continue reducing the data. For example, the wavelength solution is always required for a successful data reduction. Other steps such as the flux calibration, the extraction of the standard star, or the telluric correction can be skipped in case of, for example, missing calibration files.
