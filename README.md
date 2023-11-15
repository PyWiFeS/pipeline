## This is a version IN DEVELOPMENT 
## This info is UPDATED for the AUTOMATION BRANCH ONLY



# Pipeline
The Python data reduction pipeline for WiFeS


### NOV 2022: How is this code different from the main pyWiFeS repository?
- Updated to python3, bugs fixed, updated the headers for the new telescope setup, cleaner repository (but some defunct scripts still to be removed)
- This version is adjusted to read and reduce the data automatically all at once.  
- New `.JSON` config files: there are two templates provided to set the data reduction configuration parameters, one for each observing mode (`classic` and `nod-and-shuffle`). The pipeline chooses the right template automatically. You the user don't have to set anything in generate metadata or reduction scripts anymore. The user can make their own `.JSON` file following the same structure with their preferred setup (and update the `reduce_data.py` accordingly so it will load it).


### Known problems
- Multiprocessing not always work properly


## How to install the pipeline
- Download or clone the `pipeline` repository in the `automation` branch
- In a Python environment and from the pipeline root directory, do:
`pip install .`


## How to run the pipeline
- Put all the raw data and calibration files in the same directory 
`/Users/.../my_folder/raw_data`
- Copy the reduce data script and .json files to the above mentioned folder
`cp /Users/.../pipeline/reduction_scripts/reduce_data.py /Users/.../my_folder/`
`cp /Users/.../pipeline/reduction_scripts/*.json /Users/.../my_folder/`
- Run `reduce_data.py` giving the raw data directory as an input parameter. The pipeline will run both arms automatically and choose the observing mode checking the headers.
`python3 reduce_data.py raw_data`


### DATA REDUCED !!!!



# TO DO   
- Extract spectrum from Cube
- Splice blue and red spectra









