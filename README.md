# Pipeline
The Python data reduction pipeline for WiFeS

## Development Version - use only AUTOMATION BRANCH

**Note:** This version is in active development.

### November 2022 Update
#### How is this code different from the main pyWiFeS repository?
- Updated to Python 3, with bug fixes and headers for the new telescope setup.
- Cleaner repository (some defunct scripts to be removed).
- Adjusted to read and reduce data automatically all at once.
- New `.JSON` config files provided for each observing mode (`classic` and `nod-and-shuffle`).
  - The pipeline chooses the template automatically.
  - Users don't need to set anything in generate metadata or reduction scripts anymore.
  - Users can create their own `.JSON` file following the same structure with their preferred setup.

#### Known Problems
- Multiprocessing not always works properly.

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

## Running the Pipeline
1. Put all raw data and calibration files in the same directory: `/Users/.../my_folder/raw_data`
2. Copy the reduce data script and `.json` files to the above-mentioned folder:
    ```sh
   cp /Users/.../pipeline/reduction_scripts/reduce_data.py /Users/.../my_folder/
   cp /Users/.../pipeline/reduction_scripts/*.json /Users/.../my_folder/
   ```
3. Run `reduce_data.py`, giving the raw data directory as an input parameter. The pipeline will run both arms automatically and choose the observing mode by checking the headers.
    ```sh
   python3 reduce_data.py raw_data
   ```



**DATA REDUCED!**

### TO DO
- Extract spectrum from Cube
- Splice blue and red spectra

## Reporting Issues or Suggestions
If you encounter any issues or have suggestions for improving the pipeline, please [**open a new issue**](https://github.com/PyWiFeS/pipeline/issues) in the `issues` tab and fill out the provided template. Your feedback is valuable!






