## Data reduction

### Generate and save metadata
Run the `pipeline/reduction_scripts/extra_tools/generate_metadata_script.py` either in the directory where your data is, or passing the data directory as a command line argument like `./generate_metadata_script.py /home/mjc/wifes_data/20120922/` etc. This will generate the `save_red_metadata.py` and `save_blue_metadata.py` scripts.
Run these from the command line as `python save_red_metadata.py` etc. or edit if the automated metadata structure is not what you desire for your reduction script. These scripts will create a `my_metadata.pkl` file SOMEWHERE.

### Reduce data
#### Way 1 [`reduce_all_data.py` doesn't exist anymore?]
Copy the reduction script `reduce_all_data.py` from the `reduction_scripts/` directory and edit the steps if desired.  Run the reduction script like `./reduce_all_data.py my_metadata.pkl`, this will open the metadata file and run your desired reduction steps.
#### Way 2
Edit `reduce_red_data.py` and `reduce_blue_data.py` and reduce data with `python reduce_red_data.py my_metadata.pkl`.
