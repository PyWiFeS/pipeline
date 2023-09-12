#!/bin/bash

# =============================================================================
# 1) Put all the data in the same directory. /Users/.../raw_data
# =============================================================================


# =============================================================================
# 2) Run the script  /Users/.../pipeline/generate_metadata_script.py
# giving the data directory as an input parameter
# =============================================================================
python3 /Users/.../pipeline/generate_metadata_script.py raw_data 



# =============================================================================
# 3) Run the generated save_*_metadata.py to produce the dictionary (wifesB_20230312_metadata.pkl and wifesR_20230312_metadata.pkl)
# =============================================================================

python3 save_blue_metadata.py 
python3 save_red_metadata.py 




# =============================================================================
# 4) Create dir reduc_r and mkdir reduc_b
# =============================================================================
mkdir reduc_r reduc_b


# =============================================================================
# 5) Copy reduce data scripts
# check the reduce_*_data.py 
# then run them using as a parameter the gererated .pkl dictionary  (with the corresponding color)

# =============================================================================

cp /Users/.../pipeline/reduction_scripts/reduce_* .

python3 reduce_red_data.py wifesR_20230312_metadata.pkl

python3 reduce_blue_data.py wifesB_20230312_metadata.pkl




# DATA REDUCED !!!!





# =============================================================================
# Extract spectrum from Cube
# =============================================================================

python DEbass_extractSpec.py --redArm PATH_TO_Red*.p11.fits --blueArm PATH_TO_Blue*.p11.fits --skySub 




# =============================================================================
# Splice blue and red spectra
# =============================================================================

python3 /Users/.../slices/python3/DEbass_spliceSpec.py --blueArm OBK-282624-WiFeS-Blue-UT20230312T111056-9.p12_SN.fits --redArm /OBK-282624-WiFeS-Red--UT20230312T111056-9.p12_SN.fits --scale

python DEbass_spliceSpec.py --blueArm OBK-282624-WiFeS-Red--UT20230312T111056-9.p12_SN.fits --redArm OBK-282624-WiFeS-Blue-UT20230312T111056-9.p12_SN.fits --scale










