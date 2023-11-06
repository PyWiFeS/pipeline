# =============================================================================
# 1) Put all the raw data and calibration files in the same directory 
# /Users/.../my_folder/raw_data
# =============================================================================

# =============================================================================
# 2) Copy the reduce data script and .json files to the above menctioned folder
# =============================================================================

cp /Users/.../pipeline/reduction_scripts/reduce_data.py /Users/.../my_folder/
cp /Users/.../pipeline/reduction_scripts/*.json /Users/.../my_folder/


# =============================================================================
# 3) Run reduce_data.py giving the data directory as an input parameter
# =============================================================================

python3 reduce_data.py raw_data


# DATA REDUCED !!!!








###########
# TO DO   #
###########

# =============================================================================
# Extract spectrum from Cube
# =============================================================================

python DEbass_extractSpec.py --redArm PATH_TO_Red*.p11.fits --blueArm PATH_TO_Blue*.p11.fits --skySub 




# =============================================================================
# Splice blue and red spectra
# =============================================================================

python3 /Users/.../slices/python3/DEbass_spliceSpec.py --blueArm OBK-282624-WiFeS-Blue-UT20230312T111056-9.p12_SN.fits --redArm /OBK-282624-WiFeS-Red--UT20230312T111056-9.p12_SN.fits --scale

python DEbass_spliceSpec.py --blueArm OBK-282624-WiFeS-Red--UT20230312T111056-9.p12_SN.fits --redArm OBK-282624-WiFeS-Blue-UT20230312T111056-9.p12_SN.fits --scale










