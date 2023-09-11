#!/bin/bash
#
#SBATCH --job-name=test-pywife-reduce
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=12g
#SBATCH --cpus-per-task=8

source ./setup_env.sh

export run_dir=$1

if [ $# -ne 2 ]
  then
    echo "Arguments for a run directory name and whether to enable multithreading is required."
    exit 1
fi

if [ ! -d $1 ] 
  then
    mkdir $run_dir
fi

pushd $run_dir

is_multithreaded=$2

echo $is_multithreaded

# =============================================================================
# 1) Put all the data in the same directory. /fred/oz016/tdavies/artifacts/optical_2/raw_data
# =============================================================================
export pipeline_prefix=/fred/oz016/tdavies/projects/optical/PyWiFeS/pipeline
#export data_dir=/fred/oz016/tdavies/artifacts/optical_2/raw_data
export data_dir=/fred/oz016/tdavies/projects/optical/run1/raw_data

# =============================================================================
# 2) Run the script  /fred/oz016/tdavies/projects/optical/PyWiFeS/pipeline/generate_metadata_script.py
# giving the data directory as an input parameter
# =============================================================================

python3 $pipeline_prefix/generate_metadata_script.py $data_dir


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

cp $pipeline_prefix/reduction_scripts/reduce_* .

python3 reduce_red_data.py wifesR_20230312_metadata.pkl $is_multithreaded

python3 reduce_blue_data.py wifesB_20230312_metadata.pkl $is_multithreaded


# DATA REDUCED !!!!

# =============================================================================
# Extract spectrum from Cube
# =============================================================================

python3 DEbass/DEbass_extractSpec.py --redArm reduc_r/*.p11.fits --blueArm reduc_b/*.p11.fits --skySub 


# =============================================================================
# Splice blue and red spectra
# =============================================================================

# python3 DEbass/DEbass_spliceSpec.py --blueArm OBK-282624-WiFeS-Blue-UT20230312T111056-9.p12_SN.fits --redArm /OBK-282624-WiFeS-Red--UT20230312T111056-9.p12_SN.fits --scale

python3 DEbass/DEbass_spliceSpec.py --blueArm OBK-282624-WiFeS-Red--UT20230312T111056-9.p12_SN.fits --redArm OBK-282624-WiFeS-Blue-UT20230312T111056-9.p12_SN.fits --scale

popd
