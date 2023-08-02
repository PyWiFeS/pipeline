# =============================================================================
# 1) Poner todos los datos en el mismo directorio. En este caso /Users/felipe/2p3m-telescope/Chris-data/raw_data
# =============================================================================


# =============================================================================
# 2) Correr el script  /Users/felipe/2p3m-telescope/pipeline/generate_metadata_script.py
# proporcionando como parámetro el directorio con todos los datos 
# =============================================================================
python3 /Users/felipe/2p3m-telescope/pipeline/generate_metadata_script.py raw_data 



# =============================================================================
# 3) Run the generated save_*_metadata.py to produce dictionary (wifesB_20230312_metadata.pkl and wifesR_20230312_metadata.pkl)
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

cp /Users/felipe/2p3m-telescope/pipeline/reduction_scripts/reduce_* .

python3 reduce_red_data.py wifesR_20230312_metadata.pkl

python3 reduce_blue_data.py wifesB_20230312_metadata.pkl







# =============================================================================
# Extract spectrum from Cube
# =============================================================================

python DEbass_extractSpec.py --redArm PATH_TO_Red*.p11.fits --blueArm PATH_TO_Blue*.p11.fits --skySub 





# =============================================================================
# Splice blue and red spectra
# =============================================================================


python3 /Users/felipe/2p3m-telescope/slices/python3/DEbass_spliceSpec.py --blueArm OBK-282624-WiFeS-Blue-UT20230312T111056-9.p12_SN.fits --redArm /OBK-282624-WiFeS-Red--UT20230312T111056-9.p12_SN.fits --scale



python DEbass_spliceSpec.py --blueArm OBK-282624-WiFeS-Red--UT20230312T111056-9.p12_SN.fits --redArm OBK-282624-WiFeS-Blue-UT20230312T111056-9.p12_SN.fits --scale










