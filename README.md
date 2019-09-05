# pipeline
The Python data reduction pipeline for WiFeS

# Modified pipeline to reduce young stars
## Prepare files
`python generate_metadata_script_marusa.py config.py /priv/mulga1/marusa/2m3data/20190314`
What is the output? A new folder is created:
`/priv/mulga1/marusa/2m3reduced/20190114/reduced_r_ys` and `ys_metadata_WiFeSRed.py` is stored there. Review this file to check if everything is OK.
**NOTE**: a new subfolder `wifes`
## Reduce each night
`marusa@mash:~/reduce_young_stars/logs>python reduce_red_data_python3.py config.py /priv/mulga1/marusa/2m3data/20190314`
## Convert to ascii
## Wavelength calibration
TODO

# TODO
- Blue
- Simplify:
  - Copy all necessary files into `pipeline`
  - Make it more simple and robust.
