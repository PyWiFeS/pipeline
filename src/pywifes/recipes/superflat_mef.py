import os
from pywifes import pywifes
from pywifes.wifes_utils import wifes_recipe


# ------------------------------------------------------
# Create MEF files
# ------------------------------------------------------
@wifes_recipe
def _run_superflat_mef(metadata, gargs, prev_suffix, curr_suffix, source, **args):
    """
    Generate a Multi-Extension FITS (MEF) file for the calibration files.

    Parameters
    ----------
    metadata : str
        The metadata information.
    gargs : dict
        A dictionary containing global arguments used by the processing steps.
    prev_suffix : str
        The previous suffix.
    curr_suffix : str
        The current suffix.
    source : str
        The source of the flatfield.
        Options: 'dome', 'twi', 'wire', 'arc'.

    Optional Function Arguments
    ---------------------------
    nan_method : str
        Method for treating NaN pixels: linear interpolation in the x-direction,
        or replace with a constant value.
        Options: 'interp', 'replace'.
        Default: 'interp'.
    repl_val : float
        If 'nan_method'='replace', replace NaN pixels with the value of 'repl_val'.
        Default: 0.
    bin_x : int
        If specified, override the x-axis binning defined in the header.
        Default: None.
    bin_y : int
        If specified, override the y-axis binning defined in the header.
        Default: None.
    debug : bool
        Whether to report the parameters used in this function call.
        Default: False.

    Returns
    -------
    None
        This function does not return anything.

    Raises
    ------
    ValueError
        If the flatfield type is not recognized.

    Notes
    -----
    This function generates a Multi-Extension FITS (MEF) file for a calibration
    file. It checks for the existence of the master file and generates the MEF file
    accordingly. If the master flat file is not found, it prints a warning message
    and skips the MEF generation. The slitlet definition file is also checked and
    used if available. Finally, the MEF file is generated using the
    `pywifes.wifes_slitlet_mef` function.
    """
    if source == "dome":
        if os.path.isfile(gargs['super_dflat_fn']):
            in_fn = gargs['super_dflat_fn']
        elif os.path.isfile(gargs['super_dflat_raw']):
            in_fn = gargs['super_dflat_raw']
        else:
            print("No master dome flat found. Skipping MEF generation for dome flat.")
            return
        out_fn = gargs['super_dflat_mef']

    elif source == "twi":
        if os.path.isfile(gargs['super_tflat_fn']):
            in_fn = gargs['super_tflat_fn']
        elif os.path.isfile(gargs['super_tflat_raw']):
            in_fn = gargs['super_tflat_raw']
        else:
            print("No master twilight flat found. Skipping MEF generation for twilight flat.")
            return
        out_fn = gargs['super_tflat_mef']

    elif source == "wire":
        if os.path.isfile(gargs['super_wire_raw']):
            in_fn = gargs['super_wire_raw']
        else:
            print("No master wire frame found. Skipping MEF generation for wire.")
            return
        out_fn = gargs['super_wire_mef']

    elif source == "arc":
        if os.path.isfile(gargs['super_arc_raw']):
            in_fn = gargs['super_arc_raw']
        else:
            print("No master arc frame found. Skipping MEF generation for arc.")
            return
        out_fn = gargs['super_arc_mef']

    else:
        print(f"Calibration type '{source}' not recognised")
        raise ValueError(f"Calibration type '{source}'' not recognised")

    # check the slitlet definition file
    if os.path.isfile(gargs['slitlet_def_fn']):
        slitlet_fn = gargs['slitlet_def_fn']
    else:
        slitlet_fn = None
    # run it!
    if gargs['skip_done'] and os.path.isfile(out_fn) \
            and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
        return
    print(f"Generating MEF {source} flat")
    pywifes.wifes_slitlet_mef(
        in_fn, out_fn, data_hdu=gargs['my_data_hdu'], slitlet_def_file=slitlet_fn, **args
    )
    return
