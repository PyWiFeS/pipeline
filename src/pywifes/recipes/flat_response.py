import os
from pywifes import pywifes
from pywifes.wifes_utils import * 
# ------------------------------------------------------
# Flatfield: Response
# ------------------------------------------------------
@wifes_recipe
def _run_flat_response(metadata, gargs, prev_suffix, curr_suffix, mode="all", **args):
    """
    Generate the flatfield response function for each flat type, either dome or
    both (dome and twi).

    Parameters
    ----------
    metadata : dict
        Metadata containing information about the data FITS files.
    gargs : dict
        A dictionary containing global arguments used by the processing steps. 
    prev_suffix : str
        Suffix of the previous data.
    curr_suffix : str
        Suffix of the current data.
    mode : str, optional
        Mode for generating the response function.
        Options: 'dome', 'all'.
        Default: 'all'.

    Optional Function Arguments
    ---------------------------
    zero_var : bool
        Whether to set the VAR extension of the output response file to zero.
        Default: True.
    polydeg : int
        Order of polynomial used in fitting wavelength profile of dome flat (in
        log-flux).
    resp_min : float
        Minimum value allowed in final response function.
        Default: 0.0001.
    save_prefix : str
        Prefix for the diagnostic plot.
        Default: 'flat_response'.
    interactive_plot : bool
        Whether to interrupt processing to provide interactive plot to user.
        Default: False.
    debug : bool
        Whether to report the parameters used in this function call.
        Default: False.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the requested response mode is not recognized.
    """
    # Fit the desired style of response function
    if gargs['skip_done'] and os.path.isfile(gargs['flat_resp_fn']) \
            and os.path.getmtime(gargs['super_dflat_mef']) < os.path.getmtime(gargs['flat_resp_fn']):
        return
    info_print("Generating flatfield response function")
    if mode == "all" and not os.path.isfile(gargs['super_tflat_mef']):
        info_print("WARNING: No twilight superflat MEF found. Falling back to dome flat only.")
        mode = "dome"

    if mode == "all":
        pywifes.wifes_SG_response(
            gargs['super_dflat_mef'],
            gargs['super_tflat_mef'],
            gargs['flat_resp_fn'],
            wsol_fn=gargs['wsol_out_fn'],
            plot_dir=gargs['plot_dir_arm'],
            shape_fn=gargs['smooth_shape_fn'],
            **args
        )
    elif mode == "dome":
        pywifes.wifes_response_poly(
            gargs['super_dflat_mef'], gargs['flat_resp_fn'], wsol_fn=gargs['wsol_out_fn'], shape_fn=gargs['smooth_shape_fn'], **args
        )
    else:
        error_print("Requested response mode not recognised")
        raise ValueError("Requested response mode not recognised")
    return
