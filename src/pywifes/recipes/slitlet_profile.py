import os
from pywifes import pywifes
from pywifes.wifes_utils import wifes_recipe


# ------------------------------------------------------
# Fit slitlet profiles
# ------------------------------------------------------
@wifes_recipe
def _run_slitlet_profile(metadata, gargs, prev_suffix, curr_suffix, **args):
    """
    Fit the slitlet profiles to the flatfield.

    Parameters
    ----------
    metadata : dict
        Metadata information.
    gargs : dict
        A dictionary containing global arguments used by the processing steps. 
    prev_suffix : str
        Previous suffix.
    curr_suffix : str
        Current suffix.

    Optional Function Arguments
    ---------------------------
    shift_global : bool
        Whether to shift all slitlets by a single global value: the mean of the
        derived per-slitlet shifts.
        Default: True.
    bin_x : int
        If specified, override the x-axis binning defined in the header.
        Default: None.
    bin_y : int
        If specified, override the y-axis binning defined in the header.
        Default: None.
    interactive_plot : bool
        Whether to interrupt processing to provide interactive plot to user.
        Default: False.
    verbose : bool
        Whether to output extra messages.
        Default: False.
    debug : bool
        Whether to report the parameters used in this function call.
        Default: False.

    Returns
    -------
    None
        This function does not return anything.

    Notes
    -----
    If the interslice-cleaned super_dflat_fn file exists, it uses it as the
    flatfield_fn. Otherwise, it uses super_dflat_raw as the flatfield_fn.

    The slitlet profiles are derived using the flatfield_fn and saved to the
    output_fn.
    """
    output_fn = gargs['slitlet_def_fn']
    if os.path.isfile(gargs['super_dflat_fn']):
        flatfield_fn = gargs['super_dflat_fn']
    else:
        flatfield_fn = gargs['super_dflat_raw']
    if gargs['skip_done'] and os.path.isfile(output_fn) \
            and os.path.getmtime(flatfield_fn) < os.path.getmtime(output_fn):
        return
    pywifes.derive_slitlet_profiles(
        flatfield_fn, output_fn, data_hdu=gargs['my_data_hdu'], **args
    )
    return
