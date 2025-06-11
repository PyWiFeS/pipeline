import os
from pywifes import pywifes
from pywifes.wifes_utils import get_full_obs_list, wifes_recipe


# ------------------------------------------------------------------------
# Overscan subtraction
# ------------------------------------------------------------------------
@wifes_recipe
def _run_overscan_sub(metadata, gargs, prev_suffix, curr_suffix, poly_high_oscan=True, **args):
    """
    Subtract overscan from the input FITS files and save the results.

    Parameters
    ----------
    metadata : dict
        A dictionary containing metadata information of the FITS files.
    gargs : dict
        A dictionary containing global arguments used by the processing steps.
    prev_suffix : str
        The suffix of the previous FITS files.
    curr_suffix : str
        The suffix of the current FITS files.
    poly_high_oscan : bool
        Whether to fit the overscan with a polynomial that excludes high-count
        rows that can suffer from elevated overscan levels.
        Default: True.

    Optional Function Arguments
    ---------------------------
    detector_regions : list
        Override epoch-based detector regions with the specified values ([ymin, ymax, xmin, xmax],
        where the max values indicate the last pixel to be included).
        Default: None.
    overscan_regions : list
        Override epoch-based overscan regions with the specified values ([ymin, ymax, xmin, xmax],
        where the max values indicate the last pixel to be included).
        Default: None.
    science_regions : list
        Override epoch-based science regions with the specified values ([ymin, ymax, xmin, xmax],
        where the max values indicate the last pixel to be included).
        Default: None.
    gain : float
        Override epoch-based gain value. Units: e-/ADU.
        Default: None.
    rdnoise : float
        Override epoch-based read noise value. Units: e-.
        Default: None.
    omaskfile : str
        If 'poly_high_oscan'=True, filename of the maskfile defining the slice/interslice regions.
        Default: None, but provided automatically if `poly_high_oscan'=True.
    omask_threshold : float
        If 'poly_high_oscan'=True, threshold in per-row mean ADU relative to row with lowest mean,
        to determine whether overscan is masked.
        Default: 500.
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

    Notes
    -----
    Must override all four of (detector_regions, overscan_regions, gain, and rdnoise)
    to avoid using the epoch-based values.
    """
    full_obs_list = get_full_obs_list(metadata)
    first = True
    if not poly_high_oscan:
        first = False
        oscanmask = None
    for fn in full_obs_list:
        in_fn = os.path.join(gargs['data_dir'], "%s.fits" % fn)
        out_fn = os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, curr_suffix))
        if gargs['skip_done'] and os.path.isfile(out_fn):
            # cannot check mtime here because of fresh copy to raw_data_temp
            continue
        print(f"Subtracting Overscan for {os.path.basename(in_fn)}")
        if first:
            # Find a domeflat to generate mask for overscan
            if metadata["domeflat"]:
                dflat = os.path.join(gargs['data_dir'], "%s.fits" % metadata["domeflat"][0])
                pywifes.make_overscan_mask(dflat, omask=gargs['overscanmask_fn'], data_hdu=0)
                oscanmask = gargs['overscanmask_fn']
            else:
                oscanmask = None
            first = False
        pywifes.subtract_overscan(in_fn, out_fn, data_hdu=gargs['my_data_hdu'], omaskfile=oscanmask, **args)
    return
