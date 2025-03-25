import os
from pywifes import pywifes
from pywifes.wifes_utils import wifes_recipe


# ------------------------------------------------------
# Flat cleanup
# ------------------------------------------------------
@wifes_recipe
def _run_flat_cleanup(
    metadata,
    gargs,
    prev_suffix,
    curr_suffix,
    type=["dome", "twi"],
    offsets=[0.0, 0.0],
    **args,
):
    """
    Make the master domeflat and twilight flat corrections.

    Parameters
    ----------
    metadata : str
        The metadata information for the flat cleanup.
    gargs : dict
        A dictionary containing global arguments used by the processing steps.
    prev_suffix : str
        The previous suffix.
    curr_suffix : str
        The current suffix.
    type : list, optional
        The types of flats to correct.
        Options: List containing one or more of 'dome', 'twi'
        Default: ['dome', 'twi'].
    offsets : list, optional
        The offsets for each specified type of flat.
        Default: [0.0, 0.0].

    Optional Function Arguments
    ---------------------------
    bin_x : int
        If specified, override the x-axis binning defined in the header.
        Default: None.
    bin_y : int
        If specified, override the y-axis binning defined in the header.
        Default: None.
    data_hdu : int
        Number of primary HDU in file.
        Default: 0.
    buffer : int
        Number of y-axis pixels around each slitlet definition to exclude when
        calculating interslice background.
        Default: 0.
    radius : int
        Sigma for 2D Gaussian filtering (with 'method'='2D'). Equal in both axes.
        Default: 10.
    nsig_lim : float
        Clipping threshold in number of standard deviations above median. Computed
        independently for each slitlet.
        Default: 5.
    plot : bool
        Whether to output a diagnostic plot.
        Default: False.
    save_prefix : str
        Prefix for plot (if requested).
        Default: 'cleanup_'.
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
    """
    # check the slitlet definition file
    if os.path.isfile(gargs['slitlet_def_fn']):
        slitlet_fn = gargs['slitlet_def_fn']
    else:
        slitlet_fn = None
    if "dome" in type:
        if gargs['skip_done'] and os.path.isfile(gargs['super_dflat_fn']) \
                and os.path.getmtime(gargs['super_dflat_raw']) < os.path.getmtime(gargs['super_dflat_fn']):
            return
        if os.path.isfile(gargs['super_dflat_raw']):
            print(f"Correcting master domeflat {os.path.basename(gargs['super_dflat_raw'])}")
            pywifes.interslice_cleanup(
                gargs['super_dflat_raw'],
                gargs['super_dflat_fn'],
                slitlet_def_file=slitlet_fn,
                offset=offsets[type.index("dome")],
                method="2D",
                plot_dir=gargs['plot_dir_arm'],
                **args,
            )
        else:
            print(f"Master dome flat {os.path.basename(gargs['super_dflat_raw'])} "
                  "not found. Skipping dome flat cleanup.")

    if "twi" in type:
        if gargs['skip_done'] and os.path.isfile(gargs['super_tflat_fn']) \
                and os.path.getmtime(gargs['super_tflat_raw']) < os.path.getmtime(gargs['super_tflat_fn']):
            return
        if os.path.isfile(gargs['super_tflat_raw']):
            print(f"Correcting master twilight flat {os.path.basename(gargs['super_tflat_raw'])}")
            pywifes.interslice_cleanup(
                gargs['super_tflat_raw'],
                gargs['super_tflat_fn'],
                slitlet_def_file=slitlet_fn,
                offset=offsets[type.index("twi")],
                method="2D",
                plot_dir=gargs['plot_dir_arm'],
                **args,
            )
        else:
            print(f"Master twilight flat {os.path.basename(gargs['super_tflat_raw'])} "
                  "not found. Skipping twilight flat cleanup.")
    return
