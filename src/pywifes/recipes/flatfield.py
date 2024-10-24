import os
from pywifes import pywifes
from pywifes.wifes_utils import * 
# ------------------------------------------------------
# Flatfield: Division
# ------------------------------------------------------
@wifes_recipe
def _run_flatfield(metadata, gargs, prev_suffix, curr_suffix):
    """
    Apply flat-field correction (division) to science and standard frames.

    Parameters
    ----------
    metadata : dict
        Metadata containing information about the FITS files of the observations.
    gargs : dict
        A dictionary containing global arguments used by the processing steps. 
    prev_suffix : str
        Previous suffix of the file names (input).
    curr_suffix : str
        Current suffix of the file names (output).

    Returns
    -------
    None
    """
    sci_obs_list = get_primary_sci_obs_list(metadata)
    std_obs_list = get_primary_std_obs_list(metadata)
    info_print(f"Primary science observation list: {sci_obs_list}")
    info_print(f"Primary standard observation list: {std_obs_list}")
    for fn in sci_obs_list + std_obs_list:
        in_fn = os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, prev_suffix))
        out_fn = os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, curr_suffix))
        if gargs['skip_done'] and os.path.isfile(out_fn) \
                and os.path.getmtime(in_fn) < os.path.getmtime(out_fn) \
                and os.path.getmtime(gargs['flat_resp_fn']) < os.path.getmtime(out_fn):
            continue
        info_print(f"Flat-fielding image {os.path.basename(in_fn)}")
        pywifes.imarith_mef(in_fn, "/", gargs['flat_resp_fn'], out_fn)
        ffh = pyfits.getheader(gargs['flat_resp_fn'])
        try:
            ffin = ffh["PYWRESIN"]
        except:
            ffin = "Unknown"
        of = pyfits.open(out_fn, mode="update")
        of[0].header.set("PYWRESIN", ffin, "PyWiFeS: flatfield inputs")
        of.close()
    return
