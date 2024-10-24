import os
from pywifes import pywifes
from pywifes.wifes_utils import * 
# ------------------------------------------------------
# Sky subtraction
# ------------------------------------------------------
# since run_sky_sub_ns is called from run_sky_sub, 
# no need to have decorator here
# @wifes_recipe
def _run_sky_sub_ns(metadata, gargs, prev_suffix, curr_suffix):
    """
    Subtract sky frames from science objects in nod-and-shuffle mode.

    Parameters
    ----------
    metadata : dict
        Metadata containing information about the observations.
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
    sci_obs_list = get_sci_obs_list(metadata)
    std_obs_list = get_std_obs_list(metadata)
    ns_proc_list = sci_obs_list + std_obs_list
    for fn in ns_proc_list:
        in_fn = os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, prev_suffix))
        out_fn = os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, curr_suffix))
        sky_fn = os.path.join(gargs['out_dir'], "%s.s%s.fits" % (fn, prev_suffix))
        if gargs['skip_done'] and os.path.isfile(out_fn) \
                and os.path.getmtime(in_fn) < os.path.getmtime(out_fn) \
                and os.path.getmtime(sky_fn) < os.path.getmtime(out_fn):
            continue
        info_print(f"Subtracting N+S sky frame for {os.path.basename(in_fn)}")
        pywifes.scaled_imarith_mef(in_fn, "-", sky_fn, out_fn, scale="exptime")
    return

@wifes_recipe
def _run_sky_sub(metadata, gargs, prev_suffix, curr_suffix):
    """
    Subtract sky frames from science objects.

    Parameters
    ----------
    metadata : dict
        Metadata containing information about the observations.
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
    if gargs['obs_mode'] == "ns":
        _run_sky_sub_ns(metadata, gargs, prev_suffix, curr_suffix)
    else:
        # subtract sky frames from science objects
        for obs in metadata["sci"]:
            if len(obs["sky"]) > 0:
                if len(obs["sky"]) > 1:
                    # If multiple sky frames, scale by exposure time and
                    # median-combine. List will be recombined for every science
                    # frame in case the first sky (which defines the output
                    # filename) was associated with multiple science frames. Thus,
                    # no consideration of skip_done.
                    in_fn_list = [
                        os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, prev_suffix))
                        for fn in obs["sky"]
                    ]
                    sky_proc_fn = os.path.join(
                        gargs['out_dir'], "%s.p%s.fits" % (obs["sky"][0], curr_suffix)
                    )
                    info_print(f"Coadding sky frames into {os.path.basename(sky_proc_fn)}")
                    pywifes.imcombine_mef(in_fn_list, sky_proc_fn, scale="exptime", method="median")
                else:
                    sky_fn = obs["sky"][0]
                    sky_proc_fn = os.path.join(
                        gargs['out_dir'], "%s.p%s.fits" % (sky_fn, prev_suffix)
                    )
                for fn in obs["sci"]:
                    in_fn = os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, prev_suffix))
                    out_fn = os.path.join(
                        gargs['out_dir'], "%s.p%s.fits" % (fn, curr_suffix)
                    )
                    if gargs['skip_done'] and os.path.isfile(out_fn) \
                            and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                        continue
                    info_print(f"Subtracting sky frame for {os.path.basename(in_fn)}")
                    # subtract scaled sky framefrom science frame
                    pywifes.scaled_imarith_mef(
                        in_fn, "-", sky_proc_fn, out_fn, scale="exptime"
                    )
            else:
                for fn in obs["sci"]:
                    in_fn = os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, prev_suffix))
                    out_fn = os.path.join(
                        gargs['out_dir'], "%s.p%s.fits" % (fn, curr_suffix)
                    )
                    if gargs['skip_done'] and os.path.isfile(out_fn) \
                            and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                        continue
                    info_print(f"Copying science image {os.path.basename(in_fn)}")
                    pywifes.imcopy(in_fn, out_fn)
        # copy stdstar frames
        std_obs_list = get_std_obs_list(metadata)
        for fn in std_obs_list:
            in_fn = os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, curr_suffix))
            if gargs['skip_done'] and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                continue
            info_print(f"Copying standard star image {os.path.basename(in_fn)}")
            pywifes.imcopy(in_fn, out_fn)
    return
