import os
from pywifes import pywifes
from pywifes.wifes_utils import (
    get_sci_obs_list, get_std_obs_list, is_nodshuffle, is_subnodshuffle, set_header,
    wifes_recipe
)


# ------------------------------------------------------
# Sky subtraction
# ------------------------------------------------------
# since run_sky_sub_ns is called from run_sky_sub,
# no need to have decorator here
# @wifes_recipe
def _run_sky_sub_ns(metadata, gargs, prev_suffix, curr_suffix, separate_ns=False):
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

    Optional Function Arguments
    ---------------------------
    separate_ns : bool
        Whether to separate the two N&S positions and progress them independently.
        Default: False.

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
        if is_nodshuffle(in_fn) or is_subnodshuffle(in_fn):
            if separate_ns:
                if gargs['skip_done'] and os.path.isfile(out_fn) \
                        and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                    pass
                else:
                    print(f"Copying N&S science image {os.path.basename(in_fn)}")
                    pywifes.imcopy(in_fn, out_fn)

                out_sky_fn = os.path.join(gargs['out_dir'], "%ss.p%s.fits" % (fn, curr_suffix))
                if gargs['skip_done'] and os.path.isfile(out_sky_fn) \
                        and os.path.getmtime(sky_fn) < os.path.getmtime(out_sky_fn):
                    pass
                else:
                    print(f"Copying N&S sky image {os.path.basename(sky_fn)}")
                    pywifes.imcopy(sky_fn, out_sky_fn)
                    # Add 'sky' from N&S to science metadata
                    metadata['sci'].append({"sci": [os.path.basename(out_sky_fn).replace(f".p{curr_suffix}.fits", "")], "sky": []})
                    # Remove N&S label
                    set_header(out_fn, "WIFESOBS", "ClassicalEqual")
                    set_header(out_sky_fn, "WIFESOBS", "ClassicalEqual")
            else:
                if gargs['skip_done'] and os.path.isfile(out_fn) \
                        and os.path.getmtime(in_fn) < os.path.getmtime(out_fn) \
                        and os.path.getmtime(sky_fn) < os.path.getmtime(out_fn):
                    continue
                print(f"Subtracting N+S sky frame for {os.path.basename(in_fn)}")
                pywifes.scaled_imarith_mef(in_fn, "-", sky_fn, out_fn, scale="exptime")
        else:
            if gargs['skip_done'] and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                continue
            print(f"Copying non-N&S image {os.path.basename(in_fn)}")
            pywifes.imcopy(in_fn, out_fn)
    return


@wifes_recipe
def _run_sky_sub(metadata, gargs, prev_suffix, curr_suffix, **args):
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
        _run_sky_sub_ns(metadata, gargs, prev_suffix, curr_suffix, **args)
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
                    print(f"Coadding sky frames into {os.path.basename(sky_proc_fn)}")
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
                    print(f"Subtracting sky frame for {os.path.basename(in_fn)}")
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
                    print(f"Copying science image {os.path.basename(in_fn)}")
                    pywifes.imcopy(in_fn, out_fn)
        # copy stdstar frames
        std_obs_list = get_std_obs_list(metadata)
        for fn in std_obs_list:
            in_fn = os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(gargs['out_dir'], "%s.p%s.fits" % (fn, curr_suffix))
            if gargs['skip_done'] and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                continue
            print(f"Copying standard star image {os.path.basename(in_fn)}")
            pywifes.imcopy(in_fn, out_fn)
    return
