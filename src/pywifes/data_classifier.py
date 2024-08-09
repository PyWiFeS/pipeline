import os
from astropy.io import fits as pyfits
from itertools import chain
import json
import pandas as pd

from . import wifes_calib


def get_obs_metadata(filenames, data_dir, greedy_stds=False, coadd_mode="all", mode_save_fn=None, camera="blue"):
    """
    Retrieve metadata for observed data files.

    This function categorizes observed data files based on their image type
    ('BIAS', 'FLAT', 'SKYFLAT', 'DARK', 'ARC', 'WIRE', 'STANDARD', 'OBJECT')
    and object name. It groups standard star observations together and separates
    science observations from standard star observations.


    Parameters
    ----------
    filenames : list of str
        List of filenames of observed data files.
    data_dir : str
        Directory path where the data files are located.
    greedy_stds : bool
        Whether to treat OBJECT exposures near known standard stars as STANDARDs.
    coadd_mode : str
        Whether to group together OBJECT frames for the same target. Options:
        'all' (default), 'none', 'prompt' (user-selected).
    mode_save_fn : str
        Filename for saving grouping of images to coadd if coadd_mode = 'prompt'.
    camera : str
        Which camera, for parsing of coadd association dictionary. Allowed values: "blue", "red".

    Returns
    -------
    dict
        Dictionary containing metadata for observed data files. The dictionary
        has the following keys:

        - 'bias': List of filenames of bias frames.
        - 'domeflat': List of filenames of domeflat frames.
        - 'twiflat': List of filenames of twilight flat frames.
        - 'dark': List of filenames of dark frames.
        - 'arc': List of filenames of arc frames.
        - 'wire': List of filenames of wire frames.
        - 'sci': List of dictionaries, each containing information about science observations. Each dictionary has the following keys:
            - 'sci': List of filenames of science observations.
            - 'sky': Empty list (not used in this function).
        - 'std': List of dictionaries, each containing information about standard star observations. Each dictionary has the following keys:
            - 'sci': List of filenames of standard star observations.
            - 'name': Name of the standard star.
            - 'type': List of strings indicating the type of observation ('flux','telluric').

    """
    stdstar_list = wifes_calib.ref_fname_lookup.keys()

    # classify each obs
    bias = []
    domeflat = []
    twiflat = []
    dark = []
    arc = []
    wire = []
    stdstar = {}
    science = {}

    for filename in filenames:
        basename = filename.replace(".fits", "")

        f = pyfits.open(data_dir + filename)
        imagetype = f[0].header["IMAGETYP"].upper()
        obj_name = f[0].header["OBJECT"]
        f.close()
        # ---------------------------
        # Check if it is within a close distance to a standard star.
        # If so and greedy_stds is True, fix the object name to be the good one from the list
        try:
            near_std, std_dist = wifes_calib.find_nearest_stdstar(data_dir + filename)
            if (greedy_stds or imagetype == "STANDARD") and std_dist < 100.0:
                obj_name = near_std
        except:
            pass
        # ---------------------------
        # 1 - bias frames
        if imagetype == "BIAS" or imagetype == "ZERO":
            bias.append(basename)
        # 2 - quartz flats
        elif imagetype == "FLAT":
            domeflat.append(basename)
        # 3 - twilight flats
        elif imagetype == "SKYFLAT":
            twiflat.append(basename)
        # 4 - dark frames
        elif imagetype == "DARK":
            dark.append(basename)
        # 5 - arc frames
        elif imagetype == "ARC":
            arc.append(basename)
        # 6 - wire frames
        elif imagetype == "WIRE":
            wire.append(basename)
        # 7 - standard star
        elif imagetype == "STANDARD":
            # group standard obs together!
            if obj_name in stdstar.keys():
                stdstar[obj_name].append(basename)
            else:
                stdstar[obj_name] = [basename]

        # all else are science targets (also allow for standard stars in imagetype = OBJECT)
        elif imagetype == "OBJECT":
            if greedy_stds and obj_name in stdstar_list:
                # group standard obs together!
                if obj_name in stdstar.keys():
                    stdstar[obj_name].append(basename)
                else:
                    stdstar[obj_name] = [basename]
            else:
                # group science obs together, if desired
                if coadd_mode == "all":
                    if obj_name in science.keys():
                        science[obj_name].append(basename)
                    else:
                        science[obj_name] = [basename]
                else:
                    if obj_name in science.keys():
                        cnt = 2
                        while True:
                            if obj_name.rstrip() + f"_visit{cnt}" in science.keys():
                                cnt += 1
                            else:
                                science[obj_name.rstrip() + f"_visit{cnt}"] = [basename]
                                break
                    else:
                        science[obj_name] = [basename]

    if coadd_mode == "prompt":
        create_new = True
        old_science_full = None
        if os.path.isfile(mode_save_fn):
            try:
                with open(mode_save_fn, "r") as f:
                    old_science_full = json.load(f)
                print(f"\nRead coadd association file {os.path.basename(mode_save_fn)}")
                if camera in old_science_full:
                    old_science = old_science_full[camera]
                else:
                    print(f"Could not find arm '{camera}' in coadd association file.")
                    raise
                # Check for any new files in directory
                old_list = list(chain.from_iterable(old_science.values()))
                new_list = list(chain.from_iterable(science.values()))
                if sorted(old_list) == sorted(new_list):
                    print("Existing coadd associations:")
                    for i, (k, v) in enumerate(sorted(old_science.items())):
                        v.sort()
                        print(f"{i} - OBJECT: {k}  -  {v}")
                    print()
                    while True:
                        confirm = input("Use this set of associations? (Y/N): ")
                        if confirm.upper() == "Y" or confirm.upper() == "YES":
                            create_new = False
                            science = old_science
                            print()
                            break
                        elif confirm.upper() == "N" or confirm.upper() == "NO":
                            print()
                            break
                        else:
                            print("** Invalid input. Try again. **")
                else:
                    print("Input files have changed since definition of previous coadd association file.")
                    raise
            except:
                print("Will create a new file.")
                pass
        while create_new:
            obs_dict = {}
            print("\nSelect observation numbers to be coadded for an object (one object at a time):\n")
            for i, (k, v) in enumerate(sorted(science.items())):
                v.sort()
                print(f"{i} - OBJECT: {k}  -  {v}")
                obs_dict[str(i)] = [k, v]
            print()
            try:
                prompt_list = input("Enter comma- or space-delimited list of observation numbers to coadd "
                                    "(hit return when done): ")
                if prompt_list:
                    prompt_list = [int(x.strip()) for y in prompt_list.split(',') for x in y.split()]
                    if any(pl < 0 for pl in prompt_list) or any(pl > len(obs_dict) - 1 for pl in prompt_list):
                        print("\n** Values outside of range **")
                        raise ValueError
                    if len(prompt_list) != len(set(prompt_list)):
                        print("\n** No duplicate values allowed **")
                        raise ValueError
                    for obsnum in prompt_list[1:]:
                        science[obs_dict[str(prompt_list[0])][0]].extend(obs_dict[str(obsnum)][1])
                        del science[obs_dict[str(obsnum)][0]]
                else:
                    break
            except ValueError:
                print("\n** Received inappropriate input! **\nTry again or hit return to exit\n")
                continue
        if create_new:
            if old_science_full is None:
                old_science_full = {}
            old_science_full[camera] = science
            try:
                with open(mode_save_fn, "w") as f:
                    json.dump(old_science_full, f)
            except:
                raise IOError(f"Could not write coadd association file {mode_save_fn}")

    # #------------------
    # science dictionay

    sci_obs = []

    for obj_name in science.keys():
        # sort to ensure coaddds get identical names in each arm
        obs_list = sorted(science[obj_name])
        sci_obs.append({"sci": obs_list, "sky": []})

    # ------------------
    # stdstars dictionary
    std_obs = []

    for obj_name in stdstar.keys():
        # sort to ensure coaddds get identical names in each arm
        obs_list = sorted(stdstar[obj_name])
        std_obs.append(
            {"sci": obs_list, "name": obj_name, "type": ["flux", "telluric"]}
        )

    obs_metadata = {
        "bias": bias,
        "domeflat": domeflat,
        "twiflat": twiflat,
        "dark": dark,
        "wire": wire,
        "arc": arc,
        "sci": sci_obs,
        "std": std_obs,
    }

    return obs_metadata


def classify(data_dir, naxis2_to_process=0, greedy_stds=False, coadd_mode='all', mode_save_fn=None):
    """
    Classify FITS files in the specified directory based on the CAMERA keyword in the header. It filters files into blue and red (arms) observations, extracting metadata for each. It returns a dictionary containing metadata for the blue and red observations.


    Parameters
    ----------
    data_dir : str
        The directory containing FITS files to classify.
    naxis2_to_process : int, optional
        The value of the NAXIS2 keyword to filter files by. Defaults to 0 (no filtering).
    greedy_stds : bool
        Whether to treat OBJECT images near known standard stars as STANDARDs
    coadd_mode : str
        Whether to group together OBJECT frames for the same target. Options:
        'all' (default), 'none', 'prompt' (user-selected).
    mode_save_fn : str
        Filename for saving grouping of images to coadd if coadd_mode = 'prompt'.

    Returns
    -------
    dict
        A dictionary containing metadata for the blue and red observations.

    """

    # Get list of all fits files in directory and sort to ensure repeatability
    filenames = sorted(os.listdir(data_dir))

    # Filtering the data as per blue and red arm
    blue_filenames = []
    red_filenames = []

    for filename in filenames:
        try:
            f = pyfits.open(data_dir + filename)
            camera = f[0].header["CAMERA"]
            naxis2 = f[0].header["NAXIS2"]
            f.close()
        except:
            continue
        if naxis2_to_process != 0 and naxis2_to_process != naxis2:
            continue
        if camera == "WiFeSBlue":
            if filename in blue_filenames:
                continue
            else:
                blue_filenames.append(filename)
        if camera == "WiFeSRed":
            if filename in red_filenames:
                continue
            else:
                red_filenames.append(filename)

    blue_obs_metadata = get_obs_metadata(blue_filenames, data_dir, greedy_stds=greedy_stds,
                                         coadd_mode=coadd_mode, mode_save_fn=mode_save_fn, camera="blue")
    red_obs_metadata = get_obs_metadata(red_filenames, data_dir, greedy_stds=greedy_stds,
                                        coadd_mode=coadd_mode, mode_save_fn=mode_save_fn, camera="red")

    return {"blue": blue_obs_metadata, "red": red_obs_metadata}


def cube_matcher(paths_list):
    arm_list = []
    date_obs_list = []
    for path in paths_list:
        fits_header = pyfits.getheader(path)
        arm_list.append(fits_header["ARM"])
        date_obs_list.append(fits_header["DATE-OBS"])

    df = pd.DataFrame({"path": paths_list, "arm": arm_list, "date_obs": date_obs_list})

    matched_obs_arms = (
        df.groupby("date_obs")[["path", "arm"]]
        .apply(lambda x: x.to_dict("records"))
        .tolist()
    )
    matched_dicts = []

    for obs_arms in matched_obs_arms:
        matched_dict = {"Blue": None, "Red": None, "file_name": None}
        for obs_arm in obs_arms:
            matched_dict[obs_arm["arm"]] = obs_arm["path"]

            if obs_arm["path"] is not None:
                base = os.path.basename(obs_arm["path"])
                # Remove extention, Blue and Red- label to form the name
                file_name = (
                    os.path.splitext(base)[0].replace("Red--", "").replace("Blue-", "")
                )
                matched_dict["file_name"] = file_name
        matched_dicts.append(matched_dict)

    return matched_dicts
