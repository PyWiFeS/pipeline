import os
from astropy.io import fits as pyfits
import itertools
import pandas as pd
import pyjson5
import re
import string
import sys

from pywifes import wifes_calib


def _column_name_generator():
    # Create ever-increasing capital letter sequence, growing in length as required.
    # From https://stackoverflow.com/questions/42099312/how-to-create-an-infinite-iterator-to-generate-an-incrementing-alphabet-pattern

    for i in itertools.count(1):
        for p in itertools.product(string.ascii_uppercase, repeat=i):
            yield ''.join(p)


def get_obs_metadata(filenames, data_dir, greedy_stds=False, coadd_mode="all", mode_save_fn=None, camera="blue"):
    """
    Retrieve metadata for observed data files.

    This function categorizes observed data files based on their image type
    ('BIAS', 'DARK', 'FLAT', 'SKYFLAT', 'ARC', 'WIRE', 'STANDARD', 'OBJECT', 'SKY')
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
        - 'dark': List of filenames of dark frames.
        - 'domeflat': List of filenames of domeflat frames.
        - 'twiflat': List of filenames of twilight flat frames.
        - 'arc': List of filenames of arc frames.
        - 'wire': List of filenames of wire frames.
        - 'sci': List of dictionaries, each containing information about science observations. Each dictionary has the following keys:
            - 'sci': List of filenames of science observations.
            - 'sky': List of filenames of sky observations.
        - 'std': List of dictionaries, each containing information about standard star observations. Each dictionary has the following keys:
            - 'sci': List of filenames of standard star observations.
            - 'name': Name of the standard star.
            - 'type': List of strings indicating the type of observation ('flux','telluric').

    """
    stdstar_list = wifes_calib.ref_fname_lookup.keys()

    blue_grating = []
    red_grating = []
    blue_beamsplitter = []
    red_beamsplitter = []

    # classify each obs
    bias = []
    dark = []
    domeflat = []
    twiflat = []
    arc = []
    wire = []
    stdstar = {}
    science = {}
    sky = {}
    sky_assoc = {}

    for filename in filenames:
        basename = filename.replace(".fits", "")

        f = pyfits.open(data_dir + filename)
        imagetype = f[0].header["IMAGETYP"].upper()
        obj_name = f[0].header["OBJECT"]
        this_gratingb = f[0].header["GRATINGB"]
        this_gratingr = f[0].header["GRATINGR"]
        this_beamsplitter = f[0].header["BEAMSPLT"]
        f.close()
        # ---------------------------
        # Check if it is within a close distance to a standard star.
        # If so and greedy_stds is True, fix the object name to be the good one from the list
        std_type = None
        try:
            if (greedy_stds or imagetype == "STANDARD"):
                near_std, std_dist, temp_std_type = wifes_calib.find_nearest_stdstar(data_dir + filename, stdtype="any")
                if std_dist < 100.0:
                    obj_name = near_std
                    std_type = temp_std_type
        except:
            pass
        # ---------------------------
        # 1 - bias frames
        if imagetype == "BIAS" or imagetype == "ZERO":
            bias.append(basename)
        # 2 - dark frames
        elif imagetype == "DARK":
            dark.append(basename)
        else:
        # Open shutter observations: Check for mixtures of gratings or beamsplitters
            if camera == "blue":
                if this_gratingb not in blue_grating:
                    blue_grating.append(this_gratingb)
                if this_beamsplitter not in blue_beamsplitter:
                    blue_beamsplitter.append(this_beamsplitter)
            if camera == "red":
                if this_gratingr not in red_grating:
                    red_grating.append(this_gratingr)
                if this_beamsplitter not in red_beamsplitter:
                    red_beamsplitter.append(this_beamsplitter)

        # 3 - internal quartz lamp (dome) flats
        if imagetype == "FLAT":
            domeflat.append(basename)
        # 4 - twilight flats
        elif imagetype == "SKYFLAT":
            twiflat.append(basename)
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
                stdstar[obj_name]['flist'].append(basename)
            else:
                stdstar[obj_name] = {'flist': [basename], 'stdtype': std_type}
        # 8 - science targets (also allow for standard stars in imagetype = OBJECT)
        elif imagetype == "OBJECT":
            if greedy_stds and obj_name in stdstar_list:
                # group standard obs together!
                if obj_name in stdstar.keys():
                    stdstar[obj_name]['flist'].append(basename)
                else:
                    stdstar[obj_name] = {'flist': [basename], 'stdtype': std_type}
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
        # 9 - sky observations
        elif imagetype == "SKY":
            if coadd_mode == "all":
                # group sky frames together, if desired
                if obj_name in sky_assoc.keys():
                    sky_assoc[obj_name].append(basename)
                else:
                    sky_assoc[obj_name] = [basename]
            elif coadd_mode == "none":
                # make an assumption of one sky per science frame, associating 1st with 1st, 2nd with 2nd, etc.
                if obj_name in sky_assoc.keys():
                    cnt = 2
                    while True:
                        if obj_name.rstrip() + f"_visit{cnt}" in sky_assoc.keys():
                            cnt += 1
                        else:
                            sky_assoc[obj_name.rstrip() + f"_visit{cnt}"] = [basename]
                            break
                else:
                    sky_assoc[obj_name] = [basename]
            elif coadd_mode == "prompt":
                # otherwise keep all sky frames separate and invert key,val for manual association
                sky[basename] = obj_name

    if len(blue_beamsplitter) > 1:
        print(f"WARNING: MIXING BEAMSPLITTERS (BLUE) CAN YIELD BAD REDUCTIONS: {', '.join(bb for bb in blue_beamsplitter)}")
    if len(red_beamsplitter) > 1:
        print(f"WARNING: MIXING BEAMSPLITTERS (RED) CAN YIELD BAD REDUCTIONS: {', '.join(bb for bb in red_beamsplitter)}")
    if len(blue_grating) > 1:
        print(f"WARNING: MIXING GRATINGS (BLUE) CAN YIELD BAD REDUCTIONS: {', '.join(gg for gg in blue_grating)}")
    if len(red_grating) > 1:
        print(f"WARNING: MIXING GRATINGS (RED) CAN YIELD BAD REDUCTIONS: {', '.join(gg for gg in red_grating)}")


    if coadd_mode == "prompt":
        create_new = True
        old_science_full = None
        if os.path.isfile(mode_save_fn):
            try:
                with open(mode_save_fn, "r") as f:
                    old_science_full = pyjson5.load(f)
                print(f"\nRead coadd association file {os.path.basename(mode_save_fn)}")
                if camera in old_science_full:
                    old_science = old_science_full[camera]
                else:
                    print(f"Could not find arm '{camera}' in coadd association file.")
                    raise
                # Check for any new files in directory
                old_list = list(itertools.chain.from_iterable([old_science[x]['sci'] for x in old_science.keys()]))
                new_list = list(itertools.chain.from_iterable(science.values()))
                if sorted(old_list) == sorted(new_list):
                    print("Existing coadd associations:")
                    ilen = len(str(len(old_science)))
                    klen = len(max(science.keys(), key=len))
                    vlen = len(max(new_list, key=len))
                    print(f"{'N'.ljust(ilen)}   {'Object'.ljust(klen)}")
                    print(f"{'='*ilen}   {'='*klen}")
                    for i, (k, v) in enumerate(sorted(old_science.items())):
                        v['sci'].sort()
                        v['sky'].sort()
                        print(f"{str(i).ljust(ilen)} - {k.ljust(klen)}\n{''.ljust(ilen)} - Sci: {', '.join(vv for vv in v['sci']).ljust(vlen)}\n{''.ljust(ilen)} - Sky: {', '.join(vvv for vvv in v['sky'])}")
                        print()
                    print()
                    while True:
                        confirm = input("Use this set of associations? (Y/N): ")
                        if confirm.upper() == "Y" or confirm.upper() == "YES":
                            create_new = False
                            for oskey in old_science.keys():
                                science[oskey] = old_science[oskey]['sci']
                                sky_assoc[oskey] = old_science[oskey]['sky']
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
            ilen = len(str(len(science)))
            klen = len(max(science.keys(), key=len))
            print("\nSelect observation numbers to be coadded for an object (one object at a time):\n")
            print(f"{'N'.ljust(ilen)}   {'Object'.ljust(klen)}   Files")
            print(f"{'='*ilen}   {'='*klen}   =====")
            for i, (k, v) in enumerate(sorted(science.items())):
                v.sort()
                print(f"{str(i).ljust(ilen)} - {k.ljust(klen)} - {', '.join(vv for vv in v)}")
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
                    print()
                    break
            except ValueError:
                print("\n** Received inappropriate input! **\nTry again or hit return to exit\n")
                continue

        # Now check if there are sky exposures to associate with science targets
        while create_new:
            if len(sky) == 0:
                break
            obs_dict = {}
            sky_dict = {}
            ilen = len(str(len(science)))
            klen = len(max(science.keys(), key=len))
            vlen = len(max(list(itertools.chain.from_iterable(science.values())), key=len))
            print("\n\nSpecify one object number and any space- or comma-delimited sequence of sky exposure letters:\n(e.g. '1 B C' to associate B and C to 1, or '0' to associate none to 0)\n")
            print(f"{'N'.ljust(ilen)}   {'Object'.ljust(klen)}")
            print(f"{'='*ilen}   {'='*klen}")
            for i, (k, v) in enumerate(sorted(science.items())):
                v.sort()
                sky_string = ', '.join(ss for ss in sky_assoc[k]) if k in sky_assoc else ''
                print(f"{str(i).ljust(ilen)} - {k.ljust(klen)}\n{''.ljust(ilen)} - Sci: {', '.join(vv for vv in v).ljust(vlen)}\n{''.ljust(ilen)} - Sky: {sky_string}")
                obs_dict[str(i)] = [k, v]
            ilen = max(len(str(len(sky) // 26)) + 1, 2)
            klen = len(max(sky.keys(), key=len))
            vlen = len(max(list(itertools.chain.from_iterable(sky.values())), key=len))
            print(f"\n------ Sky Files ------\nID   {'File'.ljust(klen)}   Header label")
            print(f"==   {'='*klen}   ============")
            letset = []
            for letter, (k, v) in zip(_column_name_generator(), sorted(sky.items())):
                print(f"{letter.ljust(ilen)} - {k.ljust(klen)} - {v}")
                letset.append(letter)
                sky_dict[letter] = [k, v]
            print()
            try:
                prompt_list = input("Enter one number (science frame[s]) and a comma- or space-delimited list of "
                                    "letters (sky frame[s]) to associate (hit return when done): ")
                if prompt_list:
                    prompt_list = [x.strip().upper() for y in prompt_list.split(',') for x in y.split()]
                    prompt_list[0] = int(prompt_list[0])
                    if prompt_list[0] < 0 or prompt_list[0] > len(obs_dict) - 1:
                        print(f"\n** Science frame value {prompt_list[0]} outside of range **")
                        raise ValueError
                    if any(pl not in letset for pl in prompt_list[1:]):
                        print(f"\n** Bad selection of sky frame specifier: not all of {prompt_list[1:]} in {letset} **")
                        raise ValueError
                    if len(prompt_list) != len(set(prompt_list)):
                        print("\n** No duplicate values allowed **")
                        raise ValueError
                    sky_assoc[obs_dict[str(prompt_list[0])][0]] = []
                    for skylet in prompt_list[1:]:
                        sky_assoc[obs_dict[str(prompt_list[0])][0]].append(sky_dict[skylet][0])
                else:
                    print()
                    break
            except ValueError:
                print("\n** Received inappropriate input! **\nTry again or hit return to exit\n")
                continue
        # Save the output file
        if create_new:
            if old_science_full is None:
                old_science_full = {}
            if camera not in old_science_full:
                old_science_full[camera] = {}
            for ss in science.keys():
                if ss in sky_assoc.keys():
                    old_science_full[camera][ss] = {'sci': science[ss], 'sky': sky_assoc[ss]}
                else:
                    old_science_full[camera][ss] = {'sci': science[ss], 'sky': []}
            try:
                with open(mode_save_fn, "wb") as f:
                    pyjson5.dump(old_science_full, f)
            except:
                raise IOError(f"Could not write coadd association file {mode_save_fn}")

    # #------------------
    # science dictionay

    sci_obs = []

    for obj_name in science.keys():
        # sort to ensure coaddds get identical names in each arm
        obs_list = sorted(science[obj_name])
        if obj_name in sky_assoc.keys():
            sky_list = sorted(sky_assoc[obj_name])
        else:
            sky_list = []
        sci_obs.append({"sci": obs_list, "sky": sky_list})

    # ------------------
    # stdstars dictionary
    std_obs = []

    for obj_name in stdstar.keys():
        # sort to ensure coaddds get identical names in each arm
        obs_list = sorted(stdstar[obj_name]['flist'])
        std_obs.append(
            {"sci": obs_list, "name": obj_name, "stdtype": stdstar[obj_name]['stdtype']}
        )

    obs_metadata = {
        "bias": sorted(bias),
        "dark": sorted(dark),
        "domeflat": sorted(domeflat),
        "twiflat": sorted(twiflat),
        "wire": sorted(wire),
        "arc": sorted(arc),
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
                                         coadd_mode=coadd_mode, mode_save_fn=mode_save_fn, 
                                         camera="blue")
    red_obs_metadata = get_obs_metadata(red_filenames, data_dir, greedy_stds=greedy_stds,
                                        coadd_mode=coadd_mode, mode_save_fn=mode_save_fn, 
                                        camera="red")

    return {"blue": blue_obs_metadata, "red": red_obs_metadata}


def cube_matcher(paths_list):
    arm_list = []
    date_obs_list = []
    for path in paths_list:
        fits_header = pyfits.getheader(path)
        if "ARM" in fits_header:
            arm_list.append(fits_header["ARM"])
        elif "CAMERA" in fits_header:
            arm_list.append(re.sub("WiFeS", "", fits_header["CAMERA"]))
        else:
            arm_list.append("unknown_arm")
        date_obs_list.append(fits_header["DATE-OBS"])

    if not paths_list:
        print("Found no cubes to process. Exiting.")
        return None

    df = pd.DataFrame({"path": paths_list, "arm": arm_list, "date_obs": date_obs_list})
    df["date_obs"] = pd.to_datetime(df["date_obs"], utc=True, format="ISO8601")

    df = df.sort_values("date_obs")
    Nsec = 15  # allow Nsec seconds of difference between the arms
    df['Group'] = df["date_obs"].diff().dt.seconds.gt(Nsec).cumsum()
    matched_obs_arms = (
        df.groupby('Group')[["path", "arm"]]
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
