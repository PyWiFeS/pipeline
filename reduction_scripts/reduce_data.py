#! /usr/bin/env python3

# ------------------------------------------------------------------------
# Initial set ups and import modules
# ------------------------------------------------------------------------
import os
import pickle
from astropy.io import fits as pyfits
import gc
import datetime
import numpy as np
import json
import shutil
import glob
import argparse
import logging
from pywifes.logger_config import setup_logger, custom_print
from pywifes.data_classifier import classify, cube_matcher
from pywifes.extract_spec import detect_extract_and_save, plot_1D_spectrum
from pywifes.splice import splice_spectra, splice_cubes
from pywifes.lacosmic import lacos_wifes
from pywifes import pywifes
from pywifes import wifes_wsol
from pywifes import wifes_calib
from pywifes.pywifes import is_halfframe
from pywifes.pywifes import calib_to_half_frame
import shutil
import glob
import argparse

# Set paths
reduction_scripts_dir = os.path.dirname(__file__)
working_dir = os.getcwd()

# Setup the logger.
log_file = os.path.join(working_dir, "data_products/pywifes_logger.log")
logger = setup_logger(file=log_file, console_level=logging.WARNING, file_level=logging.INFO)

# Redirect print statements to logger with different levels
debug_print = custom_print(logger, logging.DEBUG)
info_print = custom_print(logger, logging.INFO)
warning_print = custom_print(logger, logging.WARNING)
error_print = custom_print(logger, logging.ERROR)
critical_print = custom_print(logger, logging.CRITICAL)

# Redirect warnings to logger
#logging.captureWarnings(True)

info_print("Starting PyWiFeS data reduction pipeline.")


# ------------------------------------------------------------------------
# Function definition
# ------------------------------------------------------------------------

def move_files(src_dir_path, destination_dir_path, filenames):
    """
    Move files from the source directory to the destination directory.

    Parameters
    ----------
    src_dir_path : str
        The path of the source directory.
    destination_dir_path : str
        The path of the destination directory.
    filenames : list of str
        The list of filenames to be moved.

    Raises
    ------
    Exception
        If an error occurs while moving the files.

    """
    try:
        for file in filenames:
            src_file = os.path.join(src_dir_path, file)
            dest_file = os.path.join(destination_dir_path, file)
            debug_print(f"Moving file {src_file} to {dest_file}")
            shutil.move(src_file, dest_file)
    except Exception as e:
        error_print(f"Error moving files: {e}")

def copy_files(src_dir_path, destination_dir_path, filenames):
    """
    Copy files from the source directory to the destination directory.

    Parameters
    ----------
    src_dir_path : str
        The path to the source directory.
    destination_dir_path : str
        The path to the destination directory.
    filenames : list
        A list of filenames to be copied.

    Raises
    ------
    Exception
        If there is an error while copying the files.

    """
    try:
        for file in filenames:
            src_file = os.path.join(src_dir_path, file)
            dest_file = os.path.join(destination_dir_path, file)
            #debug_print(f"Copying file {src_file} to {dest_file}")
            shutil.copy(src_file, dest_file)
    except Exception as e:
        error_print(f"Error copying files: {e}")

def get_file_names(src_dir_path, glob_pattern):
    """
    Get the names of files in a directory that match a given search query in a glob pattern.

    Parameters
    ----------
    src_dir_path : str
        The path to the source directory.

    glob_pattern : str
        The glob pattern (search query) to match the filenames.

    Returns
    -------
    list
        A list of filenames that match the glob pattern.

    """
    filepaths = glob.glob(os.path.join(src_dir_path, glob_pattern))
    names = []
    for filepath in filepaths:
        filename = os.path.basename(filepath)
        names.append(filename)
    return names

def load_config_file(filename):
    """
    Load a configuration file in JSON format.

    Parameters
    ----------
    filename : str
        The name of the configuration file.

    Returns
    -------
    dict
        A dictionary containing the loaded configuration data.

    """
    reduction_scripts_dir = os.path.dirname(__file__)
    file_path = os.path.join(reduction_scripts_dir, filename)
    info_print(f"Loading configuration file: {file_path}")
    with open(file_path, "r") as f:
        return json.load(f)

def get_full_obs_list(metadata):
    """
    Get the full observation list from the given metadata.

    Parameters:
    ----------
        metadata : dict
            A dictionary containing metadata information.

    Returns:
    -------
        list
            The full observation list.

    """
    full_obs_list = []
    base_fn_list = (
        metadata["bias"]
        + metadata["arc"]
        + metadata["wire"]
        + metadata['dark']
        + metadata["domeflat"]
        + metadata["twiflat"]
    )
    for fn in base_fn_list:
        if fn not in full_obs_list:
            full_obs_list.append(fn)
    for obs in metadata["sci"] + metadata["std"]:
        for key in obs.keys():
            if key != "type" and key != "name":
                for fn in obs[key]:
                    if fn not in full_obs_list:
                        full_obs_list.append(fn)
    debug_print(f"Full observation list: {full_obs_list}")
    return full_obs_list



def get_sci_obs_list(metadata):
    """
    Get a list of science observations from the metadata.

    Parameters
    ----------
    metadata : dict
        The metadata containing information about the observations.

    Returns
    -------
    list
        A list of science observation filenames.

    """
    sci_obs_list = []
    for obs in metadata["sci"]:
        for fn in obs["sci"]:
            if fn not in sci_obs_list:
                sci_obs_list.append(fn)
    debug_print(f"Science observation list: {sci_obs_list}")
    return sci_obs_list

def get_std_obs_list(metadata, type="all"):
    """
    Get a list of standard observations.

    Parameters:
    ----------
    - metadata : dict
        The metadata containing information about the observations.
    - type : str, optional
        The type of observations to include in the list. Default is "all".

    Returns:
    --------
    - std_obs_list : list
        A list of standard observation filenames.

    """
    std_obs_list = []
    for obs in metadata["std"]:
        for fn in obs["sci"]:
            if fn not in std_obs_list and type == "all":
                std_obs_list.append(fn)
            if fn not in std_obs_list and (type in obs["type"]):
                std_obs_list.append(fn)
    debug_print(f"Standard observation list ({type}): {std_obs_list}")
    return std_obs_list

def get_sky_obs_list(metadata):
    """
    Get a list of sky observations from the metadata.

    Parameters:
    ----------
        metadata : dict
            The metadata containing information about the observations.

    Returns:
    --------
        list
            A list of sky observation filenames.

    """
    sky_obs_list = []
    for obs in metadata["sci"] + metadata["std"]:
        if "sky" not in obs.keys():
            continue
        for fn in obs["sky"]:
            if fn not in sky_obs_list:
                sky_obs_list.append(fn)
    info_print(f"Sky observation list: {sky_obs_list}")            
    return sky_obs_list

def get_associated_calib(metadata, this_fn, type):
    """
    Get the associated calibration file for a given data file.

    Parameters
    ----------
    metadata : dict
        The metadata dictionary containing information about the observations.
    this_fn : str
        The filename of the data file for which to find the associated calibration file.
    type : str
        The type of calibration file to search for.

    Returns
    -------
    str or bool
        The filename of the associated calibration file if found, or False if not found.
    """
    for obs in metadata["sci"] + metadata["std"]:
        if "sky" in obs.keys():
            sky = obs["sky"]
        else:
            sky = []
        for fn in obs["sci"] + sky:
            if fn == this_fn:
                if type in obs.keys():
                    if obs[type] != "":
                        return obs[type]
    return False

def get_primary_sci_obs_list(metadata):
    """
    Get the list of primary science observations from the metadata.

    Parameters
    ----------
    metadata : dict
        The metadata containing information about the observations.

    Returns
    -------
    list
        The list of primary science observations.

    """
    sci_obs_list = [obs["sci"][0] for obs in metadata["sci"]]
    debug_print(f"Primary science observation list: {sci_obs_list}")
    return sci_obs_list

def get_primary_std_obs_list(metadata, type="all"):
    """
    Get the list of primary standard observations based on the given metadata and type.

    Parameters:
    ----------
        metadata : dict
            The metadata containing information about the observations.
        type : str, optional
            The type of standard star observations to include in the list.
            Possible values are "all", "telluric", or "flux". Defaults to "all".

    Returns:
    --------
        list
            The list of primary standard observations.

    Raises:
    -------
        ValueError 
            If the standard star type is not understood.

    """
    if type == "all":
        std_obs_list = [obs["sci"][0] for obs in metadata["std"]]
    elif type == "telluric" or type == "flux":
        std_obs_list = []
        for obs in metadata["std"]:
            if obs["sci"][0] not in std_obs_list and (type in obs["type"]):
                std_obs_list.append(obs["sci"][0])
    else:
        error_print("Standard star type not understood!")
        error_print("PyWiFeS Data Reduction pipeline will crash now ...")
        raise ValueError("Standard star type not understood")
    debug_print(f"Primary standard observation list ({type}): {std_obs_list}")
    return std_obs_list

# ------------------------------------------------------------------------


def main():
    start_time = datetime.datetime.now()
    info_print(f"Pipeline started at {start_time}")

    # ------------------------------------------------------------------------
    # DEFINE THE PROCESSING STEPS
    # ------------------------------------------------------------------------

    # ------------------------------------------------------------------------
    # Overscan subtraction
    # ------------------------------------------------------------------------
    def run_overscan_sub(metadata, prev_suffix, curr_suffix):
        """
        Subtract overscan from the input FITS files and save the results.

        Parameters
        ----------
        metadata : dict
            A dictionary containing metadata information of the FITS files.
        prev_suffix : str
            The suffix of the previous FITS files.
        curr_suffix : str
            The suffix of the current FITS files.

        Returns
        -------
        None
        """
        full_obs_list = get_full_obs_list(metadata)
        for fn in full_obs_list:
            in_fn = os.path.join(temp_data_dir, "%s.fits" % fn)
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            if skip_done and os.path.isfile(out_fn):
                continue
            info_print(f"Subtracting Overscan for {in_fn.split('/')[-1]}")
            pywifes.subtract_overscan(in_fn, out_fn, data_hdu=0)
        return

    # ------------------------------------------------------
    # repair bad pixels!
    # ------------------------------------------------------
    def run_bpm_repair(metadata, prev_suffix, curr_suffix):
        """
        Repairs bad pixels in the input FITS files and saves the repaired files.

        Parameters
        ----------
        metadata : dict
            A dictionary containing metadata information of the FITS files.
        prev_suffix : str
            The suffix of the previous version of the FITS files.
        curr_suffix : str
            The suffix of the current version of the FITS files.

        Returns
        -------
        None
        """
        full_obs_list = get_full_obs_list(metadata)
        for basename in full_obs_list:
            input_filename = f"{basename}.p{prev_suffix}.fits"
            output_filename = f"{basename}.p{curr_suffix}.fits"
            input_filepath = os.path.join(out_dir, input_filename)
            output_filepath = os.path.join(out_dir, output_filename)
            if skip_done and os.path.isfile(output_filepath):
                continue
            info_print(f"Repairing {arm} bad pixels for {input_filename}")
            if arm == "red":
                pywifes.repair_red_bad_pix(
                    input_filepath, output_filepath, data_hdu=0
                )
            if arm == "blue":
                pywifes.repair_blue_bad_pix(
                    input_filepath, output_filepath, data_hdu=0
                )

    # ------------------------------------------------------
    # Generate super-bias
    # ------------------------------------------------------
    def run_superbias(metadata, prev_suffix, curr_suffix, method="row_med", **args):
        '''
        Generate superbias for the entire dataset and for each science frame.
        
        Parameters:
        - metadata: dict
            Metadata containing information about the dataset.
        - prev_suffix: str
            Previous suffix used in the filenames.
        - curr_suffix: str
            Current suffix to be used in the filenames.
        - method: str, optional
            Method to generate the superbias. 
            Fit a smart surface to the bias or take the median of each row.
            Default is "row_med".
        - **args: dict
            Additional arguments to be passed to the generate_wifes_bias_fit function.
        
        Returns:
        None
        '''
        bias_list = [
            os.path.join(out_dir, "%s.p%s.fits" % (x, prev_suffix))
            for x in metadata["bias"]
        ]
        info_print("Calculating Global Superbias")
        pywifes.imcombine(bias_list, superbias_fn, data_hdu=0)
        if method == "fit" or method == "row_med":
            pywifes.generate_wifes_bias_fit(
                superbias_fn,
                superbias_fit_fn,
                data_hdu=0,
                method=method,
                plot_dir=plot_dir_arm,
                arm=arm, 
                **args,
            )
        else:
            pywifes.imcopy(superbias_fn, superbias_fit_fn)
        # generate local superbiases for any science frames
        sci_obs_list = get_sci_obs_list(metadata)
        std_obs_list = get_std_obs_list(metadata)
        info_print('***************************************************')
        info_print(f"Science observation list: {sci_obs_list}")
        info_print(f"Standard observation list: {std_obs_list}")
        info_print('***************************************************')
        for fn in sci_obs_list + std_obs_list:
            local_biases = get_associated_calib(metadata, fn, "bias")
            if local_biases:
                local_bias_fn = get_associated_calib(metadata, fn, "bias")[0]
                info_print(f"Calculating Local Superbias for {local_bias_fn}")
                local_superbias = os.path.join(
                    out_dir, "%s.fits" % (local_bias_fn + ".lsb")
                )
                local_superbias_fit = os.path.join(
                    out_dir, "%s.fits" % (local_bias_fn + ".lsb_fit")
                )
                if os.path.isfile(local_superbias_fit):
                    continue
                # step 1 - coadd biases
                local_biases_filename = [
                    os.path.join(out_dir, "%s.p%s.fits" % (x, prev_suffix))
                    for x in local_biases
                ]
                pywifes.imcombine(
                    local_biases_filename, local_superbias, data_hdu=0
                )
                # step 2 - generate fit
                if method == "fit" or method == "row_med":
                    pywifes.generate_wifes_bias_fit(
                        local_superbias,
                        local_superbias_fit,
                        data_hdu=0,
                        method=method,
                        **args,
                    )
                else:
                    pywifes.imcopy(local_superbias, local_superbias_fit)
        return

    # ----------------------------------------------------
    # Subtract bias
    # ----------------------------------------------------
    def run_bias_sub(metadata, prev_suffix, curr_suffix, method="sub", **args):
        """
        Subtract bias from all the input data FITS files.

        Parameters
        ----------
        metadata : dict
            Metadata dictionary containing information about the FITS files of the observations.
        prev_suffix : str
            Previous suffix of the data files.
        curr_suffix : str
            Current suffix of the data files.
        method : str, optional
            Method for bias subtraction. Default is subtraction "sub".
            The other option is "copy" for copying the files without subtraction.
        **args : dict
            Additional keyword arguments.

        Returns
        -------
        None
        """
        full_obs_list = get_full_obs_list(metadata)
        for fn in full_obs_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            if skip_done and os.path.isfile(out_fn):
                continue
            local_biases = get_associated_calib(metadata, fn, "bias")
            if local_biases:
                local_bias_fn = get_associated_calib(metadata, fn, "bias")[0]
                local_superbias = os.path.join(
                    out_dir, "%s.fits" % (local_bias_fn + ".lsb")
                )
                bias_fit_fn = os.path.join(
                    out_dir, "%s.fits" % (local_bias_fn + ".lsb_fit")
                )
                bias_type = "local"
            else:
                bias_fit_fn = superbias_fit_fn
                bias_type = "global"

            # subtract it!
            info_print(f"Subtracting {bias_type} superbias for {os.path.basename(in_fn)}")
            if method == "copy":
                pywifes.imcopy(in_fn, out_fn)
            else:
                pywifes.imarith(in_fn, "-", bias_fit_fn, out_fn, data_hdu=0)
        return

    # ------------------------------------------------------
    # Generate super-flat
    # ------------------------------------------------------
    def run_superflat(
        metadata, prev_suffix, curr_suffix, source, scale=None, method="median"
    ):
        """
        Generate a co-add superflat for a given type, either "dome" or "twi".

        Parameters
        ----------
        metadata : dict
            Metadata containing information about the FITS files from which extract the flats.
        prev_suffix : str
            Previous suffix of the input files.
        curr_suffix : str
            Current suffix of the output files.
        source : str
            Type of flatfield frame ("dome" or "twi").
        scale : float, optional
            Scale factor for the flatfield generation.
            Default is None.
        method : str, optional
            Method for combining the flatfield images.
            Default is "median".

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If the flatfield type is not recognized.

        """
        if source == "dome":
            flat_list = [
                os.path.join(out_dir, "%s.p%s.fits" % (x, prev_suffix))
                for x in metadata["domeflat"]
            ]
            info_print(f"List of {source} flats: {flat_list}")
            out_fn = super_dflat_raw
        elif source == "twi":
            flat_list = [
                os.path.join(out_dir, "%s.p%s.fits" % (x, prev_suffix))
                for x in metadata["twiflat"]
            ]
            info_print(f"List of {source} flats: {flat_list}")
            out_fn = super_tflat_raw
        else:
            error_print("Flatfield type not recognized")
            raise ValueError("Flatfield type not recognized")
        if not flat_list:
            warning_print(f"No {source} flats found. Skipping the superflat generation for {source}.")
            return
        info_print(f"Generating co-add {source} flat")
        pywifes.imcombine(
            flat_list, out_fn, data_hdu=0, scale=scale, method=method
        )
        return

    # ------------------------------------------------------
    # Flat cleanup
    # ------------------------------------------------------
    def run_flat_cleanup(
        metadata,
        prev_suffix,
        curr_suffix,
        type=["dome", "twi"],
        offsets=[0.0, 0.0],
        **args,
    ):
        '''
        Make the master domeflat and twilight flat corrections.

        Parameters
        ----------
        metadata : str
            The metadata information for the flat cleanup.
        prev_suffix : str
            The previous suffix.
        curr_suffix : str
            The current suffix.
        type : list, optional
            The types of flats to correct (default is ["dome", "twi"]).
        offsets : list, optional
            The offsets for each type of flat (default is [0.0, 0.0]).
        **args : dict, optional
            Additional keyword arguments.

        Returns
        -------
        None
        
        '''
        # check the slitlet definition file
        if os.path.isfile(slitlet_def_fn):
            slitlet_fn = slitlet_def_fn
        else:
            slitlet_fn = None
        if "dome" in type:
            if os.path.isfile(super_dflat_raw):
                info_print(f"Correcting master domeflat {os.path.basename(super_dflat_fn)}")
                pywifes.interslice_cleanup(
                    super_dflat_raw,
                    super_dflat_fn,
                    slitlet_fn,
                    offset=offsets[type.index("dome")],
                    method="2D",
                    plot_dir=plot_dir_arm,
                    **args,
                )
            else:
                warning_print(f"Master dome flat {os.path.basename(super_dflat_raw)} not found. Skipping dome flat cleanup.")

        if "twi" in type:
            if os.path.isfile(super_tflat_raw):
                info_print(f"Correcting master twilight flat {os.path.basename(super_tflat_fn)}")
                pywifes.interslice_cleanup(
                    super_tflat_raw,
                    super_tflat_fn,
                    slitlet_fn,
                    offset=offsets[type.index("twi")],
                    method="2D",
                    plot_dir=plot_dir_arm,
                    **args,
                )
            else:
                warning_print(f"Master twilight flat {os.path.basename(super_tflat_raw)} not found. Skipping twilight flat cleanup.")
        return

    # ------------------------------------------------------
    # Fit slitlet profiles
    # ------------------------------------------------------
    def run_slitlet_profile(metadata, prev_suffix, curr_suffix, **args):
        ''' 
        Fit the slitlet profiles to the flatfield.

        Parameters
        ----------
        metadata : dict
            Metadata information.
        prev_suffix : str
            Previous suffix.
        curr_suffix : str
            Current suffix.
        **args : dict
            Additional keyword arguments.

        Returns
        -------
        None
            This function does not return anything.

        Notes
        -----
        If the super_dflat_fn file exists, it uses it as the flatfield_fn. Otherwise, it uses super_dflat_raw as the flatfield_fn.

        The slitlet profiles are derived using the flatfield_fn and saved to the output_fn.

        Examples
        --------
        >>> run_slitlet_profile(metadata, prev_suffix, curr_suffix, arg1=value1, arg2=value2)
        '''
        if os.path.isfile(super_dflat_fn):
            flatfield_fn = super_dflat_fn
        else:
            flatfield_fn = super_dflat_raw
        output_fn = slitlet_def_fn
        pywifes.derive_slitlet_profiles(
            flatfield_fn, output_fn, data_hdu=0, **args
        )
        return

    # ------------------------------------------------------
    # Create MEF files
    # ------------------------------------------------------
    def run_superflat_mef(metadata, prev_suffix, curr_suffix, source):
        """
        Generate a Multi-Extension FITS (MEF) file for dome or twilight flat.

        Parameters
        ----------
        metadata : str
            The metadata information.
        prev_suffix : str
            The previous suffix.
        curr_suffix : str
            The current suffix.
        source : str
            The source of the flatfield (dome or twi).

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
        This function generates a Multi-Extension FITS (MEF) file for dome or twilight flat.
        It checks for the existence of the master dome or twilight flat file and generates the MEF file accordingly.
        If the master flat file is not found, it prints a warning message and skips the MEF generation.
        The slitlet definition file is also checked and used if available.
        Finally, the MEF file is generated using the `pywifes.wifes_slitlet_mef` function.

        """
        if source == "dome":
            if os.path.isfile(super_dflat_fn):
                in_fn = super_dflat_fn
            elif os.path.isfile(super_dflat_raw):
                in_fn = super_dflat_raw
            else:
                warning_print(f"No master dome flat found. Skipping MEF generation for dome flat.")
                return
            out_fn = super_dflat_mef

        elif source == "twi":
            if os.path.isfile(super_tflat_fn):
                in_fn = super_tflat_fn
            elif os.path.isfile(super_tflat_raw):
                in_fn = super_tflat_raw
            else:
                warning_print(f"No master twilight flat found. Skipping MEF generation for twilight flat.")
                return
            out_fn = super_tflat_mef
        else:
            error_print("Flatfield type not recognized")
            raise ValueError("Flatfield type not recognized")
        # check the slitlet definition file
        if os.path.isfile(slitlet_def_fn):
            slitlet_fn = slitlet_def_fn
        else:
            slitlet_fn = None
        # run it!
        info_print(f"Generating MEF {source} flat")
        pywifes.wifes_slitlet_mef(
            in_fn, out_fn, data_hdu=0, slitlet_def_file=slitlet_fn
        )
        return

    def run_slitlet_mef(metadata, prev_suffix, curr_suffix, ns=False):
        """
        Create a Multi-Extension Fits (MEF) file for each observation in the metadata.

        Parameters
        ----------
        metadata : dict
            Metadata containing information about the observations.
        prev_suffix : str
            Previous suffix of the fits file.
        curr_suffix : str
            Current suffix of the fits file.
        ns : bool, optional
            Flag indicating whether to process non-standard observations, by default False.

        Returns
        -------
        None
        """
        full_obs_list = get_full_obs_list(metadata)
        sci_obs_list = get_sci_obs_list(metadata)
        std_obs_list = get_std_obs_list(metadata)
        sky_obs_list = get_sky_obs_list(metadata)
        ns_proc_list = sci_obs_list + std_obs_list

        # check the slitlet definition file
        if os.path.isfile(slitlet_def_fn):
            slitlet_fn = slitlet_def_fn
        else:
            slitlet_fn = None

        for fn in full_obs_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))

            if skip_done and os.path.isfile(out_fn):
                continue

            info_print(f"Creating MEF file for {os.path.basename(in_fn)}")

            if ns and fn in ns_proc_list:
                sky_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, curr_suffix))
                pywifes.wifes_slitlet_mef_ns(
                    in_fn,
                    out_fn,
                    sky_fn,
                    data_hdu=0,
                    slitlet_def_file=slitlet_fn,
                )
            else:
                pywifes.wifes_slitlet_mef(
                    in_fn, out_fn, data_hdu=0, slitlet_def_file=slitlet_fn
                )

            gc.collect()

        return

    # ------------------------------------------------------
    # Wavelength solution
    # ------------------------------------------------------
    def run_wave_soln(metadata, prev_suffix, curr_suffix, **args):
        '''
        Wavelength Solution:
        Generate the master arc solution, based on generic arcs at first.
        Then looks for the local wavelength solutions for science or standards (sky not required at this stage).
        Check if the file has a dedicated arc associated with it.
        If two arcs are present, find a solution for both to later interpolate between them.
        Restrict it to the first two arcs in the list (in case the feature is
        being unknowingly used).

        Parameters
        ----------
        metadata : dict
            Metadata information for the data.
        prev_suffix : str
            Previous suffix of the file.
        curr_suffix : str
            Current suffix of the file.
        **args : dict
            Additional keyword arguments.

        Returns
        -------
        None
            This function does not return anything.

        '''
        wsol_in_fn = os.path.join(
            out_dir, "%s.p%s.fits" % (metadata["arc"][0], prev_suffix)
        )
        info_print(f"Deriving master wavelength solution from {os.path.basename(wsol_in_fn)}")
        wifes_wsol.derive_wifes_wave_solution(wsol_in_fn, wsol_out_fn, plot_dir=plot_dir_arm, **args)
        sci_obs_list = get_sci_obs_list(metadata)
        std_obs_list = get_std_obs_list(metadata)
       
        for fn in sci_obs_list + std_obs_list:
            local_arcs = get_associated_calib(metadata, fn, "arc")

            if local_arcs:
                for i in range(np.min([2, np.size(local_arcs)])):
                    local_arc_fn = os.path.join(
                        out_dir, "%s.p%s.fits" % (local_arcs[i], prev_suffix)
                    )
                    
                    local_wsol_out_fn = os.path.join(
                        out_dir, "%s.wsol.fits" % (local_arcs[i])
                    )
                    
                    if os.path.isfile(local_wsol_out_fn):
                        continue
                    info_print(f"Deriving local wavelength solution for {local_arcs[i]}")
                    
                    wifes_wsol.derive_wifes_wave_solution(
                        local_arc_fn, 
                        local_wsol_out_fn,
                        plot_dir=plot_dir_arm 
                        **args
                    )
                    
        return

    # ------------------------------------------------------
    # Wire solution
    # ------------------------------------------------------
    def run_wire_soln(metadata, prev_suffix, curr_suffix):
        ''' 
        Derive global wire solution first, then local wire solutions for any specific observations.

        Parameters
        ----------
        metadata : dict
            Metadata containing information about the observations.
        prev_suffix : str
            Previous suffix used in the file names.
        curr_suffix : str
            Current suffix to be used in the file names.

        Returns
        -------
        None
        '''
        # Global wire solution
        wire_in_fn = os.path.join(
            out_dir, "%s.p%s.fits" % (metadata["wire"][0], prev_suffix)
        )
        info_print(f"Deriving global wire solution from {os.path.basename(wire_in_fn)}")
        pywifes.derive_wifes_wire_solution(wire_in_fn, wire_out_fn)
        
        # Wire solutions for any specific observations
        sci_obs_list = get_sci_obs_list(metadata)
        std_obs_list = get_std_obs_list(metadata)
        for fn in sci_obs_list + std_obs_list:
            # Check if the file has a dedicated wire associated with it ...
            # Only for Science and Std stars (sky not required at this stage)
            local_wires = get_associated_calib(metadata, fn, "wire")
            if local_wires:
                local_wire_fn = os.path.join(
                    out_dir, "%s.p%s.fits" % (local_wires[0], prev_suffix)
                )
                local_wire_out_fn = os.path.join(
                    out_dir, "%s.wire.fits" % (out_dir, local_wires[0])
                )
                if os.path.isfile(local_wire_out_fn):
                    continue
                info_print(f"Deriving local wire solution for {local_wires[0]}")
                pywifes.derive_wifes_wire_solution(local_wire_fn, local_wire_out_fn)
        return

    # ------------------------------------------------------
    # Cosmic Rays
    # ------------------------------------------------------
    def run_cosmic_rays(
        metadata,
        prev_suffix,
        curr_suffix,
        ns=False,
        multithread=False,
        max_processes=-1,
    ):
        ''' 
        Clean cosmic rays on all science and standard frames.

        Parameters
        ----------
        metadata : dict
            Metadata containing information about the observations.
        prev_suffix : str
            Previous suffix of the FITS file to apply the correction (input).
        curr_suffix : str
            Current suffix of the FITS file with the correction alredy applied (output).
        ns : bool, optional
            Flag indicating whether to clean cosmics in the non-science frames (default is False).
        multithread : bool, optional
            Flag indicating whether to use multithreading for cosmic ray cleaning (default is False).
        max_processes : int, optional
            Maximum number of processes to use for multithreading (default is -1, which uses all available processes).

        Returns
        -------
        None
        '''
        # now run ONLY ON SCIENCE TARGETS AND STANDARDS
        sci_obs_list = get_sci_obs_list(metadata)
        sky_obs_list = get_sky_obs_list(metadata)
        std_obs_list = get_std_obs_list(metadata)
        for fn in sci_obs_list + sky_obs_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            info_print(f"Cleaning cosmics in {os.path.basename(in_fn)}")
            if skip_done and os.path.isfile(out_fn):
                continue
            lacos_wifes(
                in_fn,
                out_fn,
                wsol_fn=wsol_out_fn,
                niter=3,
                sig_clip=10.0,
                obj_lim=10.0,
                sig_frac=0.2,
                is_multithread=multithread,
                max_processes=max_processes,
            )
            if ns:
                in_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, prev_suffix))
                out_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, curr_suffix))
                info_print(f"Cleaning cosmics in {os.path.basename(in_fn)}")
                lacos_wifes(
                    in_fn,
                    out_fn,
                    wsol_fn=wsol_out_fn,
                    niter=3,
                    sig_clip=10.0,
                    obj_lim=10.0,
                    sig_frac=0.2,
                    is_multithread=multithread,
                    max_processes=max_processes,
                )
            gc.collect()
        for fn in std_obs_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            if skip_done and os.path.isfile(out_fn):
                continue
            info_print(f"Cleaning cosmics in standard star {os.path.basename(in_fn)}")
            lacos_wifes(
                in_fn,
                out_fn,
                wsol_fn=wsol_out_fn,
                niter=3,
                sig_clip=10.0,
                obj_lim=10.0,
                sig_frac=0.2,
                is_multithread=multithread,
                max_processes=max_processes,
            )
            if ns:
                in_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, prev_suffix))
                out_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, curr_suffix))
                info_print(f"Cleaning cosmics in standard star {os.path.basename(in_fn)}")
                lacos_wifes(
                    in_fn,
                    out_fn,
                    wsol_fn=wsol_out_fn,
                    niter=3,
                    sig_clip=10.0,
                    obj_lim=10.0,
                    sig_frac=0.2,
                    is_multithread=multithread,
                    max_processes=max_processes,
                )
            gc.collect()
        return

    # ------------------------------------------------------
    # Sky subtraction
    # ------------------------------------------------------
    def run_sky_sub_ns(metadata, prev_suffix, curr_suffix):
        '''
        Subtract sky frames from science objects in nod-and-shuffle mode.

        Parameters
        ----------
        metadata : dict
            Metadata containing information about the observations.
        prev_suffix : str
            Previous suffix of the file names (input).
        curr_suffix : str
            Current suffix of the file names (output).

        Returns
        -------
        None
        '''
        sci_obs_list = get_sci_obs_list(metadata)
        std_obs_list = get_std_obs_list(metadata)
        ns_proc_list = sci_obs_list + std_obs_list
        for fn in ns_proc_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            sky_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, prev_suffix))
            info_print(f"Subtracting N+S sky frame for {os.path.basename(in_fn)}")
            pywifes.scaled_imarith_mef(in_fn, "-", sky_fn, out_fn, scale="exptime")
        return

    def run_sky_sub(metadata, prev_suffix, curr_suffix, ns=False):
        """
        Subtract sky frames from science objects.

        Parameters
        ----------
        metadata : dict
            Metadata containing information about the observations.
        prev_suffix : str
            Previous suffix of the file names (input).
        curr_suffix : str
            Current suffix of the file names (output).
        ns : bool, optional
            Flag indicating whether to run sky subtraction in non-and-shuffle mode.
            Default is False.

        Returns
        -------
        None
        """
        if ns:
            run_sky_sub_ns(metadata, prev_suffix, curr_suffix)
        else:
            # subtract sky frames from science objects
            for obs in metadata["sci"]:
                if len(obs["sky"]) > 0:
                    sky_fn = obs["sky"][0]
                    sky_proc_fn = os.path.join(
                        out_dir, "%s.p%s.fits" % (sky_fn, prev_suffix)
                    )
                    for fn in obs["sci"]:
                        in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
                        out_fn = os.path.join(
                            out_dir, "%s.p%s.fits" % (fn, curr_suffix)
                        )
                        info_print(f"Subtracting sky frame for {os.path.basename(in_fn)}")
                        # subtract scaled sky framefrom science frame
                        pywifes.scaled_imarith_mef(
                            in_fn, "-", sky_proc_fn, out_fn, scale="exptime"
                        )
                else:
                    for fn in obs["sci"]:
                        in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
                        out_fn = os.path.join(
                            out_dir, "%s.p%s.fits" % (fn, curr_suffix)
                        )
                        info_print(f"Copying science image {os.path.basename(in_fn)}")
                        pywifes.imcopy(in_fn, out_fn)
            # copy stdstar frames
            std_obs_list = get_std_obs_list(metadata)
            for fn in std_obs_list:
                in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
                out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
                if skip_done and os.path.isfile(out_fn):
                    continue
                info_print(f"Copying standard star image {os.path.basename(in_fn)}")
                pywifes.imcopy(in_fn, out_fn)
        return

    # ------------------------------------------------------
    # Image coaddition for science and standards
    # ------------------------------------------------------
    def run_obs_coadd(metadata, prev_suffix, curr_suffix, method="sum", scale=None):
        '''
        Coadd science and standard frames.

        Parameters
        ----------
        metadata : dict
            Dictionary containing metadata information of the FITS files.
        prev_suffix : str
            Previous suffix of the file name (input).
        curr_suffix : str
            Current suffix of the file name (output).
        method : str, optional
            Method used for coadding the frames. Options are "median" or "sum" (default is "sum").
        scale : float, optional
            Scale factor applied to the frames during coaddition (default is None).

        Returns
        -------
        None
        '''
        for obs in metadata["sci"] + metadata["std"]:
            # If just one, then copy it
            if len(obs["sci"]) == 1:
                fn = obs["sci"][0]
                in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
                out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
                if skip_done and os.path.isfile(out_fn):
                    continue
                info_print(f"Copying image {os.path.basename(in_fn)}")
                pywifes.imcopy(in_fn, out_fn)
            # Coadd sci frames
            else:
                in_fn_list = [
                    os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
                    for fn in obs["sci"]
                ]
                out_fn = os.path.join(
                    out_dir, "%s.p%s.fits" % (obs["sci"][0], curr_suffix)
                )
                info_print(f"Coadding images for {os.path.basename(in_fn_list[0])}")
                pywifes.imcombine_mef(in_fn_list, out_fn, scale=scale, method=method)
        return

    # ------------------------------------------------------
    # Flatfield: Response
    # ------------------------------------------------------
    def run_flat_response(metadata, prev_suffix, curr_suffix, mode="all"):
        '''
        Generate the flatfield response function for each flat type, either dome or both (dome and twi).

        Parameters
        ----------
        metadata : dict
            Metadata containing information about the data FITS files.
        prev_suffix : str
            Suffix of the previous data.
        curr_suffix : str
            Suffix of the current data.
        mode : str, optional
            Mode for generating the response function ("dome" or "all"). Default is "all".

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If the requested response mode is not recognized.
        '''

        # Fit the desired style of response function
        info_print("Generating flatfield response function")
        if mode == "all":
            pywifes.wifes_2dim_response(
                super_dflat_mef, 
                super_tflat_mef, 
                flat_resp_fn, 
                wsol_fn=wsol_out_fn, 
                plot=True, 
                plot_dir=plot_dir_arm,
            )
        elif mode == "dome":
            pywifes.wifes_response_poly(
                super_dflat_mef, flat_resp_fn, wsol_fn=wsol_out_fn
            )
        else:
            error_print("Requested response mode not recognized")
            raise ValueError("Requested response mode not recognized")
        return

    # ------------------------------------------------------
    # Flatfield: Division
    # ------------------------------------------------------
    def run_flatfield(metadata, prev_suffix, curr_suffix):
        ''' 
        Apply flat-field correction (division) to science and standard frames.

        Parameters
        ----------
        metadata : dict
            Metadata containing information about the FITS files of the observations.
        prev_suffix : str
            Previous suffix of the file names (input).
        curr_suffix : str
            Current suffix of the file names (output).

        Returns
        -------
        None
        '''
        sci_obs_list = get_primary_sci_obs_list(metadata)
        std_obs_list = get_primary_std_obs_list(metadata)
        info_print(f"Primary science observation list: {sci_obs_list}")
        info_print(f"Primary standard observation list: {std_obs_list}")
        for fn in sci_obs_list + std_obs_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            if skip_done and os.path.isfile(out_fn):
                continue
            info_print(f"Flat-fielding image {os.path.basename(in_fn)}")
            pywifes.imarith_mef(in_fn, "/", flat_resp_fn, out_fn)
        return

    # ------------------------------------------------------
    # Data Cube Generation
    # ------------------------------------------------------
    def run_cube_gen(metadata, prev_suffix, curr_suffix, **args):
        ''' 
        Generate data cubes for science and standard frames.

        Parameters
        ----------
        metadata : dict
            Metadata dictionary containing information about the FITS of the observations.
        prev_suffix : str
            Previous suffix used in the file names (input).
        curr_suffix : str
            Current suffix to be used in the file names (output).
        **args : dict
            Additional keyword arguments to be passed to the `generate_wifes_cube` function. 
            Includes multithread setup or the wavelegnth range.

        Returns
        -------
        None
        '''
        sci_obs_list = get_primary_sci_obs_list(metadata)
        std_obs_list = get_primary_std_obs_list(metadata)
        for fn in sci_obs_list + std_obs_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            if skip_done and os.path.isfile(out_fn):
                continue
            info_print(f"Generating Data Cube for {os.path.basename(in_fn)}")
            # decide whether to use global or local wsol and wire files
            local_wires = get_associated_calib(metadata, fn, "wire")
            if local_wires:
                wire_fn = os.path.join(out_dir, "%s.wire.fits" % (local_wires[0]))
                info_print(f"(Note: using {os.path.basename(wire_fn)} as wire file)")
            else:
                wire_fn = wire_out_fn
            local_arcs = get_associated_calib(metadata, fn, "arc")
            if local_arcs:
                # Do I have two arcs? Do they surround the Science file?
                # Implement linear interpolation:
                if len(local_arcs) == 2:
                    # First, get the Science time
                    f = pyfits.open(in_fn)
                    sci_header = f[0].header
                    sci_time = sci_header["DATE-OBS"]
                    # Now get the arc times
                    arc_times = ["", ""]
                    for i in range(2):
                        local_wsol_out_fn_extra = os.path.join(
                            out_dir, f"{local_arcs[i]}.wsol.fits_extra.pkl")
                        with open(local_wsol_out_fn_extra, "rb") as f:
                            try:
                                f_pickled = pickle.load(f, protocol=2)
                            except:
                                f_pickled = pickle.load(f)  
                        f.close()
                        arc_times[i] = f_pickled[-1][0]

                    # Make sure the Science is between the arcs:
                    t0 = datetime.datetime(
                        np.int(arc_times[0].split("-")[0]),
                        np.int(arc_times[0].split("-")[1]),
                        np.int(arc_times[0].split("-")[2].split("T")[0]),
                        np.int(arc_times[0].split("T")[1].split(":")[0]),
                        np.int(arc_times[0].split(":")[1]),
                        np.int(arc_times[0].split(":")[2].split(".")[0]),
                    )
                    t1 = datetime.datetime(
                        np.int(sci_time.split("-")[0]),
                        np.int(sci_time.split("-")[1]),
                        np.int(sci_time.split("-")[2].split("T")[0]),
                        np.int(sci_time.split("T")[1].split(":")[0]),
                        np.int(sci_time.split(":")[1]),
                        np.int(sci_time.split(":")[2].split(".")[0]),
                    )
                    t2 = datetime.datetime(
                        np.int(arc_times[1].split("-")[0]),
                        np.int(arc_times[1].split("-")[1]),
                        np.int(arc_times[1].split("-")[2].split("T")[0]),
                        np.int(arc_times[1].split("T")[1].split(":")[0]),
                        np.int(arc_times[1].split(":")[1]),
                        np.int(arc_times[1].split(":")[2].split(".")[0]),
                    )
                    ds1 = (t1 - t0).total_seconds()
                    ds2 = (t2 - t1).total_seconds()
                    if ds1 > 0 and ds2 > 0:
                        # Interpolate betweent the two arcs
                        file_camera = sci_header["CAMERA"]
                        if file_camera == "WiFeSRed":
                            w1 = ds1 / (ds1 + ds2)
                            w2 = ds2 / (ds1 + ds2)
                        else:  # file_camera == 'WiFeSBlue'
                            w1 = ds2 / (ds1 + ds2)
                            w2 = ds1 / (ds1 + ds2)

                        # Open the arc solution files
                        fn0 = os.path.join(out_dir, f"{local_arcs[0]}.wsol.fits")
                        fn1 = os.path.join(out_dir, f"{local_arcs[1]}.wsol.fits")
                        fits0 = pyfits.open(fn0)
                        fits1 = pyfits.open(fn1)

                        for i in range(1, len(fits0)):
                            fits0[i].data = w1 * fits0[i].data + w2 * fits1[i].data

                        wsol_fn = os.path.join(out_dir, "%s.wsol.fits" % (fn))
                        fits0.writeto(wsol_fn, overwrite=True)

                        info_print("(2 arcs found)")
                        info_print(f"(Note: using {w1:.2f}x{local_arcs[0]}.wsol.fits + {w2:.2f}x{local_arcs[1]}.wsol.fits as wsol file)")
                            
                    else:
                        # Arcs do not surround the Science frame
                        # Revert to using the first one instead
                        wsol_fn = os.path.join(out_dir, f"{local_arcs[0]}.wsol.fits")
                        info_print("(2 arcs found, but they do not bracket the Science frame!)")
                        info_print(f"(Note: using {os.path.basename(wsol_fn)} as wsol file)")
                else:
                    # IF Either 1 or more than two arcs present, only use the first one.
                    wsol_fn = os.path.join(out_dir, f"{local_arcs[0]}.wsol.fits")
                    info_print(f"(Note: using {os.path.basename(wsol_fn)} as wsol file)")
            else:
                wsol_fn = wsol_out_fn

            # All done, let's generate the cube
            pywifes.generate_wifes_cube(
                in_fn,
                out_fn,
                wire_fn=wire_fn,
                wsol_fn=wsol_fn,
                ny_orig=76,
                offset_orig=2.0,
                **args,
            )
        return

    # ------------------------------------------------------
    # Standard star extraction
    # ------------------------------------------------------
    def run_extract_stars(metadata, prev_suffix, curr_suffix, type="all", **args):
        '''
        Extract standard stars spectrum.

        Parameters
        ----------
        metadata : dict
            The metadata containing information about the FITS files of the observations.
        prev_suffix : str
            The suffix of the previous data files.
        curr_suffix : str
            The suffix of the current data files.
        type : str, optional
            The type of standard stars to extract. Default is "all".
        **args
            Additional keyword arguments.

        Returns
        -------
        None
        '''
        # For each std, extract spectrum as desired
        std_obs_list = get_primary_std_obs_list(metadata, type=type)
        for fn in std_obs_list:
            in_fn = os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
            out_fn = os.path.join(out_dir, f"{fn}.x{prev_suffix}.dat")
            info_print(f"Extract {type} standard star from {os.path.basename(in_fn)}")
            wifes_calib.extract_wifes_stdstar(
                in_fn, save_fn=out_fn, save_mode="ascii", **args
            )
        return
    
    # ------------------------------------------------------
    # Sensitivity Function fit
    # ------------------------------------------------------
    def run_derive_calib(metadata, prev_suffix, curr_suffix, method="poly", **args):
        ''' 
        Derive the sensitivity function from the extracted standard stars.

        Parameters
        ----------
        metadata : dict
            Metadata containing information about the standard stars.
        prev_suffix : str
            Previous suffix of the file names.
        curr_suffix : str
            Current suffix of the file names.
        method : str, optional
            Method for deriving the sensitivity function. Options "smooth_SG" or "poly". 
            Default is "poly".
        **args : dict
            Additional keyword arguments for the calibration function.

        Returns
        -------
        None
        '''
        std_obs_list = get_primary_std_obs_list(metadata, type="flux")

        std_cube_list = [
            os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
            for fn in std_obs_list
        ]
        extract_list = [
            os.path.join(out_dir, f"{fn}.x{prev_suffix}.dat")
            for fn in std_obs_list
        ]
        info_print("Deriving sensitivity function")
        best_calib = wifes_calib.derive_wifes_calibration(
            std_cube_list, calib_fn, extract_in_list=extract_list, method=method, plot_dir=plot_dir_arm, **args
        )
        return

    # ------------------------------------------------------
    # Flux Calibration
    # ------------------------------------------------------
    def run_flux_calib(metadata, prev_suffix, curr_suffix, mode="pywifes", **args):
            '''
            Flux calibrate all science and standard observations.

            Parameters
            ----------
            metadata : dict
                Metadata containing information about the observations.
            prev_suffix : str
                Previous suffix of the file name (input).
            curr_suffix : str
                Current suffix of the file name (output).
            mode : str, optional
                Calibration mode. Options "iraf" or "pywifes". 
                Default is "pywifes".
            **args : dict
                Additional keyword arguments.

            Returns
            -------
            None
            '''
            sci_obs_list = get_primary_sci_obs_list(metadata)
            std_obs_list = get_primary_std_obs_list(metadata)
            for fn in sci_obs_list + std_obs_list:
                in_fn = os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
                out_fn = os.path.join(out_dir, f"{fn}.p{curr_suffix}.fits")
                info_print(f"Flux-calibrating cube {os.path.basename(in_fn)}")
                wifes_calib.calibrate_wifes_cube(in_fn, out_fn, calib_fn, mode)
            return

    # ------------------------------------------------------
    # Telluric - derive
    # ------------------------------------------------------
    def run_derive_telluric(metadata, prev_suffix, curr_suffix, **args):
        '''
        Derive the telluric correction from the standard star.

        Parameters
        ----------
        metadata : dict
            Metadata containing information about the observations.
        prev_suffix : str
            Previous suffix used in the file names.
        curr_suffix : str
            Current suffix to be used in the file names.
        **args : dict
            Additional keyword arguments for plotting options.

        Returns
        -------
        None
        '''
        std_obs_list = get_primary_std_obs_list(metadata, "telluric")
        std_cube_list = [
            os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
            for fn in std_obs_list
        ]
        extract_list = [
            os.path.join(out_dir, f"{fn}.x{prev_suffix}.dat")
            for fn in std_obs_list
        ]
        info_print("Deriving telluric correction")
        wifes_calib.derive_wifes_telluric(
            std_cube_list, tellcorr_fn, extract_in_list=extract_list, plot_dir=plot_dir_arm, **args
        )
        return

    def run_telluric_corr(metadata, prev_suffix, curr_suffix, **args):
            '''
            Apply telluric correction for all science and standard observations.

            Parameters
            ----------
            metadata : dict
                Metadata containing information about the observations.
            prev_suffix : str
                Previous suffix of the file name (input).
            curr_suffix : str
                Current suffix of the file name (output). 
            **args : dict
                Additional keyword arguments.

            Returns
            -------
            None
            '''
            sci_obs_list = get_primary_sci_obs_list(metadata)
            std_obs_list = get_primary_std_obs_list(metadata)
            for fn in sci_obs_list + std_obs_list:
                in_fn = os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
                out_fn = os.path.join(out_dir, f"{fn}.p{curr_suffix}.fits")
                info_print(f"Correcting telluric in {os.path.basename(in_fn)}")
                wifes_calib.apply_wifes_telluric(in_fn, out_fn, tellcorr_fn)
            return

    def run_save_3dcube(metadata, prev_suffix, curr_suffix, **args):
        '''
        Save 3D Data Cube for all science and standard observations.
        Final data product.

        Parameters
        ----------
        metadata : dict
            Metadata containing information about the observations.
        prev_suffix : str
            Previous suffix of the file name (input).
        curr_suffix : str
            Current suffix of the file name.
            Saves final "*.cube.fits" file (output).
        **args : dict
            Additional keyword arguments to be passed to the `generate_wifes_3dcube` function.

        Returns
        -------
        None
        '''
        sci_obs_list = get_primary_sci_obs_list(metadata)
        std_obs_list = get_primary_std_obs_list(metadata)

        # Check if is half-frame from the first sci image
        sci_filename = temp_data_dir +  sci_obs_list[0] + ".fits"

        halfframe = is_halfframe(sci_filename)
        
        # now generate cubes
        for fn in sci_obs_list + std_obs_list:
            in_fn = os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
            out_fn = os.path.join(out_dir, f"{fn}.{curr_suffix}.fits")
            info_print(f"Saving 3D Data Cube for {os.path.basename(in_fn)}")
            pywifes.generate_wifes_3dcube(in_fn, out_fn, halfframe=halfframe, **args)
        return

    # --------------------------------------------
    # INICIATE THE SCRIPT
    # --------------------------------------------

    # Initialize ArgumentParser with a description
    parser = argparse.ArgumentParser(
        description="The Python data reduction pipeline for WiFeS."
    )

    # The raw data directory is a required positional argument
    parser.add_argument("user_data_dir", type=str, help="Path to the raw data directory.")

    # Option for specifying the path to the red parameters JSON file
    parser.add_argument(
        "--red-params",
        type=str,
        help="Optional: Path to the configuration JSON file containing parameters for reducing the blue arm.",
    )

    # Option for specifying the path to the blue parameters JSON file
    parser.add_argument(
        "--blue-params",
        type=str,
        help="Optional: Path to the configuration JSON file containing parameters for reducing the blue arm.",
    )

    # Option for triggering the reduction from master calibrations
    parser.add_argument(
        "--from-master",
        type=str,
        const="./data_products/master_calib",
        nargs="?",
        help="Optional: Path to the master calibrations directory. If not provided, the default path will be used: .",
    )

    # Option for specifying to skip already completed steps
    parser.add_argument(
        "-skip-done",
        action="store_true",
        help="Optional: Skip already completed steps.",
    )

    parser.add_argument(
        "-just-calib",
        action="store_true",
        help="Optional: Only basics master calibration files will produced.",
    )

    args = parser.parse_args()

    # Validate and process the user_data_dir
    user_data_dir = os.path.abspath(args.user_data_dir)
    if not user_data_dir.endswith("/"):
        user_data_dir += "/"
    info_print(f"Processing data in directory: {user_data_dir}")

    # Handling reduction parameters.
    params_path = {
        "blue": None,
        "red": None,
    }

    # Red
    if args.red_params:
        params_path["red"] = os.path.abspath(args.red_params)
        info_print(f"Using red parameters from: {params_path['red']}")

    # Blue
    if args.blue_params:
        params_path["blue"] = os.path.abspath(args.blue_params)
        info_print(f"Using blue parameters from: {params_path['blue']}")

    # Reduction from master calibration frames
    from_master = args.from_master

    # Only basics master calibration.
    just_calib = args.just_calib

    # Set to skip already done files
    skip_done = args.skip_done

    # Set paths
    reduction_scripts_dir = os.path.dirname(__file__)
    working_dir = os.getcwd()

    # Creates a temporary data directory containning all raw data for reduction.
    temp_data_dir = os.path.join(working_dir, f"data_products/intermediate/raw_data_temp/")
    os.makedirs(temp_data_dir, exist_ok=True)

    all_fits_names = get_file_names(user_data_dir, "*.fits")
    # Copy raw data  from user's direcory into temporaty raw directory.
    copy_files(user_data_dir, temp_data_dir, all_fits_names)


    # Creates a directory for plot.
    plot_dir = os.path.join(working_dir, f"data_products/plots/")
    os.makedirs(plot_dir, exist_ok=True)


    # Classify all raw data (red and blue arm)
    obs_metadatas = classify(temp_data_dir)

    # Set grism_key dictionary due to different keyword names for red and blue arms.
    grism_key = {
        "blue": "GRATINGB",
        "red": "GRATINGR",
    }

    # Set the directory for master calibration files (default is ./data_products/master_calib/)
    # Define a list of steps to skip if the reduction is being performed using master calibration files
    if from_master:
        master_dir = os.path.abspath(from_master)
        extra_skip_steps = [
            "superbias",
            "superflat",
            "slitlet_profile",
            "flat_cleanup",
            "superflat_mef",
            "wave_soln",
            "wire_soln",
            "flat_response",
        ]

    else:
        # Master calibration files firectory
        master_dir = os.path.join(working_dir, "data_products/master_calib/")
        os.makedirs(master_dir, exist_ok=True)
        # No extra skiped steps in principal.
        extra_skip_steps = []

    info_print(f"Processing using master calibrations from: '{master_dir}'.")

    for arm in obs_metadatas.keys():
        
        try:
            # ------------------------------------------------------------------------
            #      LOAD JSON FILE WITH USER DATA REDUCTION SETUP
            # ------------------------------------------------------------------------
            obs_metadata = obs_metadatas[arm]
            # Determine the grism and observing mode used in the first image of science, standard, or arc of the respective arm.
            # Skip reductions steps if no objects no standar star observations are present.

            if obs_metadata["sci"]:
                reference_filename = obs_metadata["sci"][0]["sci"][0] + ".fits"
            elif obs_metadata["std"]:
                reference_filename = obs_metadata["std"][0]["sci"][0] + ".fits"
            elif obs_metadata["arc"]:
                reference_filename = obs_metadata["arc"][0] + ".fits"
            else:
                error_print("No science, standard, or arc files found in metadata.")
                raise ValueError("No science, standard, or arc files found in metadata.")

            # Check observing mode 
            if pywifes.is_nodshuffle(temp_data_dir + reference_filename):
                obs_mode = "ns"

            elif pywifes.is_subnodshuffle(temp_data_dir + reference_filename):
                obs_mode = "ns"
            else:
                obs_mode = "class"

            # Check if is half-frame
            halfframe = is_halfframe(temp_data_dir + reference_filename)
            if halfframe:
                obs_metadata = calib_to_half_frame(obs_metadata,temp_data_dir) 


            # Grism
            grism = pyfits.getheader(temp_data_dir + reference_filename)[grism_key[arm]]

            # Set the JSON file path and read it.
            if just_calib and (params_path[arm] is None):
                json_path = f"./pipeline_params/just-calib/{arm}/params_{obs_mode}_{grism}.json"
            elif params_path[arm] is None:
                json_path = f"./pipeline_params/{arm}/params_{obs_mode}_{grism}.json"
            else:
                json_path = params_path[arm]

            # Load the JSON file
            proc_steps = load_config_file(json_path)

            # Create data products directory structure
            out_dir = os.path.join(working_dir, f"data_products/intermediate/{arm}")
            os.makedirs(out_dir, exist_ok=True)

            calib_prefix = f"wifes_{arm}"

            # Creates a directory for diagnisis plots (one per arm).
            plot_dir_arm = os.path.join(plot_dir, arm)
            os.makedirs(plot_dir_arm, exist_ok=True)



            # ------------------------------------------------------------------------
            # Define names for master calibration files and set their path.
            # ------------------------------------------------------------------------
            # Bias Master Files
            superbias_fn = os.path.join(master_dir, "%s_superbias.fits" % calib_prefix)
            superbias_fit_fn = os.path.join(master_dir, "%s_superbias_fit.fits" % calib_prefix)

            # Flat Master Files
            # Dome
            super_dflat_raw = os.path.join(master_dir, "%s_super_domeflat_raw.fits" % calib_prefix)
            super_dflat_fn = os.path.join(master_dir, "%s_super_domeflat.fits" % calib_prefix)
            super_dflat_mef = os.path.join(master_dir, "%s_super_domeflat_mef.fits" % calib_prefix)
            # Twilight
            super_tflat_raw = os.path.join(master_dir, "%s_super_twiflat_raw.fits" % calib_prefix)
            super_tflat_fn = os.path.join(master_dir, "%s_super_twiflat.fits" % calib_prefix)
            super_tflat_mef = os.path.join(master_dir, "%s_super_twiflat_mef.fits" % calib_prefix)
            # Slitlet definition
            slitlet_def_fn = os.path.join(master_dir, "%s_slitlet_defs.pkl" % calib_prefix)
            wsol_out_fn = os.path.join(master_dir, "%s_wave_soln.fits" % calib_prefix)
            wire_out_fn = os.path.join(master_dir, "%s_wire_soln.fits" % calib_prefix)
            flat_resp_fn = os.path.join(master_dir, "%s_resp_mef.fits" % calib_prefix)
            calib_fn = os.path.join(master_dir, "%s_calib.pkl" % calib_prefix)
            tellcorr_fn = os.path.join(master_dir, "%s_tellcorr.pkl" % calib_prefix)

            # When reducing from master calibration files, if the tellic_correction file is already among the master calibrations, skip its generation.
            if from_master and os.path.exists(tellcorr_fn):
                extra_skip_steps.append("derive_calib")

            # ------------------------------------------------------------------------
            # Run proccessing steps
            # ------------------------------------------------------------------------
            info_print(f"________________________________________________________________")
            info_print(f"Starting processing of {arm} arm")
            info_print(f"________________________________________________________________")

            prev_suffix = None
            for step in proc_steps[arm]:
                step_name = step["step"]
                step_run = step["run"]
                step_suffix = step["suffix"]
                step_args = step["args"]
                func_name = "run_" + step_name
                func = locals()[func_name]

                # When master calibrations are in use, the steps listed in skip_steps will be skipped.
                if step_name in extra_skip_steps:
                    info_print('======================')
                    info_print(f"Skipping step: {step_name}")
                    info_print('======================')
                    continue

                if step_run:
                    info_print('======================')
                    info_print(step_name)
                    info_print('======================')

                    func(
                        obs_metadata,
                        prev_suffix=prev_suffix,
                        curr_suffix=step_suffix,
                        **step_args,
                    )
                    if step_suffix != None:
                        prev_suffix = step_suffix

                else:
                    pass

            # Extra Data Reduction Quality Plots:

            # Dome Flats    
            try:
                title = 'Dome Flatfield'            
                output_plot = os.path.join(plot_dir_arm, "raw_domeflat_check.png")

                flat_image_path = os.path.join(master_dir, f"wifes_{arm}_super_domeflat_raw.fits")
                slitlet_path = os.path.join(master_dir, f"wifes_{arm}_slitlet_defs.pkl")
                flatfiled_plot(flat_image_path, slitlet_path, title, output_plot)
            except:
                pass

            # Twilight Flats
            try:
                title = 'Twilight Flatfield'            
                output_plot = os.path.join(plot_dir_arm, "raw_twiflat_check.png")
                flat_image_path = os.path.join(master_dir, f"wifes_{arm}_super_twiflat_raw.fits")
                slitlet_path = os.path.join(master_dir, f"wifes_{arm}_slitlet_defs.pkl")
                flatfiled_plot(flat_image_path, slitlet_path, title, output_plot)
            except:
                pass

        except Exception as exc:
            warning_print("________________________________________________________________")
            warning_print(f"{arm} arm skipped, an error occurred during processing: '{exc}'.")        
            warning_print("________________________________________________________________")


    # Delete temporary directory containing raw data.
    shutil.rmtree(temp_data_dir)

    # ----------------------------------------------------------
    # Move reduce cube to the data_products directory
    # ----------------------------------------------------------
    if just_calib:
        info_print("Only basics master calibration files have been produced.")
    else:
        destination_dir = os.path.join(working_dir, "data_products")
        info_print(f"Moving reduced 3D cubes to {destination_dir}.")

        # Red
        red_cubes_path = os.path.join(working_dir, "data_products/intermediate/red/")
        red_cubes_file_name = get_file_names(red_cubes_path, "*.cube.fits")
        # Move reduced cubes to the data_product
        move_files(red_cubes_path, destination_dir, red_cubes_file_name)

        # Blue
        blue_cubes_path = os.path.join(working_dir, "data_products/intermediate/blue/")
        blue_cubes_file_name = get_file_names(blue_cubes_path, "*.cube.fits")
        # Move reduced cubes to the data_product
        move_files(blue_cubes_path, destination_dir, blue_cubes_file_name)

        # ----------------------------------------------------------
        # Find and list all reduced cubes in the destination directory
        # ----------------------------------------------------------
        reduced_cubes_paths = [
            os.path.join(destination_dir, file_name) for file_name in blue_cubes_file_name
        ] + [os.path.join(destination_dir, file_name) for file_name in red_cubes_file_name]

        # ----------------------------------------------------------
        # Match cubes from the same observation based on DATE-OBS
        # ----------------------------------------------------------
        matched_cubes = cube_matcher(reduced_cubes_paths)

        # ----------------------------------------------------------
        # Read extraction parameters from JSON file
        # ----------------------------------------------------------
        extract_params = load_config_file(f"./pipeline_params/params_extract_{obs_mode}.json")

        # ----------------------------------------------------------
        # Loop over matched cubes list
        # ----------------------------------------------------------
        for match_cubes in matched_cubes:
            # ----------
            # Extraction
            # ----------
            info_print('======================')
            info_print('Extracting spectra')
            info_print('======================')
            blue_cube_path = match_cubes["Blue"]
            red_cube_path = match_cubes["Red"]
            plot_name = match_cubes["file_name"].replace(".cube", "_detection_plot.png")
            plot_path = os.path.join(plot_dir,plot_name)
            plot = extract_params["plot"]
            # Run auto-extraction
            detect_extract_and_save(
                blue_cube_path,
                red_cube_path,
                destination_dir,
                r_arcsec=extract_params["r_arcsec"],
                border_width=extract_params["border_width"],
                sky_sub=extract_params["sky_sub"],
                plot=plot,
                plot_path=plot_path,
            )


            # ------------------------------------
            # Splice only paired cubes and spectra
            # ------------------------------------

            if blue_cube_path is not None and red_cube_path is not None:
                info_print('======================')
                info_print('Splicing blue and red cubes')
                info_print('======================')
                
                blue_cube_name = os.path.basename(blue_cube_path)
                red_cube_name = os.path.basename(red_cube_path)

                # Get filename of form `xxx-Splice-UTxxx.cube.fits`
                spliced_cube_name = blue_cube_name.replace("Blue", "Splice")
                spliced_cube_path = os.path.join(destination_dir, spliced_cube_name)

                # Splice cubes
                splice_cubes(match_cubes["Blue"], match_cubes["Red"], spliced_cube_path)

                # Find blue spectra files matching the pattern 'xxx-Blue-UTxxx.spec.ap*'
                pattern_blue = os.path.join(
                    destination_dir, blue_cube_name.replace("cube", "spec.ap*")
                )
                blue_specs = glob.glob(pattern_blue)

                # Find red spectra files matching the pattern 'xxx-Red-UTxxx.spec.ap*'
                pattern_red = os.path.join(
                    destination_dir, red_cube_name.replace("cube", "spec.ap*")
                )
                red_specs = glob.glob(pattern_red)

                # Splice spectra
                for blue_spec, red_spec in zip(blue_specs, red_specs):
                    # Generate filename for spliced spectrum 'xxx-Splice-UTxxx.spec.apx.fits'
                    spliced_spectrum_name = os.path.basename(blue_spec).replace(
                        "Blue", "Splice"
                    )
                    output = os.path.join(
                        working_dir, destination_dir, spliced_spectrum_name
                    )
                    splice_spectra(blue_spec, red_spec, output)


            # Plot extracted spectra:
            if plot:
                if blue_cube_path is not None and red_cube_path is not None:

                    blue_pattern = blue_cube_path.replace("cube.fits", "spec.ap*")
                    blue_specs = glob.glob(blue_pattern)

                    red_pattern = red_cube_path.replace("cube.fits", "spec.ap*")
                    red_specs = glob.glob(red_pattern)

                    for blue_spec, red_spec in zip(blue_specs, red_specs):
                        spliced_spec = blue_spec.replace("Blue", "Splice")
                        
                        plot_1D_spectrum(blue_spec,plot_dir=plot_dir)
                        plot_1D_spectrum(red_spec,plot_dir=plot_dir)
                        plot_1D_spectrum(spliced_spec,plot_dir=plot_dir)


                elif blue_cube_path is not None:

                    blue_pattern = blue_cube_path.replace("cube.fits", "spec.ap*")
                    blue_specs = glob.glob(blue_pattern)


                    for blue_spec in blue_specs:

                        plot_1D_spectrum(blue_spec,plot_dir=plot_dir)



                elif red_cube_path is not None:

                    red_pattern = red_cube_path.replace("cube.fits", "spec.ap*")
                    red_specs = glob.glob(red_pattern)


                    for red_spec in red_specs:
                                            
                        plot_1D_spectrum(red_spec,plot_dir=plot_dir)



    # ----------------------------------------------------------
    # Print total running time
    # ----------------------------------------------------------
    end_time = datetime.datetime.now()
    duration = end_time - start_time
    messagge = "All done in %.01f seconds." % duration.total_seconds()
    info_print(messagge)
    print('\U0001F52D',messagge,'\u2B50')

if __name__ == "__main__":
    main()
