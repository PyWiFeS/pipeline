#! /usr/bin/env python3

import argparse
from astropy.io import fits as pyfits
import datetime
import gc
import glob
import logging
import numpy
import os
import pickle
import pyjson5
import re
import shutil

from pywifes import pywifes
from pywifes import wifes_calib
from pywifes import wifes_wsol
from pywifes.data_classifier import classify, cube_matcher
from pywifes.extract_spec import detect_extract_and_save, plot_1D_spectrum
from pywifes.lacosmic import lacos_wifes
from pywifes.logger_config import setup_logger, custom_print
from pywifes.splice import splice_spectra, splice_cubes
from pywifes.wifes_utils import is_halfframe, is_nodshuffle, is_subnodshuffle, is_taros

# Set paths
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
# logging.captureWarnings(True)

info_print("Starting PyWiFeS data reduction pipeline thread.")


# ------------------------------------------------------------------------
# Function definition
# ------------------------------------------------------------------------

def move_files(src_dir_path, destination_dir_path, filenames):
    try:
        for file in filenames:
            src_file = os.path.join(src_dir_path, file)
            dest_file = os.path.join(destination_dir_path, file)
            debug_print(f"Moving file {src_file} to {dest_file}")
            shutil.move(src_file, dest_file)
    except Exception as e:
        error_print(f"Error moving files: {e}")


def copy_files(src_dir_path, destination_dir_path, filenames):
    try:
        for file in filenames:
            src_file = os.path.join(src_dir_path, file)
            dest_file = os.path.join(destination_dir_path, file)
            # Handle common file compression for raw data
            if re.search("\\.fz$", src_file):
                dest_file = dest_file.rstrip(".fz")
                if os.path.isfile(dest_file) \
                        and os.path.getmtime(src_file) < os.path.getmtime(dest_file):
                    continue
                temph = pyfits.open(src_file)
                pyfits.writeto(dest_file, data=temph[1].data, header=temph[1].header,
                               output_verify="fix", overwrite=True)
                temph.close()
            elif re.search("\\.gz$", src_file):
                dest_file = dest_file.rstrip(".gz")
                if os.path.isfile(dest_file) \
                        and os.path.getmtime(src_file) < os.path.getmtime(dest_file):
                    continue
                temph = pyfits.open(src_file)
                pyfits.writeto(dest_file, data=temph[0].data, header=temph[0].header,
                               output_verify="fix", overwrite=True)
                temph.close()
            else:
                if os.path.isfile(dest_file) \
                        and os.path.getmtime(src_file) < os.path.getmtime(dest_file):
                    continue
                shutil.copy(src_file, dest_file)
    except Exception as e:
        error_print(f"Error copying files: {e}")


def get_file_names(src_dir_path, glob_pattern):
    filepaths = glob.glob(os.path.join(src_dir_path, glob_pattern))
    names = []
    for filepath in filepaths:
        filename = os.path.basename(filepath)
        names.append(filename)
    return names


def load_config_file(filename):
    reduction_scripts_dir = os.path.dirname(__file__)
    file_path = os.path.join(reduction_scripts_dir, filename)
    info_print(f"Loading configuration file: {file_path}")
    with open(file_path, "r") as f:
        return pyjson5.load(f)

# ------------------------------------------------------------------------


def main():
    start_time = datetime.datetime.now()
    info_print(f"Pipeline started at {start_time}")

    # ------------------------------------------------------------------------
    # METADATA WRANGLING FUNCTIONS
    # ------------------------------------------------------------------------
    def get_full_obs_list(metadata, exclude=None):
        keep_type = ["bias", "arc", "wire", "dark", "domeflat", "twiflat"]
        # Allow types to be excluded
        if exclude is not None:
            keep_type = [t for t in keep_type if t not in exclude]
        full_obs_list = []
        base_fn_list = [fn for t in keep_type for fn in metadata[t]]
        for fn in base_fn_list:
            if fn not in full_obs_list:
                full_obs_list.append(fn)

        keep_type = ["sci", "std"]
        # Allow types to be excluded
        if exclude is not None:
            keep_type = [t for t in keep_type if t not in exclude]
        for obs in [fdict for t in keep_type for fdict in metadata[t]]:
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
        - metadata (dict): The metadata containing information about the observations.
        - type (str, optional): The type of observations to include in the list. Default is "all".

        Returns:
        - std_obs_list (list): A list of standard observation filenames.

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
            metadata (dict): The metadata containing information about the observations.
            type (str, optional): The type of standard star observations to include in the list.
                Possible values are "all", "telluric", or "flux". Defaults to "all".

        Returns:
            list: The list of primary standard observations.

        Raises:
            ValueError: If the standard star type is not understood.

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
    # DEFINE THE PROCESSING STEPS
    # ------------------------------------------------------------------------
    # Overscan subtraction
    # ------------------------------------------------------------------------
    def run_overscan_sub(metadata, prev_suffix, curr_suffix, poly_high_oscan=True, **args):
        full_obs_list = get_full_obs_list(metadata)
        first = True
        if not poly_high_oscan:
            first = False
            oscanmask = None
        for fn in full_obs_list:
            in_fn = os.path.join(temp_data_dir, "%s.fits" % fn)
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            if skip_done and os.path.isfile(out_fn):
                # cannot check mtime here because of fresh copy to raw_data_temp
                continue
            info_print(f"Subtracting Overscan for {os.path.basename(in_fn)}")
            if first:
                # Find a domeflat to generate mask for overscan
                if metadata["domeflat"]:
                    dflat = os.path.join(temp_data_dir, "%s.fits" % metadata["domeflat"][0])
                    pywifes.make_overscan_mask(dflat, omask=overscanmask_fn, data_hdu=my_data_hdu)
                    oscanmask = overscanmask_fn
                else:
                    oscanmask = None
                first = False
            pywifes.subtract_overscan(in_fn, out_fn, data_hdu=my_data_hdu, omaskfile=oscanmask, **args)
        return

    # ------------------------------------------------------
    # repair bad pixels!
    # ------------------------------------------------------
    def run_bpm_repair(metadata, prev_suffix, curr_suffix, **args):
        full_obs_list = get_full_obs_list(metadata)
        for basename in full_obs_list:
            input_filename = f"{basename}.p{prev_suffix}.fits"
            output_filename = f"{basename}.p{curr_suffix}.fits"
            input_filepath = os.path.join(out_dir, input_filename)
            output_filepath = os.path.join(out_dir, output_filename)
            if skip_done and os.path.isfile(output_filepath) \
                    and os.path.getmtime(input_filepath) < os.path.getmtime(output_filepath):
                continue
            info_print(f"Repairing {arm} bad pixels for {input_filename}")
            if arm == "red":
                pywifes.repair_red_bad_pix(
                    input_filepath, output_filepath, data_hdu=my_data_hdu, **args
                )
            if arm == "blue":
                pywifes.repair_blue_bad_pix(
                    input_filepath, output_filepath, data_hdu=my_data_hdu, **args
                )

    # ------------------------------------------------------
    # Generate super-bias
    # ------------------------------------------------------
    def run_superbias(metadata, prev_suffix, curr_suffix, method="row_med", **args):
        '''Generate superbias for the entire dataset and for each science frame.
        Fit a smart surface to the bias or take the median of each row.
        '''
        bias_list = [os.path.join(out_dir, "%s.p%s.fits" % (x, prev_suffix))
                     for x in metadata["bias"]
                     ]
        if not (skip_done and os.path.isfile(superbias_fn) and os.path.isfile(superbias_fit_fn)
                and os.path.getmtime(superbias_fn) < os.path.getmtime(superbias_fit_fn)):
            info_print("Calculating Global Superbias")
            pywifes.imcombine(bias_list, superbias_fn, data_hdu=my_data_hdu,
                              kwstring="BIASN", commstring="bias")
            if method == "fit" or method == "row_med":
                pywifes.generate_wifes_bias_fit(
                    superbias_fn,
                    superbias_fit_fn,
                    data_hdu=my_data_hdu,
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
                local_superbias = os.path.join(
                    out_dir, "%s.fits" % (local_bias_fn + ".lsb")
                )
                local_superbias_fit = os.path.join(
                    out_dir, "%s.fits" % (local_bias_fn + ".lsb_fit")
                )
                if skip_done and os.path.isfile(local_superbias) and os.path.isfile(local_superbias_fit) \
                        and os.path.getmtime(local_superbias) < os.path.getmtime(local_superbias_fit):
                    continue
                info_print(f"Calculating Local Superbias for {local_bias_fn}")
                # step 1 - coadd biases
                local_biases_filename = [
                    os.path.join(out_dir, "%s.p%s.fits" % (x, prev_suffix))
                    for x in local_biases
                ]
                pywifes.imcombine(
                    local_biases_filename, local_superbias, data_hdu=my_data_hdu,
                    kwstring="LOCBN", commstring="local bias"
                )
                # step 2 - generate fit
                if method == "fit" or method == "row_med":
                    pywifes.generate_wifes_bias_fit(
                        local_superbias,
                        local_superbias_fit,
                        data_hdu=my_data_hdu,
                        method=method,
                        **args,
                    )
                else:
                    pywifes.imcopy(local_superbias, local_superbias_fit)
        return

    # ----------------------------------------------------
    # Subtract bias
    # ----------------------------------------------------
    def run_bias_sub(metadata, prev_suffix, curr_suffix, method="subtract"):
        full_obs_list = get_full_obs_list(metadata)
        for fn in full_obs_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                continue
            local_biases = get_associated_calib(metadata, fn, "bias")
            if local_biases:
                local_bias_fn = get_associated_calib(metadata, fn, "bias")[0]
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
            elif method == "subtract":
                pywifes.imarith(in_fn, "-", bias_fit_fn, out_fn, data_hdu=my_data_hdu)
            else:
                raise ValueError(f"Unknown bias_sub method '{method}'. Options: 'subtract', 'copy'.")
        return

    # ------------------------------------------------------
    # Generate super-flat/wire/arc
    # ------------------------------------------------------
    def run_superflat(
        metadata, prev_suffix, curr_suffix, source, **args
    ):
        kwstring = None
        commstring = None
        outvarimg = None
        save_prefix = 'imcombine_inputs'
        if source == "dome":
            out_fn = super_dflat_raw
            flat_list = [
                os.path.join(out_dir, "%s.p%s.fits" % (x, prev_suffix))
                for x in metadata["domeflat"]
            ]
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(flat_list[0]) < os.path.getmtime(out_fn):
                return
            info_print(f"List of {source} flats: {flat_list}")
            kwstring = "FLATN"
            commstring = "lamp flat"
            outvarimg = os.path.join(master_dir, f"wifes_{arm}_super_domeflat_raw_var.fits")
            save_prefix = f"dflat_{save_prefix}"
        elif source == "twi":
            out_fn = super_tflat_raw
            flat_list = [
                os.path.join(out_dir, "%s.p%s.fits" % (x, prev_suffix))
                for x in metadata["twiflat"]
            ]
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(flat_list[0]) < os.path.getmtime(out_fn):
                return
            info_print(f"List of {source} flats: {flat_list}")
            kwstring = "TWIN"
            commstring = "twilight flat"
            save_prefix = f"twi_{save_prefix}"
        elif source == "wire":
            out_fn = super_wire_raw
            flat_list = [
                os.path.join(out_dir, "%s.p%s.fits" % (x, prev_suffix))
                for x in metadata["wire"]
            ]
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(flat_list[0]) < os.path.getmtime(out_fn):
                return
            info_print(f"List of wire frames: {flat_list}")
            kwstring = "WIREN"
            commstring = "wire"
        elif source == "arc":
            out_fn = super_arc_raw
            flat_list = [
                os.path.join(out_dir, "%s.p%s.fits" % (x, prev_suffix))
                for x in metadata["arc"]
            ]
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(flat_list[0]) < os.path.getmtime(out_fn):
                return
            info_print(f"List of arc frames: {flat_list}")
            kwstring = "ARCN"
            commstring = "arc"
        else:
            error_print(f"Calibration type '{source}' not recognised")
            raise ValueError(f"Calibration type '{source}' not recognised")
        if not flat_list:
            warning_print(f"No {source} flats found. Skipping the superflat generation for {source}.")
            return
        info_print(f"Generating co-add {source} flat")
        pywifes.imcombine(
            flat_list, out_fn, data_hdu=my_data_hdu, kwstring=kwstring, commstring=commstring,
            outvarimg=outvarimg, plot_dir=plot_dir_arm, save_prefix=save_prefix, **args
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
        '''
        # check the slitlet definition file
        if os.path.isfile(slitlet_def_fn):
            slitlet_fn = slitlet_def_fn
        else:
            slitlet_fn = None
        if "dome" in type:
            if skip_done and os.path.isfile(super_dflat_fn) \
                    and os.path.getmtime(super_dflat_raw) < os.path.getmtime(super_dflat_fn):
                return
            if os.path.isfile(super_dflat_raw):
                info_print(f"Correcting master domeflat {os.path.basename(super_dflat_raw)}")
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
            if skip_done and os.path.isfile(super_tflat_fn) \
                    and os.path.getmtime(super_tflat_raw) < os.path.getmtime(super_tflat_fn):
                return
            if os.path.isfile(super_tflat_raw):
                info_print(f"Correcting master twilight flat {os.path.basename(super_tflat_raw)}")
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
        '''
        output_fn = slitlet_def_fn
        if os.path.isfile(super_dflat_fn):
            flatfield_fn = super_dflat_fn
        else:
            flatfield_fn = super_dflat_raw
        if skip_done and os.path.isfile(output_fn) \
                and os.path.getmtime(flatfield_fn) < os.path.getmtime(output_fn):
            return
        pywifes.derive_slitlet_profiles(
            flatfield_fn, output_fn, data_hdu=my_data_hdu, **args
        )
        return

    # ------------------------------------------------------
    # Create MEF files
    # ------------------------------------------------------
    def run_superflat_mef(metadata, prev_suffix, curr_suffix, source, **args):
        if source == "dome":
            if os.path.isfile(super_dflat_fn):
                in_fn = super_dflat_fn
            elif os.path.isfile(super_dflat_raw):
                in_fn = super_dflat_raw
            else:
                warning_print("No master dome flat found. Skipping MEF generation for dome flat.")
                return
            out_fn = super_dflat_mef

        elif source == "twi":
            if os.path.isfile(super_tflat_fn):
                in_fn = super_tflat_fn
            elif os.path.isfile(super_tflat_raw):
                in_fn = super_tflat_raw
            else:
                warning_print("No master twilight flat found. Skipping MEF generation for twilight flat.")
                return
            out_fn = super_tflat_mef

        elif source == "wire":
            if os.path.isfile(super_wire_raw):
                in_fn = super_wire_raw
            else:
                warning_print("No master wire frame found. Skipping MEF generation for wire.")
                return
            out_fn = super_wire_mef

        elif source == "arc":
            if os.path.isfile(super_arc_raw):
                in_fn = super_arc_raw
            else:
                warning_print("No master arc frame found. Skipping MEF generation for arc.")
                return
            out_fn = super_arc_mef

        else:
            error_print(f"Calibration type '{source}' not recognised")
            raise ValueError(f"Calibration type '{source}'' not recognised")

        # check the slitlet definition file
        if os.path.isfile(slitlet_def_fn):
            slitlet_fn = slitlet_def_fn
        else:
            slitlet_fn = None
        # run it!
        if skip_done and os.path.isfile(out_fn) \
                and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
            return
        info_print(f"Generating MEF {source} flat")
        pywifes.wifes_slitlet_mef(
            in_fn, out_fn, data_hdu=my_data_hdu, slitlet_def_file=slitlet_fn, **args
        )
        return

    def run_slitlet_mef(metadata, prev_suffix, curr_suffix, **args):
        # Do not need MEF versions of individual frames for those with supercals (and with no "locals")
        excl_list = ["bias", "dark", "domeflat", "twiflat"]
        full_obs_list = get_full_obs_list(metadata, exclude=excl_list)
        sci_obs_list = get_sci_obs_list(metadata)
        std_obs_list = get_std_obs_list(metadata)
        ns_proc_list = sci_obs_list + std_obs_list
        # check the slitlet definition file
        if os.path.isfile(slitlet_def_fn):
            slitlet_fn = slitlet_def_fn
        else:
            slitlet_fn = None
        for fn in full_obs_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                continue
            info_print(f"Creating MEF file for {os.path.basename(in_fn)}")
            if obs_mode == "ns" and fn in ns_proc_list:
                sky_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, curr_suffix))
                pywifes.wifes_slitlet_mef_ns(
                    in_fn,
                    out_fn,
                    sky_fn,
                    data_hdu=my_data_hdu,
                    slitlet_def_file=slitlet_fn,
                    **args
                )
            else:
                pywifes.wifes_slitlet_mef(
                    in_fn, out_fn, data_hdu=my_data_hdu, slitlet_def_file=slitlet_fn, **args
                )
            gc.collect()
        return

    # ------------------------------------------------------
    # Wavelength solution
    # ------------------------------------------------------
    def run_wave_soln(metadata, prev_suffix, curr_suffix, **args):
        '''
        Generate the master arc solution, based on generic arcs at first.
        Then looks for the local wavelength solutions for science or standards (sky not required at this stage).
        Check if the file has a dedicated arc associated with it.
        If two arcs are present, find a solution for both to later interpolate between them.
        Restrict it to the first two arcs in the list (in case the feature is
        being unknowingly used).
        '''
        # Global arc solution
        if os.path.isfile(super_arc_mef):
            wsol_in_fn = super_arc_mef
        else:
            wsol_in_fn = os.path.join(
                out_dir, "%s.p%s.fits" % (metadata["arc"][0], prev_suffix)
            )
        if not (skip_done and os.path.isfile(wsol_out_fn)
                and os.path.getmtime(wsol_in_fn) < os.path.getmtime(wsol_out_fn)):
            info_print(f"Deriving master wavelength solution from {os.path.basename(wsol_in_fn)}")
            wifes_wsol.derive_wifes_wave_solution(wsol_in_fn, wsol_out_fn, plot_dir=plot_dir_arm, **args)

        # Arc solutions for any specific obsevations
        sci_obs_list = get_sci_obs_list(metadata)
        std_obs_list = get_std_obs_list(metadata)

        for fn in sci_obs_list + std_obs_list:
            local_arcs = get_associated_calib(metadata, fn, "arc")

            if local_arcs:
                for i in range(numpy.min([2, numpy.size(local_arcs)])):
                    local_arc_fn = os.path.join(
                        out_dir, "%s.p%s.fits" % (local_arcs[i], prev_suffix)
                    )

                    local_wsol_out_fn = os.path.join(
                        out_dir, "%s.wsol.fits" % (local_arcs[i])
                    )

                    if skip_done and os.path.isfile(local_wsol_out_fn) \
                            and os.path.getmtime(local_arc_fn) < os.path.getmtime(local_wsol_out_fn):
                        continue
                    info_print(f"Deriving local wavelength solution for {local_arcs[i]}")

                    wifes_wsol.derive_wifes_wave_solution(
                        local_arc_fn,
                        local_wsol_out_fn,
                        plot_dir=plot_dir_arm,
                        **args
                    )

        return

    # ------------------------------------------------------
    # Wire solution
    # ------------------------------------------------------
    def run_wire_soln(metadata, prev_suffix, curr_suffix, **args):
        '''
        Global wire solution first, then local wire solutions for any specific observations.
        '''
        # Global wire solution
        if os.path.isfile(super_wire_mef):
            wire_in_fn = super_wire_mef
        else:
            wire_in_fn = os.path.join(
                out_dir, "%s.p%s.fits" % (metadata["wire"][0], prev_suffix)
            )
        if not (skip_done and os.path.isfile(wire_out_fn)
                and os.path.getmtime(wire_in_fn) < os.path.getmtime(wire_out_fn)):
            info_print(f"Deriving global wire solution from {os.path.basename(wire_in_fn)}")
            pywifes.derive_wifes_wire_solution(wire_in_fn,
                                               wire_out_fn,
                                               plot_dir=plot_dir_arm,
                                               **args)

        # Wire solutions for any specific obsevations
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
                    out_dir, "%s.wire.fits" % (local_wires[0])
                )
                if skip_done and os.path.isfile(local_wire_out_fn) \
                        and os.path.getmtime(local_wire_fn) < os.path.getmtime(local_wire_out_fn):
                    continue
                info_print(f"Deriving local wire solution for {local_wires[0]}")
                pywifes.derive_wifes_wire_solution(local_wire_fn,
                                                   local_wire_out_fn,
                                                   plot_dir=plot_dir_arm,
                                                   save_prefix='local_wire_fit_params',
                                                   **args)
        return

    # ------------------------------------------------------
    # Flatfield: Response
    # ------------------------------------------------------
    def run_flat_response(metadata, prev_suffix, curr_suffix, mode="all", **args):
        '''
        Generate the flatfield response function.
        '''
        # Fit the desired style of response function
        if skip_done and os.path.isfile(flat_resp_fn) \
                and os.path.getmtime(super_dflat_mef) < os.path.getmtime(flat_resp_fn):
            return
        info_print("Generating flatfield response function")
        if mode == "all" and not os.path.isfile(super_tflat_mef):
            info_print("WARNING: No twilight superflat MEF found. Falling back to dome flat only.")
            mode = "dome"

        if mode == "all":
            pywifes.wifes_2dim_response(
                super_dflat_mef,
                super_tflat_mef,
                flat_resp_fn,
                wsol_fn=wsol_out_fn,
                plot=True,
                plot_dir=plot_dir_arm,
                **args
            )
        elif mode == "dome":
            pywifes.wifes_response_poly(
                super_dflat_mef, flat_resp_fn, wsol_fn=wsol_out_fn
            )
        else:
            error_print("Requested response mode not recognised")
            raise ValueError("Requested response mode not recognised")
        return

    # ------------------------------------------------------
    # Cosmic Rays
    # ------------------------------------------------------
    def run_cosmic_rays(
        metadata,
        prev_suffix,
        curr_suffix,
        multithread=False,
        max_processes=-1,
    ):
        '''
        Clean cosmic rays on all science and standard frames.
        '''
        # now run ONLY ON SCIENCE TARGETS AND STANDARDS
        sci_obs_list = get_sci_obs_list(metadata)
        sky_obs_list = get_sky_obs_list(metadata)
        std_obs_list = get_std_obs_list(metadata)
        for fn in sci_obs_list + sky_obs_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                continue
            info_print(f"Cleaning cosmics in {os.path.basename(in_fn)}")
            in_hdr = pyfits.getheader(in_fn)
            rdnoise = 5.0 if 'RDNOISE' not in in_hdr else in_hdr['RDNOISE']
            lacos_wifes(
                in_fn,
                out_fn,
                rdnoise=rdnoise,
                wsol_fn=wsol_out_fn,
                niter=3,
                sig_clip=10.0,
                obj_lim=10.0,
                sig_frac=0.2,
                is_multithread=multithread,
                max_processes=max_processes,
            )
            if obs_mode == "ns":
                in_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, prev_suffix))
                out_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, curr_suffix))
                info_print(f"Cleaning cosmics in {os.path.basename(in_fn)}")
                in_hdr = pyfits.getheader(in_fn)
                rdnoise = 5.0 if 'RDNOISE' not in in_hdr else in_hdr['RDNOISE']
                lacos_wifes(
                    in_fn,
                    out_fn,
                    rdnoise=rdnoise,
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
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                continue
            info_print(f"Cleaning cosmics in standard star {os.path.basename(in_fn)}")
            in_hdr = pyfits.getheader(in_fn)
            rdnoise = 5.0 if 'RDNOISE' not in in_hdr else in_hdr['RDNOISE']
            lacos_wifes(
                in_fn,
                out_fn,
                rdnoise=rdnoise,
                wsol_fn=wsol_out_fn,
                niter=3,
                sig_clip=10.0,
                obj_lim=10.0,
                sig_frac=0.2,
                is_multithread=multithread,
                max_processes=max_processes,
            )
            if obs_mode == "ns":
                in_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, prev_suffix))
                out_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, curr_suffix))
                info_print(f"Cleaning cosmics in standard star {os.path.basename(in_fn)}")
                in_hdr = pyfits.getheader(in_fn)
                rdnoise = 5.0 if 'RDNOISE' not in in_hdr else in_hdr['RDNOISE']
                lacos_wifes(
                    in_fn,
                    out_fn,
                    rdnoise=rdnoise,
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
        Subtract sky frames from science objects.
        '''
        sci_obs_list = get_sci_obs_list(metadata)
        std_obs_list = get_std_obs_list(metadata)
        ns_proc_list = sci_obs_list + std_obs_list
        for fn in ns_proc_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            sky_fn = os.path.join(out_dir, "%s.s%s.fits" % (fn, prev_suffix))
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn) \
                    and os.path.getmtime(sky_fn) < os.path.getmtime(out_fn):
                continue
            info_print(f"Subtracting N+S sky frame for {os.path.basename(in_fn)}")
            pywifes.scaled_imarith_mef(in_fn, "-", sky_fn, out_fn, scale="exptime")
        return

    def run_sky_sub(metadata, prev_suffix, curr_suffix):
        if obs_mode == "ns":
            run_sky_sub_ns(metadata, prev_suffix, curr_suffix)
        else:
            # subtract sky frames from science objects
            for obs in metadata["sci"]:
                if len(obs["sky"]) > 0:
                    if len(obs["sky"]) > 1:
                        # If multiple sky frames, scale by exposure time and median-combine.
                        # List will be recombined for every science frame in case the first sky
                        # (which defines the output filename) was associated with multiple science
                        # frames. Thus, no consideration of skip_done.
                        in_fn_list = [
                            os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
                            for fn in obs["sky"]
                        ]
                        sky_proc_fn = os.path.join(
                            out_dir, "%s.p%s.fits" % (obs["sky"][0], curr_suffix)
                        )
                        info_print(f"Coadding sky frames into {os.path.basename(sky_proc_fn)}")
                        pywifes.imcombine_mef(in_fn_list, sky_proc_fn, scale="exptime", method="median")
                    else:
                        sky_fn = obs["sky"][0]
                        sky_proc_fn = os.path.join(
                            out_dir, "%s.p%s.fits" % (sky_fn, prev_suffix)
                        )
                    for fn in obs["sci"]:
                        in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
                        out_fn = os.path.join(
                            out_dir, "%s.p%s.fits" % (fn, curr_suffix)
                        )
                        if skip_done and os.path.isfile(out_fn) \
                                and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                            continue
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
                        if skip_done and os.path.isfile(out_fn) \
                                and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                            continue
                        info_print(f"Copying science image {os.path.basename(in_fn)}")
                        pywifes.imcopy(in_fn, out_fn)
            # copy stdstar frames
            std_obs_list = get_std_obs_list(metadata)
            for fn in std_obs_list:
                in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
                out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
                if skip_done and os.path.isfile(out_fn) \
                        and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
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
        '''
        for obs in metadata["sci"] + metadata["std"]:
            # If just one, then copy it
            if len(obs["sci"]) == 1:
                fn = obs["sci"][0]
                in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
                out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
                if skip_done and os.path.isfile(out_fn) \
                        and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
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
                if skip_done and os.path.isfile(out_fn) \
                        and os.path.getmtime(in_fn_list[0]) < os.path.getmtime(out_fn):
                    continue
                info_print(f"Coadding images for {os.path.basename(in_fn_list[0])}")
                pywifes.imcombine_mef(in_fn_list, out_fn, scale=scale, method=method)
        return

    # ------------------------------------------------------
    # Flatfield: Division
    # ------------------------------------------------------
    def run_flatfield(metadata, prev_suffix, curr_suffix):
        '''
        Flat-field correction to science and standard frames.
        '''
        sci_obs_list = get_primary_sci_obs_list(metadata)
        std_obs_list = get_primary_std_obs_list(metadata)
        info_print(f"Primary science observation list: {sci_obs_list}")
        info_print(f"Primary standard observation list: {std_obs_list}")
        for fn in sci_obs_list + std_obs_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn) \
                    and os.path.getmtime(flat_resp_fn) < os.path.getmtime(out_fn):
                continue
            info_print(f"Flat-fielding image {os.path.basename(in_fn)}")
            pywifes.imarith_mef(in_fn, "/", flat_resp_fn, out_fn)
            ffh = pyfits.getheader(flat_resp_fn)
            try:
                ffin = ffh["PYWRESIN"]
            except:
                ffin = "Unknown"
            of = pyfits.open(out_fn, mode="update")
            of[0].header.set("PYWRESIN", ffin, "PyWiFeS: flatfield inputs")
            of.close()
        return

    # ------------------------------------------------------
    # Data Cube Generation
    # ------------------------------------------------------
    def run_cube_gen(metadata, prev_suffix, curr_suffix, **args):
        '''
        Generate data cubes for science and standard frames.
        '''
        sci_obs_list = get_primary_sci_obs_list(metadata)
        std_obs_list = get_primary_std_obs_list(metadata)
        for fn in sci_obs_list + std_obs_list:
            in_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, prev_suffix))
            out_fn = os.path.join(out_dir, "%s.p%s.fits" % (fn, curr_suffix))
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
                        int(arc_times[0].split("-")[0]),
                        int(arc_times[0].split("-")[1]),
                        int(arc_times[0].split("-")[2].split("T")[0]),
                        int(arc_times[0].split("T")[1].split(":")[0]),
                        int(arc_times[0].split(":")[1]),
                        int(arc_times[0].split(":")[2].split(".")[0]),
                    )
                    t1 = datetime.datetime(
                        int(sci_time.split("-")[0]),
                        int(sci_time.split("-")[1]),
                        int(sci_time.split("-")[2].split("T")[0]),
                        int(sci_time.split("T")[1].split(":")[0]),
                        int(sci_time.split(":")[1]),
                        int(sci_time.split(":")[2].split(".")[0]),
                    )
                    t2 = datetime.datetime(
                        int(arc_times[1].split("-")[0]),
                        int(arc_times[1].split("-")[1]),
                        int(arc_times[1].split("-")[2].split("T")[0]),
                        int(arc_times[1].split("T")[1].split(":")[0]),
                        int(arc_times[1].split(":")[1]),
                        int(arc_times[1].split(":")[2].split(".")[0]),
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
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn) \
                    and os.path.getmtime(wire_fn) < os.path.getmtime(out_fn) \
                    and os.path.getmtime(wsol_fn) < os.path.getmtime(out_fn):
                continue

            info_print(f"Generating Data Cube for {os.path.basename(in_fn)}")
            # All done, let's generate the cube
            pywifes.generate_wifes_cube(
                in_fn,
                out_fn,
                wire_fn=wire_fn,
                wsol_fn=wsol_fn,
                **args
            )
        return

    # ------------------------------------------------------
    # Standard star extraction
    # ------------------------------------------------------
    def run_extract_stars(metadata, prev_suffix, curr_suffix, type="all", **args):
        '''
        Extract standard stars spectrum.
        '''
        # For each std, extract spectrum as desired
        std_obs_list = get_primary_std_obs_list(metadata, type=type)
        for fn in std_obs_list:
            in_fn = os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
            out_fn = os.path.join(out_dir, f"{fn}.x{prev_suffix}.dat")
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                continue
            info_print(f"Extract {type} standard star from {os.path.basename(in_fn)}")
            wifes_calib.extract_wifes_stdstar(
                in_fn, save_fn=out_fn, save_mode="ascii", **args
            )
        return

    # Sensitivity Function fit
    def run_derive_calib(metadata, prev_suffix, curr_suffix, method="poly", **args):
        '''
        Derive the sensitivity function from the extracted standard stars.
        kwargs:
            - method: 'poly' or 'spline'
            - order: polynomial order
            - knots: number of knots for the spline
        '''
        std_obs_list = get_primary_std_obs_list(metadata, type="flux")

        if len(std_obs_list) == 0:
            info_print("No flux standard stars to derive calibration. Skipping.")
            return

        std_cube_list = [
            os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
            for fn in std_obs_list
        ]
        extract_list = [
            os.path.join(out_dir, f"{fn}.x{prev_suffix}.dat")
            for fn in std_obs_list
        ]
        if skip_done and os.path.isfile(calib_fn) \
                and os.path.getmtime(std_cube_list[0]) < os.path.getmtime(calib_fn):
            return
        info_print("Deriving sensitivity function")
        wifes_calib.derive_wifes_calibration(
            std_cube_list, calib_fn, extract_in_list=extract_list, method=method, plot_dir=plot_dir_arm, **args
        )
        return

    # ------------------------------------------------------
    # Flux Calibration
    # ------------------------------------------------------
    def run_flux_calib(metadata, prev_suffix, curr_suffix, mode="pywifes", **args):
        '''
        Flux calibrate all science and standard observations.
        '''
        sci_obs_list = get_primary_sci_obs_list(metadata)
        std_obs_list = get_primary_std_obs_list(metadata)
        for fn in sci_obs_list + std_obs_list:
            in_fn = os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
            out_fn = os.path.join(out_dir, f"{fn}.p{curr_suffix}.fits")
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                continue
            info_print(f"Flux-calibrating cube {os.path.basename(in_fn)}")
            wifes_calib.calibrate_wifes_cube(in_fn, out_fn, calib_fn, mode, **args)
        return

    # ------------------------------------------------------
    # Telluric - derive
    # ------------------------------------------------------
    def run_derive_telluric(metadata, prev_suffix, curr_suffix, **args):
        '''
        Derive the telluric correction from the standard star.
        '''
        std_obs_list = get_primary_std_obs_list(metadata, "telluric")
        if len(std_obs_list) == 0:
            info_print("No telluric standard stars found. Skipping.")
            return

        std_cube_list = [
            os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
            for fn in std_obs_list
        ]
        extract_list = [
            os.path.join(out_dir, f"{fn}.x{prev_suffix}.dat")
            for fn in std_obs_list
        ]
        if skip_done and os.path.isfile(tellcorr_fn) \
                and os.path.getmtime(std_cube_list[0]) < os.path.getmtime(tellcorr_fn):
            return
        info_print("Deriving telluric correction")
        wifes_calib.derive_wifes_telluric(
            std_cube_list, tellcorr_fn, extract_in_list=extract_list, plot_dir=plot_dir_arm, **args
        )
        return

    def run_telluric_corr(metadata, prev_suffix, curr_suffix):
        '''
        Apply telluric correction for all science and standard observations.
        '''
        sci_obs_list = get_primary_sci_obs_list(metadata)
        std_obs_list = get_primary_std_obs_list(metadata)
        for fn in sci_obs_list + std_obs_list:
            in_fn = os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
            out_fn = os.path.join(out_dir, f"{fn}.p{curr_suffix}.fits")
            if skip_done and os.path.isfile(out_fn) \
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn):
                continue
            info_print(f"Correcting telluric in {os.path.basename(in_fn)}")
            wifes_calib.apply_wifes_telluric(in_fn, out_fn, tellcorr_fn)
        return

    def run_save_3dcube(metadata, prev_suffix, curr_suffix, **args):
        '''
        Save 3D Data Cube for all science and standard observations.
        '''
        sci_obs_list = get_primary_sci_obs_list(metadata)
        std_obs_list = get_primary_std_obs_list(metadata)

        # Check if is half-frame from the first sci image
        sci_filename = temp_data_dir + sci_obs_list[0] + ".fits"

        halfframe = is_halfframe(sci_filename)
        taros = is_taros(sci_filename)

        # now generate cubes
        for fn in sci_obs_list + std_obs_list:
            in_fn = os.path.join(out_dir, f"{fn}.p{prev_suffix}.fits")
            out_fn = os.path.join(out_dir, f"{fn}.{curr_suffix}.fits")
            if skip_done and (
                (
                    os.path.isfile(out_fn)
                    and os.path.getmtime(in_fn) < os.path.getmtime(out_fn)
                ) or (
                    os.path.isfile(os.path.join(working_dir, "data_products", f"{fn}.{curr_suffix}.fits"))
                    and os.path.getmtime(in_fn) < os.path.getmtime(os.path.join(working_dir, "data_products", f"{fn}.{curr_suffix}.fits"))
                )
            ):
                continue
            info_print(f"Saving 3D Data Cube for {os.path.basename(in_fn)}")
            pywifes.generate_wifes_3dcube(in_fn, out_fn, halfframe=halfframe, taros=taros, **args)
        return

    # --------------------------------------------
    # INITIATE THE SCRIPT
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

    # Option for specifying to auto-extract and splice the output datacubes
    parser.add_argument(
        "-extract-and-splice",
        action="store_true",
        help="Optional: Auto-extract and splice the datacubes.",
    )

    # Option for triggering the reduction from master calibrations
    parser.add_argument(
        "--from-master",
        type=str,
        const="./data_products/master_calib",
        nargs="?",
        help="Optional: Path to the master calibrations directory. If not provided, the default path will be used: .",
    )

    # Option for specifying to only produce the master calibration files
    parser.add_argument(
        "-just-calib",
        action="store_true",
        help="Optional: Only basics master calibration files will produced.",
    )

    # Option for specifying to skip already completed steps
    parser.add_argument(
        "-skip-done",
        action="store_true",
        help="Optional: Skip already completed steps.",
    )

    # Option for treating OBJECT images near standard stars as standards even if IMAGETYP != STANDARD
    parser.add_argument(
        "-greedy-stds",
        action="store_true",
        help="Optional: Treat OBJECT as STANDARD if near known standard.",
    )

    # Option to avoid coadding all science frames with the same OBJECT name.
    # Useful for time-series, position shifts, etc.
    parser.add_argument(
        "--coadd-mode",
        choices=["all", "none", "prompt"],
        default="all",
        help="Optional: coadd 'all' (default), 'none', 'prompt' (user-selected, blue/red separate) frames of a given OBJECT",
        required=False,
        type=str,
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

    # Only basics master calibration
    just_calib = args.just_calib

    # Set to skip already done files
    skip_done = args.skip_done

    # Auto-extract and splice the datacubes
    extract_and_splice = args.extract_and_splice

    # Creates a temporary data directory containning all raw data for reduction.
    temp_data_dir = os.path.join(working_dir, "data_products/intermediate/raw_data_temp/")
    os.makedirs(temp_data_dir, exist_ok=True)

    all_fits_names = get_file_names(user_data_dir, "*.fits*")
    # Copy raw data from user's direcory into temporaty raw directory.
    copy_files(user_data_dir, temp_data_dir, all_fits_names)

    # Creates a directory for plot.
    plot_dir = os.path.join(working_dir, "data_products/plots/")
    os.makedirs(plot_dir, exist_ok=True)

    # Classify all raw data (red and blue arm)
    mode_save_fn = os.path.join(working_dir, "data_products/coadd_mode.json5")
    obs_metadatas = classify(temp_data_dir,
                             greedy_stds=args.greedy_stds,
                             coadd_mode=args.coadd_mode,
                             mode_save_fn=mode_save_fn)

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

    info_print(f"Processing using master calibrations {'to' if just_calib else 'from'}: '{master_dir}'.")

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
            obs_mode = "classical"
            if is_nodshuffle(temp_data_dir + reference_filename) \
                    or is_subnodshuffle(temp_data_dir + reference_filename):
                obs_mode = "ns"

            # Check if is half-frame
            halfframe = is_halfframe(temp_data_dir + reference_filename)
            if halfframe:
                obs_metadata = pywifes.calib_to_half_frame(obs_metadata, temp_data_dir)

            # Grism
            grism = pyfits.getheader(temp_data_dir + reference_filename)[grism_key[arm]]

            # Set the JSON file path and read it.
            if params_path[arm] is None:
                json_path = f"./pipeline_params/{arm}/params_{grism}.json5"
            else:
                json_path = params_path[arm]

            # Load the JSON file
            proc_steps = load_config_file(json_path)

            # Remove irrelevant steps for just_calib usage
            if just_calib:
                json_idx = []
                for this_snum, this_step in enumerate(proc_steps[arm]):
                    if this_step['step'] not in ['sky_sub', 'obs_coadd', 'flatfield',
                                                 'cube_gen', 'extract_stars', 'derive_calib',
                                                 'flux_calib', 'derive_telluric', 'telluric_corr',
                                                 'save_3dcube']:
                        json_idx.append(this_snum)
                proc_steps[arm] = [proc_steps[arm][i] for i in json_idx]

            # Create data products directory structure
            out_dir = os.path.join(working_dir, f"data_products/intermediate/{arm}")
            os.makedirs(out_dir, exist_ok=True)

            calib_prefix = f"wifes_{arm}"

            # WiFeS specific parameter
            my_data_hdu = 0

            # Creates a directory for diagnisis plots (one per arm).
            plot_dir_arm = os.path.join(plot_dir, arm)
            os.makedirs(plot_dir_arm, exist_ok=True)

            # ------------------------------------------------------------------------
            # Define names for master calibration files and set their path.
            # ------------------------------------------------------------------------
            # Bias Master Files
            overscanmask_fn = os.path.join(master_dir, "%s_overscanmask.fits" % calib_prefix)
            superbias_fn = os.path.join(master_dir, "%s_superbias.fits" % calib_prefix)
            superbias_fit_fn = os.path.join(master_dir, "%s_superbias_fit.fits" % calib_prefix)

            # Dome Master Files
            super_dflat_raw = os.path.join(master_dir, "%s_super_domeflat_raw.fits" % calib_prefix)
            super_dflat_fn = os.path.join(master_dir, "%s_super_domeflat.fits" % calib_prefix)
            super_dflat_mef = os.path.join(master_dir, "%s_super_domeflat_mef.fits" % calib_prefix)

            # Twilight Master Files
            super_tflat_raw = os.path.join(master_dir, "%s_super_twiflat_raw.fits" % calib_prefix)
            super_tflat_fn = os.path.join(master_dir, "%s_super_twiflat.fits" % calib_prefix)
            super_tflat_mef = os.path.join(master_dir, "%s_super_twiflat_mef.fits" % calib_prefix)

            # Wire Master Files
            super_wire_raw = os.path.join(master_dir, "%s_super_wire_raw.fits" % calib_prefix)
            super_wire_mef = os.path.join(master_dir, "%s_super_wire_mef.fits" % calib_prefix)

            # Arc Master Files
            super_arc_raw = os.path.join(master_dir, "%s_super_arc_raw.fits" % calib_prefix)
            super_arc_mef = os.path.join(master_dir, "%s_super_arc_mef.fits" % calib_prefix)

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
            info_print("________________________________________________________________")
            info_print(f"Starting processing of {arm} arm")
            info_print("________________________________________________________________")

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
                    if step_suffix is not None:
                        prev_suffix = step_suffix

                else:
                    pass

            info_print(f"Successfully completed arm {arm}\n")

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

        if extract_and_splice:
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
            extract_params = load_config_file("./pipeline_params/params_extract.json5")

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
                plot_path = os.path.join(plot_dir, plot_name)
                plot = extract_params["plot"]
                taros = is_taros(blue_cube_path)
                # Run auto-extraction
                detect_extract_and_save(
                    blue_cube_path,
                    red_cube_path,
                    destination_dir,
                    r_arcsec=extract_params["r_arcsec"],
                    border_width=extract_params["border_width"],
                    sky_sub=False if obs_mode == "ns" else True,
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
                    if taros:
                        spliced_cube_name = blue_cube_name.replace("T2m3wb", "T2m3wSplice")
                    else:
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
                        if taros:
                            spliced_spectrum_name = os.path.basename(blue_spec).replace(
                                "T2m3wb", "T2m3wSplice"
                            )
                        else:
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
                            if taros:
                                spliced_spec = blue_spec.replace("T2m3wb", "T2m3wSplice")
                            else:
                                spliced_spec = blue_spec.replace("Blue", "Splice")

                            plot_1D_spectrum(blue_spec, plot_dir=plot_dir)
                            plot_1D_spectrum(red_spec, plot_dir=plot_dir)
                            plot_1D_spectrum(spliced_spec, plot_dir=plot_dir)

                    elif blue_cube_path is not None:

                        blue_pattern = blue_cube_path.replace("cube.fits", "spec.ap*")
                        blue_specs = glob.glob(blue_pattern)

                        for blue_spec in blue_specs:

                            plot_1D_spectrum(blue_spec, plot_dir=plot_dir)

                    elif red_cube_path is not None:

                        red_pattern = red_cube_path.replace("cube.fits", "spec.ap*")
                        red_specs = glob.glob(red_pattern)

                        for red_spec in red_specs:

                            plot_1D_spectrum(red_spec, plot_dir=plot_dir)

    # ----------------------------------------------------------
    # Print total running time
    # ----------------------------------------------------------
    end_time = datetime.datetime.now()
    duration = end_time - start_time
    messagge = "All done in %.01f seconds." % duration.total_seconds()
    info_print(messagge)
    print('\U0001F52D', messagge, '\u2B50')


if __name__ == "__main__":
    main()
