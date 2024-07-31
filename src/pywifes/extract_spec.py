import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from photutils.detection import find_peaks
from photutils.centroids import centroid_com
from photutils.aperture import EllipticalAperture
from photutils.aperture import EllipticalAnnulus
import matplotlib.colors as mcolors
from matplotlib.patheffects import withStroke
from astropy.wcs import WCS
import os
import re

# Suppress the NoDetectionsWarning as we have set a warning for no detection
import warnings
from photutils.utils.exceptions import NoDetectionsWarning

warnings.filterwarnings("ignore", category=NoDetectionsWarning)

# Logger
from pywifes.logger_config import custom_print
import logging


# Redirect print statements to logger
logger = logging.getLogger("PyWiFeS")
print = custom_print(logger)


def extract_aperture_name(spec_name):
    match = re.search(r"ap(\d)", spec_name)
    if match:
        aperture_number = match.group(1)
        return f"Aperture {aperture_number}"
    return None


def extract_and_save(
    cube_path, sci, var, source_ap, ap_index, sky_ap, output_dir, sci_hdr, var_hdr
):
    if cube_path is not None:
        # Extraction
        flux, var = spect_extract(sci, var, source_ap, sky_ap=sky_ap)

        # Write out the results
        base = os.path.basename(cube_path)
        base_output = base.replace("cube.fits", f"spec.ap{ap_index}.fits")
        print(f"Saving extracted spectra for {base}")
        output_path = os.path.join(output_dir, base_output)
        write_1D_spec(flux, var, sci_hdr, var_hdr, output_path)


def plot_arm(ax, cube_path, sci, title, source_apertures, sky_aps, border_width):
    ax.set_title(title)
    if cube_path is not None:
        plot_apertures(
            sci, source_apertures, sky_aps=sky_aps, border_width=border_width
        )
        ax.set_xlabel("Right Ascension")
        ax.set_ylabel("Declination ")
    else:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        text = f"No {title.lower()} arm data"
        ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=12, color="black")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")


def collapse_cube(*data_cubes):
    """
    Collapse one or more data cubes by computing the median along the wavelength dimension.

    Parameters
    ----------
    *data_cubes : array_like
        Variable-length argument list of data cubes to be collapsed.
        Each data cube should have the same spatial shape.

    Returns
    -------
    median_image : ndarray
        2D array representing the median-collapsed image along the wavelength dimension.

    Notes
    -----
    This function takes one or more data cubes and collapses them by computing the median
    along the wavelength dimension. The resulting 2D array represents the
    median image.
    """
    # Filter out None elements from the data_cubes
    valid_data_cubes = [cube for cube in data_cubes if cube is not None]

    if len(valid_data_cubes) == 1:
        joined_cubes = valid_data_cubes[0]
    else:
        joined_cubes = np.concatenate(valid_data_cubes, axis=0)

    median = np.nanmedian(joined_cubes, axis=0)

    return median


def sec_image(ap, image):
    """
    Extract pixel values from a section within the aperture.

    Parameters
    ----------
    aperture : `photutils.aperture` object
        Aperture object representing the region of interest.

    image : array_like
        2D array representing the image data.

    Returns
    -------
    sec_data : ndarray
        1D array containing the weighted pixel values from the section within the aperture.

    Notes
    -----
    This function extracts pixel values from a section within the specified aperture
    from the input image. The extraction is performed using a aperture mask,
    where pixel values outside the aperture are zero. The resulting
    1D array contains the pixel values from the section within the aperture.
    """
    mask = ap.to_mask(method="center")
    sec_data = mask.get_values(image)

    return sec_data


def spect_extract(sci_cube, var_cube, source_ap, sky_ap=None):
    """
    Extract and calculate spectral flux and variance within the given aperture for each layer or wavelenght step.

    Parameters
    ----------
    sci_cube : array_like
        3D array representing the scientific data cube.

    var_cube : array_like
        3D array representing the variance data cube.

    source_ap : `photutils.aperture` object
        Aperture object representing the source region to be extracted.

    sky_ap : `photutils.aperture` object, optional
        Aperture object representing the sky region to be extracted.
        Default is None.

    Returns
    -------
    flux : ndarray
        1D array containing the spectral flux values.

    variance : ndarray
        1D array containing the spectral variance values.

    Notes
    -----
    This function extracts spectral flux and variance values within the specified
    source aperture from the input data cubes. If a sky aperture is provided,
    it calculates and subtracts the corresponding sky values.

    """
    # Extract the pixels inside the aperture for all the wavelentghs
    fl = []
    var = []

    for layer in range(sci_cube.shape[0]):
        sci_section = sec_image(source_ap, sci_cube[layer])
        var_section = sec_image(source_ap, var_cube[layer])
        area = source_ap.area_overlap(sci_cube[layer], method="center")

        average = np.average(
            sci_section,
            weights=np.reciprocal(var_section),
        )

        error = np.reciprocal(np.sum(np.reciprocal(var_section)))

        if sky_ap is not None:
            sky_section = sec_image(sky_ap, sci_cube[layer])
            sky_var_section = sec_image(sky_ap, var_cube[layer])

            sky_average = np.average(
                sky_section,
                weights=np.reciprocal(sky_var_section),
            )

            sky_error = np.reciprocal(np.sum(np.reciprocal(sky_var_section)))

            fl.append((average - sky_average) * area)
            var.append((error + sky_error) * area)

        else:
            fl.append(average * area)
            var.append(error * area)

    return np.array(fl), np.array(var)


def write_1D_spec(sci_data, var_data, sci_cube_header, var_cube_header, output):
    """
    This function writes the 1D spectrum data and its variance to a FITS file.
    It updates the WCS headers to match the 1D spectrum header.

    Parameters
    ----------
    sci_data : numpy.ndarray
        1D array containing the science data (flux) of the spectrum.
    var_data : numpy.ndarray
        1D array containing the variance data of the spectrum.
    sci_cube_header : astropy.io.fits.Header
        Header from the science cube.
    var_cube_header : astropy.io.fits.Header
        Header object from the variance data.
    output : str
        Output file path for the FITS file containing the 1D spectrum.

    Returns
    -------
    None
        This function does not return any value. It writes the corresponding fits.

    """

    # Create a blank HDUList for the 1D spectrum
    hdulist = fits.HDUList(fits.PrimaryHDU())

    # Update WCS headers to match the 1D spectrum header
    sci_cube_header_copy = sci_cube_header.copy()
    var_cube_header_copy = var_cube_header.copy()
    headers = [sci_cube_header_copy, var_cube_header_copy]

    for header in headers:
        # Update axis 1 WCS information to match wavelength solution (axis 3)
        header["CDELT1"] = header["CDELT3"]
        header["CRPIX1"] = header["CRPIX3"]
        header["CRVAL1"] = header["CRVAL3"]
        header["CUNIT1"] = "Angstrom"
        header["CTYPE1"] = header["CTYPE3"]
        header["NAXIS"] = 1
        header["NAXIS1"] = header["NAXIS3"]

        # Remove spatial WCS keys
        keys_to_remove = [
            "NAXIS2",
            "NAXIS3",
            "CRVAL2",
            "CRVAL3",
            "CTYPE2",
            "CTYPE3",
            "CDELT2",
            "CDELT3",
            "CRPIX2",
            "CRPIX3",
            "PC1_1",
            "PC1_2",
            "PC2_1",
            "PC2_2",
        ]
        for key in keys_to_remove:
            header.remove(key)

    # Set science data and headers in the blank HDUList
    hdulist[0].data = sci_data
    hdulist[0].header = sci_cube_header_copy

    # Add variance data as an additional HDU
    hdu_fluxvar = fits.ImageHDU(data=var_data, header=var_cube_header_copy)
    hdulist.append(hdu_fluxvar)

    # Write the FITS file containing the 1D spectrum
    hdulist.writeto(output, overwrite=True)
    hdulist.close()


def plot_apertures(data_cube, source_apertures, sky_aps=None, border_width=0):
    """
    Plot apertures on a collapsed data cube image.

    Parameters
    ----------
    data_cube : numpy.ndarray
        The 3D data cube containing the image data.

    source_apertures : list of photutils.aperture objects
        List of source apertures to plot on the image.

    sky_aps : list of photutils.aperture objects, optional
        List of sky apertures to plot on the image. Default is None.

    border_width : int, optional
        Width of the border to exclude from the image when calculating the image contrast.
        Default is 0.

    Returns
    -------
    None
        This function does not return any value. It plots the apertures on the image.

    """

    # Collapse cube in the waevelenght dimesion for obtaing a median image
    collapsed_cube = collapse_cube(data_cube)

    # Improbe ccontrats on the image usefull for see faint sources
    vmin, vmax = np.percentile(
        collapsed_cube[border_width:-border_width, border_width:-border_width], (5, 95)
    )
    # Show collapsed image
    plt.imshow(collapsed_cube, vmin=vmin, vmax=vmax)

    # Color of the overlay area
    cmap = mcolors.ListedColormap(["white"])
    alpha = 0.5

    for index, source_ap in enumerate(source_apertures):
        ap_index = index + 1
        # Plot a overlaped transparent area
        mask = source_ap.to_mask(method="center").to_image(np.shape(collapsed_cube))
        mask[mask == 0] = np.nan
        plt.imshow(mask, alpha=alpha, cmap=cmap)
        # Plot theoretical aperture contourngit
        source_ap.plot(color="white", lw=0.8, ls="--")
        # Plot the aperture number
        plt.text(
            source_ap.positions[0],
            source_ap.positions[1],
            ap_index,
            color="white",
            ha="center",
            va="center",
            fontsize=12,
            path_effects=[withStroke(linewidth=2, foreground="black")],
        )

    if sky_aps is not None:
        for sky_ap in sky_aps:
            # Plot a overlaped transparent area
            mask = sky_ap.to_mask(method="center").to_image(np.shape(collapsed_cube))
            mask[mask == 0] = np.nan
            plt.imshow(mask, alpha=alpha, cmap=cmap)
            # Plot theoretical aperture contourn
            sky_ap.plot(color="white", lw=0.8, ls="--")


def read_cube_data(cube_path):
    cube_data = {
        "sci": None,
        "sci_hdr": None,
        "var": None,
        "var_hdr": None,
        "wcs": None,
        "binning_x": None,
        "binning_y": None,
    }

    if cube_path is not None:
        sci, sci_hdr = fits.getdata(cube_path, 0, header=True)
        var, var_hdr = fits.getdata(cube_path, 1, header=True)
        wcs = WCS(sci_hdr).celestial
        binning_x, binning_y = np.int_(sci_hdr["CCDSUM"].split())
        cube_data["sci"] = sci
        cube_data["sci_hdr"] = sci_hdr
        cube_data["var"] = var
        cube_data["var_hdr"] = var_hdr
        cube_data["wcs"] = wcs
        cube_data["binning_x"] = binning_x
        cube_data["binning_y"] = binning_y

    return cube_data


def detect_extract_and_save(
    blue_cube_path=None,
    red_cube_path=None,
    output_dir=None,
    r_arcsec=2,
    border_width=3,
    sky_sub=False,
    plot=False,
    plot_path="detected_apertures_plot.png",
):
    """
    Detects sources in the input cubes, extracts and saves their spectra. Optionally, it can plot the extracted spurces and the detected apertures around them.

    Parameters:
    -----------
    blue_cube_path : str, optional
        Path to the blue cube file.
    red_cube_path : str, optional
        Path to the red cube file.
    output_dir : str, optional
        Directory where the extracted spectra will be saved.
    r_arcsec : float, optional
        Radius of the circular aperture in arcseconds.
    border_width : int, optional
        Width of the border to be excluded from the statistics.
    sky_sub : bool, optional
        Flag indicating whether to perform sky subtraction.
    plot : bool, optional
        Flag indicating whether to generate a plot of the detected apertures.
    plot_path : str, optional
        Path to save the plot.

    Returns:
    --------
    None
    """
    # reading the data from cubes
    # Blue arm
    blue_cube_data = read_cube_data(blue_cube_path)
    blue_sci = blue_cube_data["sci"]
    blue_sci_hdr = blue_cube_data["sci_hdr"]
    blue_var = blue_cube_data["var"]
    blue_var_hdr = blue_cube_data["var_hdr"]
    blue_wcs = blue_cube_data["wcs"]
    if blue_cube_data["sci"] is not None:
        binning_x = blue_cube_data["binning_x"]
        binning_y = blue_cube_data["binning_y"]
        object = blue_sci_hdr["OBJECT"]

    # Red arm
    red_cube_data = read_cube_data(red_cube_path)
    red_sci = red_cube_data["sci"]
    red_sci_hdr = red_cube_data["sci_hdr"]
    red_var = red_cube_data["var"]
    red_var_hdr = red_cube_data["var_hdr"]
    red_wcs = red_cube_data["wcs"]
    if red_cube_data["sci"] is not None:
        binning_x = red_cube_data["binning_x"]
        binning_y = red_cube_data["binning_y"]
        object = red_sci_hdr["OBJECT"]

    # Calculate pixel scale from binning
    pixel_scale_x = binning_x  # arcsec/pix
    pixel_scale_y = binning_y / 2  # arcsec/pix

    # Average all the data for detecting the source
    collapsed_cube = collapse_cube(blue_sci, red_sci)

    # Automatic source detection in the collapsed cubes (red + blue)
    npeaks = 3  # Number of peaks to be detected
    # Avoid the edges of size border_width on the statistics
    collapsed_cube_no_edge = collapsed_cube[
        border_width:-border_width, border_width:-border_width
    ]
    mean, median, std = sigma_clipped_stats(collapsed_cube_no_edge, sigma=5.0)
    threshold = median + (1.5 * std)

    detection = find_peaks(
        collapsed_cube,
        threshold,
        border_width=border_width,
        npeaks=npeaks,
        centroid_func=centroid_com,
    )

    if detection is None:
        print("No source detected.")

    else:
        # Sort detections as brighter peaks first
        detection.sort("peak_value", reverse=True)

        # Detected sources positions
        positions = np.transpose((detection["x_peak"], detection["y_peak"]))

        # Set the annulus
        a = r_arcsec / pixel_scale_x
        b = r_arcsec / pixel_scale_y

        # Creates source apertures in the detected positions
        source_apertures = EllipticalAperture(positions, a=a, b=b)

        if sky_sub:
            a_in = a * 3
            a_out = a * 4

            b_in = b * 3
            b_out = b * 4

            sky_aps = EllipticalAnnulus(
                positions,
                a_in=a_in,
                a_out=a_out,
                b_in=b_in,
                b_out=b_out,
            )

        else:
            sky_aps = None

        # Flux extraction for the aperture at all the wavelenghts
        for index, source_ap in enumerate(source_apertures):
            if sky_aps is not None:
                sky_ap = sky_aps[index]

            else:
                sky_ap = None
            ap_index = index + 1

            # Extracting blue cube
            extract_and_save(
                blue_cube_path,
                blue_sci,
                blue_var,
                source_ap,
                ap_index,
                sky_ap,
                output_dir,
                blue_sci_hdr,
                blue_var_hdr,
            )

            # Extracting red cube
            extract_and_save(
                red_cube_path,
                red_sci,
                red_var,
                source_ap,
                ap_index,
                sky_ap,
                output_dir,
                red_sci_hdr,
                red_var_hdr,
            )

        if plot:
            plt.suptitle(object)
            # Plot Red
            ax1 = plt.subplot(1, 2, 2, projection=red_wcs)
            plot_arm(
                ax1,
                red_cube_path,
                red_sci,
                "Red arm",
                source_apertures,
                sky_aps,
                border_width,
            )

            # Plot Blue
            ax0 = plt.subplot(1, 2, 1, projection=blue_wcs)
            plot_arm(
                ax0,
                blue_cube_path,
                blue_sci,
                "Blue arm",
                source_apertures,
                sky_aps,
                border_width,
            )

            plt.tight_layout()
            fig_output = os.path.join(output_dir, plot_path)
            plt.savefig(fig_output, bbox_inches="tight", dpi=300)
            plt.close()


# TODO Replace with a proper SpecUtils loader


class SingleSpec(object):
    """
    Class representing a single spectrum for analysis
    """

    def __init__(self, fitsFILE):
        self.flux, self.header = fits.getdata(fitsFILE, 0, header=True)
        self.wl = (
            np.arange(self.header["NAXIS1"]) - self.header["CRPIX1"] + 1
        ) * self.header["CDELT1"] + self.header["CRVAL1"]
        # Temporary
        self.fluxvar = fits.getdata(fitsFILE, 1, header=False)
        self.minWL = np.min(self.wl)
        self.maxWL = np.max(self.wl)
        return


def plot_1D_spectrum(spec_path, plot_dir):
    """
    Plot a 1D spectrum with error bars and save the plot as a PNG file in the provided directory.

    Parameters
    ----------
    spec_path : str
        The path to the spectrum file.
    plot_dir : str
        The directory where the plot will be saved.

    Returns
    -------
    None
        This function does not return anything.

    """
    spec = SingleSpec(spec_path)
    flux = spec.flux
    wl = spec.wl
    fluxvar = spec.fluxvar
    # Calculate the error as the square root of the flux variance
    flux_error = np.sqrt(fluxvar)

    # Plot the spectrum with error bars
    spec_name = os.path.basename(spec_path)
    plot_name = spec_name.replace(".fits", ".png")
    plot_path = os.path.join(plot_dir, plot_name)
    aperture_name = extract_aperture_name(spec_name)

    fig = plt.figure(figsize=(12, 5))

    # Plot the error bars
    plt.errorbar(
        wl,
        flux,
        yerr=flux_error,
        fmt="none",
        ecolor="grey",
        elinewidth=0.5,
        capsize=1.5,
    )

    # Plot the step plot for the spectrum values
    plt.step(wl, flux, where="mid", color="b")

    # Customize the plot
    plt.title(f"{spec_name} \n" + aperture_name)
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Flux")
    plt.grid(True)

    plt.savefig(plot_path, dpi=300)
    plt.close()
