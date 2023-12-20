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

    if len(valid_data_cubes)==1:
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


def spect_extract(sci_cube, var_cube, source_ap, sky_ap):
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

    #  Extract the pixels inside the apertur for all the wavelengh
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
            
        err = np.reciprocal(np.sum(np.reciprocal(var_section)))

        if sky_ap is not None:
            sky_section = sec_image(sky_ap, sci_cube[layer])
            sky_var_section = sec_image(sky_ap, var_cube[layer])

            sky_average = np.average(
                sky_section,
                weights=np.reciprocal(sky_var_section),
                ) 

            sky_err = np.reciprocal(np.sum(np.reciprocal(sky_var_section)))

            fl.append((average - sky_average) * area)
            var.append((err + sky_err) * area)  

        else:
            fl.append(average * area)
            var.append(err * area)  

    return np.array(fl), np.array(var)


def writeFITS(sci_cube, var_cube, sci_header, var_header, output):
    hdulist = fits.HDUList(fits.PrimaryHDU())
    hdulist[0].data = sci_cube
    hdulist[0].header = sci_header
    hdulist[0].header["CRPIX1"] = 1.0
    hdulist[0].header["CRVAL1"] = hdulist[0].header["CRVAL3"]
    hdulist[0].header["CDELT1"] = hdulist[0].header["CDELT3"]

    hdr_fluxvar = fits.Header()
    hdr_fluxvar = var_header
    hdr_fluxvar["CRPIX1"] = 1.0
    hdr_fluxvar["CRVAL1"] = hdr_fluxvar["CRVAL3"]
    hdr_fluxvar["CDELT1"] = hdr_fluxvar["CDELT3"]

    hdu_fluxvar = fits.ImageHDU(data=var_cube, header=var_header)
    hdulist.append(hdu_fluxvar)

    hdulist.writeto(output, overwrite=True)
    hdulist.close()

    return


def plot_apertures(data_cube, source_aps, sky_aps=None, border_width=0):
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

    for index, source_ap in enumerate(source_aps):
        ap_index = index + 1
        # Plot a overlaped transparent area
        mask = source_ap.to_mask(method="center").to_image(np.shape(collapsed_cube))
        mask[mask == 0] = np.nan
        plt.imshow(mask, alpha=alpha, cmap=cmap)
        # Plot theoretical aperture contourn
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



def auto_extract(
    blue_cube_path=None, red_cube_path=None, a=1, b=1, border_width = 3, sky_sub=False, check_plot=False
):
    # Load in the data

    # Blue arm
    if blue_cube_path is not None:
        blue_sci, b_sci_hdr = fits.getdata(blue_cube_path, 0, header=True)
        b_var, b_var_hdr = fits.getdata(blue_cube_path, 1, header=True)
        blue_wcs = WCS(b_sci_hdr).celestial

    else:
        blue_sci = None
        blue_wcs = None

    # Red arm
    if red_cube_path is not None:
        red_sci, r_sci_hdr = fits.getdata(red_cube_path, 0, header=True)
        r_var, r_var_hdr = fits.getdata(red_cube_path, 1, header=True)
        red_wcs = WCS(r_sci_hdr).celestial
    else:
        red_sci = None
        red_wcs = None


    # Average all the data for detecting the source
    collapsed_cube = collapse_cube(blue_sci, red_sci)

    # Automatic source detection in the collapsed cubes (red + blue)
    npeaks = 3  # Number of peaks to be detected
    mean, median, std = sigma_clipped_stats(collapsed_cube, sigma=5.0)
    threshold = median + (3 * std)

    detection = find_peaks(
        collapsed_cube,
        threshold,
        border_width=border_width,
        npeaks=npeaks,
        centroid_func=centroid_com,
    )

    if detection is None:
        raise ValueError("No source detected.")
    
    # Sort detections as brighter peaks first
    detection.sort('peak_value', reverse=True)
    

    # Detected sources positions
    positions = np.transpose((detection["x_peak"], detection["y_peak"]))

    # Creates source apertures in the detected positions
    source_aps = EllipticalAperture(positions, a=a, b=b)

    # Set the annulus
    if sky_sub:

        a_in = a * 3
        a_out = a * 4

        b_in = b * 3
        b_out = b * 4

        sky_aps = EllipticalAnnulus(positions, a_in=a_in, a_out=a_out, b_out=b_out, b_in=b_in)

    else:
        sky_aps = None


    # Flux extraction for the aperture at all the wavelenghts
    for index, source_ap in enumerate(source_aps):
        if sky_aps is not None:
            sky_ap = sky_aps[index]

        else:
            sky_ap = None
        ap_index = index + 1

        if blue_cube_path is not None:
            # Extraction
            blue_flux, blue_var = spect_extract(blue_sci, b_var, source_ap, sky_ap=sky_ap)

            # Write out the results
            blue_output = blue_cube_path.replace("p11.fits", "ap%s.p12.fits" % ap_index)
            writeFITS(blue_flux, blue_var, b_sci_hdr, b_var_hdr, blue_output)

        if red_cube_path is not None:
            # Extraction
            red_flux, red_var = spect_extract(red_sci, r_var, source_ap, sky_ap=sky_ap)
            # Write out the results
            red_output = red_cube_path.replace("p11.fits", "ap%s.p12.fits" % ap_index)
            writeFITS(red_flux, red_var, r_sci_hdr, r_var_hdr, red_output)

    if check_plot:

        plt.close("all")

        # Plot Red
        ax1 = plt.subplot(1, 2, 2, projection=red_wcs)
        plt.title("Red arm")
        if red_cube_path is not None:
            plot_apertures(red_sci, source_aps, sky_aps=sky_aps, border_width=border_width)
                    #  Axis labels    
            plt.xlabel('Right Ascension')
            plt.ylabel('Declination ')

        else:
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            # Add text to the center of the figure
            text = "No red arm data"
            plt.text(0.5, 0.5, text, ha='center', va='center', fontsize=12, color='black')

            # Hide ticks and labels on both axes
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.axis('off')


        # Plot Blue
        ax0 = plt.subplot(1, 2, 1, projection=blue_wcs)
        plt.title("Blue arm")
        if blue_cube_path is not None:
            plot_apertures(blue_sci, source_aps, sky_aps=sky_aps, border_width=border_width)
            plt.xlabel('Right Ascension')
            plt.ylabel('Declination ')
    
        else:
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            # Add text to the center of the figure
            text = "No blue arm data"
            plt.text(0.5, 0.5, text, ha='center', va='center', fontsize=12, color='black')

            # Hide ticks and labels on both axes
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.axis('off')


        plt.tight_layout()
        plt.savefig("detected_apertures_plot.pdf", bbox_inches="tight", dpi=300)


