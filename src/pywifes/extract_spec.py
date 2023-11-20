import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from photutils.detection import find_peaks
from photutils.centroids import centroid_com
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
import matplotlib.colors as mcolors
from matplotlib.patheffects import withStroke


def collapse_cube(*data_cubes):
    joined_cubes = np.concatenate(data_cubes, axis=0)
    median = np.nanmedian(joined_cubes, axis=0)

    return median


def sec_image(ap, image):
    mask = ap.to_mask(method="center")
    mask_data = mask.data
    sec = mask.cutout(image)
    sec_weight = sec * mask_data
    sec_data = sec_weight[mask_data > 0]

    return sec_data


def plot_apertures(red_cube, blue_cube, source_aps, sky_aps=None):
    blue_collapse = collapse_cube(blue_cube)
    ave_collapse = collapse_cube(red_cube)

    plt.close("all")

    # plt.figure(1)
    # Plot Red
    plt.subplot(1, 2, 2)
    plt.title("Red arm")
    vmin, vmax = np.percentile(ave_collapse, (5, 95))
    plt.imshow(ave_collapse, vmin=vmin, vmax=vmax)

    # Color of the overlay area
    cmap = mcolors.ListedColormap(["white"])
    alpha = 0.5

    for index, source_ap in enumerate(source_aps):
        ap_index = index + 1
        # Plot a overlaped transparent area
        mask = source_ap.to_mask(method="center").to_image(np.shape(ave_collapse))
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
            mask = sky_ap.to_mask(method="center").to_image(np.shape(ave_collapse))
            mask[mask == 0] = np.nan
            plt.imshow(mask, alpha=alpha, cmap=cmap)
            # Plot theoretical aperture contourn
            sky_ap.plot(color="white", lw=0.8, ls="--")

    # Plot blue
    plt.subplot(1, 2, 1)
    plt.title("Blue arm")
    vmin, vmax = np.percentile(blue_collapse, (5, 95))
    plt.imshow(blue_collapse, vmin=vmin, vmax=vmax)
    for source_ap in source_aps:
        # Plot a overlaped transparent area
        mask = source_ap.to_mask(method="center").to_image(np.shape(blue_collapse))
        mask[mask == 0] = np.nan
        plt.imshow(mask, alpha=alpha, cmap=cmap)
        # Plot theoretical aperture contourn
        source_ap.plot(color="white", lw=0.8, ls="--")

    if sky_aps is not None:
        for sky_ap in sky_aps:
            # Plot a overlaped transparent area
            mask = sky_ap.to_mask(method="center").to_image(np.shape(blue_collapse))
            mask[mask == 0] = np.nan
            plt.imshow(mask, alpha=alpha, cmap=cmap)
            # Plot theoretical aperture contourn
            sky_ap.plot(color="white", lw=0.8, ls="--")
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

    plt.tight_layout()
    plt.savefig("check_apertures_plot.pdf", bbox_inches="tight", dpi=300)


def spect_extract(sci_cube, var_cube, source_ap, sky_ap):
    #  Extract the pixels inside the apertur for all the wavelengh
    fl = []
    var = []

    for layer in range(sci_cube.shape[0]):
        sci_section = sec_image(source_ap, sci_cube[layer])
        var_section = sec_image(source_ap, var_cube[layer])

        area = source_ap.area_overlap(sci_cube[0], method="center")

        average = np.average(
            sci_section
        )  # , weights=np.reciprocal(var_section)) TODO: deal with nan var to calculate recirpocal
        err = np.reciprocal(np.sum(np.reciprocal(var_section)))

        if sky_ap is not None:
            sky_section = sec_image(sky_ap, sci_cube[layer])
            sky_var_section = sec_image(sky_ap, var_cube[layer])
            sky_average = np.average(
                sky_section
            )  # , weights=np.reciprocal(sky_var_section)) TODO: deal with nan var to calculate recirpocal

            sky_err = np.reciprocal(np.sum(np.reciprocal(sky_var_section)))

            fl.append((average - sky_average) * area)
            var.append((err + sky_err) * area)  # TODO values of teh order of 10x-31 ???

        else:
            fl.append(average * area)
            var.append(err * area)  # TODO values of teh order of 10x-31 ???

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


def auto_extract(blue_cube_path, red_cube_path, sky_sub=False, check_plot=False):
    # Load in the data

    # Blue arm
    blue_sci, b_sci_hdr = fits.getdata(blue_cube_path, 0, header=True)
    b_var, b_var_hdr = fits.getdata(blue_cube_path, 1, header=True)

    # Red arm
    red_sci, r_sci_hdr = fits.getdata(red_cube_path, 0, header=True)
    r_var, r_var_hdr = fits.getdata(red_cube_path, 1, header=True)

    # Average all the data
    collapsed_cube = collapse_cube(blue_sci, red_sci)

    # Automatic source detection in the collapsed cubes (red + blue)
    npeaks = 3  # Number of peaks to be detected
    border_width = 4  # Exclueded pixels atd the edge
    mean, median, std = sigma_clipped_stats(collapsed_cube, sigma=3.0)
    threshold = median + (3.0 * std)
    detection = find_peaks(
        collapsed_cube,
        threshold,
        border_width=border_width,
        npeaks=npeaks,
        centroid_func=centroid_com,
    )

    # Detected sources positions
    positions = np.transpose((detection["x_peak"], detection["y_peak"]))

    # Set the aperture center at the detected positon and the size of the aperture
    r = 5
    source_aps = CircularAperture(positions, r=r)

    # Set the size of the annulus if needed
    if sky_sub:
        r_in = r + 1
        r_out = r_in + 1
        sky_aps = CircularAnnulus(positions, r_in=r_in, r_out=r_out)

    else:
        sky_aps = None

    # Flux extraction for the aperture at all the wavelenghts

    for index, (source_ap, sky_ap) in enumerate(zip(source_aps, sky_aps)):
        ap_index = index + 1

        blue_flux, blue_var = spect_extract(blue_sci, b_var, source_ap, sky_ap=sky_ap)
        red_flux, red_var = spect_extract(red_sci, r_var, source_ap, sky_ap=sky_ap)

        # Write out the results
        blue_output = blue_cube_path.replace("p11.fits", "ap%s.p12.fits" % ap_index)
        red_output = red_cube_path.replace("p11.fits", "ap%s.p12.fits" % ap_index)

        writeFITS(blue_flux, blue_var, b_sci_hdr, b_var_hdr, blue_output)
        writeFITS(red_flux, red_var, r_sci_hdr, r_var_hdr, red_output)

    if check_plot:
        plot_apertures(red_sci, blue_sci, source_aps, sky_aps)


# ###################################################################################################################

# def updateHeader(args): # TODO do we need updates this values in the header?

#     version="anyversion"
#     metaData=DEbass.getMetadataVersion()
#     reducedBy=args.reducedBy
#     reducedBy=args.observedBy
#     redDate=DEbass.getUTC()
#     if args.ToO is not None:
#         ToO=True
#         ToOTime = args.ToO
#     else:
#         ToO=False

#     # Update red arm
#     red=fits.open(args.redArm,mode='update')
#     red[0].header['REDUCBY']=reducedBy
#     red[0].header['OBSERVBY']=observedBy
#     red[0].header['REDDATE']=redDate
#     red[0].header['PIPELINE']=version
#     red[0].header['METADATA']=metaData
#     red[0].header['TOO']=ToO
#     if ToO:
#         red[0].header['TOOTIME']=args.ToO
#     red.close()

#     # Update blue arm
#     blue=fits.open(args.blueArm,mode='update')
#     blue[0].header['REDUCBY']=reducedBy
#     blue[0].header['OBSERVBY']=observedBy
#     blue[0].header['REDDATE']=redDate
#     blue[0].header['PIPELINE']=version
#     blue[0].header['METADATA']=metaData
#     blue[0].header['TOO']=ToO
#     if ToO:
#         blue[0].header['TOOTIME']=args.ToO
#     blue.close()

#     return
