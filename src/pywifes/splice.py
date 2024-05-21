import astropy.io.fits as fits
import numpy as np
import scipy.sparse as sp
from pywifes.logger_config import custom_print
import logging

# Redirect print statements to logger
logger = logging.getLogger('PyWiFeS')
print = custom_print(logger)

def calculate_wavelength_array(
    blue_CRVAL1, blue_CDELT1, blue_CRPIX1, red_CRVAL1, red_CDELT1, red_NAXIS1
):
    """
    Calculates the wavelength array for the overlap between blue and red spectra.

    Parameters:
    - blue_CRVAL1 (float): Reference value of the first wavelength for the blue spectrum.
    - blue_CDELT1 (float): Wavelength increment per pixel for the blue spectrum.
    - blue_CRPIX1 (float): Reference pixel for the blue spectrum.
    - red_CRVAL1 (float): Reference value of the first wavelength for the red spectrum.
    - red_CDELT1 (float): Wavelength increment per pixel for the red spectrum.
    - red_NAXIS1 (int): Number of pixels in the red spectrum.

    Returns:
    - np.ndarray: Wavelength array calculated from the given parameters.
    """
    # Calculate the number of wavelength points needed
    nwl = int((red_CRVAL1 + red_CDELT1 * red_NAXIS1 - blue_CRVAL1) / blue_CDELT1)
    # Generate the wavelength array
    wl = (np.arange(nwl) - blue_CRPIX1 + 1) * blue_CDELT1 + blue_CRVAL1
    return wl


def is_evenly_sampled(x_original):
    dx = x_original[1:] - x_original[:-1]
    return abs((dx.max() - dx.min()) / np.median(dx)) < 1e-3


def a_lanczos(x_original, x, a=2, missing=np.nan):
    """
    Returns transformation matrix for Lanczos (sinc) interpolation

    Inputs
        x_original:  x-values at which the original function is sampled
        x:  x-values to which we want to interpolate
        a:  integer order of the kernel
        missing:  value to use for missing data (default:  np.nan)
    Outputs
        A:  sparse representation of the transformation matrix
    """
    if not is_evenly_sampled(x_original):
        print("a_lanczos:  warning -- x_original not evenly sampled!")
        print("Lanczos interpolation may not give good results.")
    # Figure out which x values are in range
    x = np.atleast_1d(x)
    ok = np.all([x >= x_original[0], x < x_original[-1]], axis=0)
    # Find bins and fractional positions within bins, and interpolate
    new_fractional_pos = np.array(
        (len(x_original) - 1) * (x - x_original[0]) / (x_original[-1] - x_original[0])
    )
    new_pos = np.array([int(f) for f in new_fractional_pos])
    # Write down the matrix for Lanczos interpolation
    epsilon = np.finfo(np.double).eps

    # The following line of code results in a deprecation waring
    # A sparce matrix in link list format
    A = sp.lil_matrix((len(x), len(x_original)))
    for i in range(len(x)):
        # If x[i] is in interpolation range...
        if ok[i]:
            # Generate a range of new_pos within the original signal
            ik = range(
                max(0, new_pos[i] - a + 1), min(len(x_original) - 1, new_pos[i] + a + 1)
            )
            # Find L(x-i) for each of these
            xk = np.pi * (new_fractional_pos[i] - np.array(ik) + epsilon)
            w = a * np.sin(xk) * np.sin(xk / a) / xk**2
            A[i, ik] = w / w.sum()
        elif x[i] == x_original[-1]:
            A[i, -1] = 1
        else:
            A[i, 0] = missing
    return A


class SingleSpec(object):
    """
    Class representing a single spectrum for analysis
    """

    def __init__(self, fits_path):
        self.flux, self.header = fits.getdata(fits_path, 0, header=True)
        self.wl = (
            np.arange(self.header["NAXIS1"]) - self.header["CRPIX1"] + 1
        ) * self.header["CDELT1"] + self.header["CRVAL1"]

        self.fluxvar = fits.getdata(fits_path, 1, header=False)
        self.min_wl = np.min(self.wl)
        self.max_wl = np.max(self.wl)
        return


def join_spectra(blueSpec, redSpec):

    if redSpec.min_wl > blueSpec.max_wl:
        return None, None

    else:
        # Generate the resampling matrices
        blue_CRVAL1 = blueSpec.header["CRVAL1"]
        blue_CDELT1 = blueSpec.header["CDELT1"]
        blue_CRPIX1 = blueSpec.header["CRPIX1"]
        blue_NAXIS1 = blueSpec.header["NAXIS1"]

        red_CRVAL1 = redSpec.header["CRVAL1"]
        red_CDELT1 = redSpec.header["CDELT1"]
        red_NAXIS1 = redSpec.header["NAXIS1"]

        wl = calculate_wavelength_array(
            blue_CRVAL1, blue_CDELT1, blue_CRPIX1, red_CRVAL1, red_CDELT1, red_NAXIS1
        )

        # One does not need to interplolate the blue spectra if the waveelength is
        # set to the blue spectrum. However the code is robust to this.

        AB = a_lanczos(blueSpec.wl, wl, 3).tocsr()
        AR = a_lanczos(redSpec.wl, wl, 3).tocsr()

        # Blue
        flux_B = np.array(AB * blueSpec.flux).ravel()
        diag_B = sp.dia_matrix(
            ([blueSpec.fluxvar], [0]), shape=[blue_NAXIS1, blue_NAXIS1]
        )
        fluxvar_B = np.array((AB * diag_B * AB.T).sum(axis=1)).ravel()

        # Red
        flux_R = np.array(AR * redSpec.flux).ravel()
        diag_R = sp.dia_matrix(([redSpec.fluxvar], [0]), shape=[red_NAXIS1, red_NAXIS1])
        fluxvar_R = np.array((AR * diag_R * AR.T).sum(axis=1)).ravel()

        BUFFER = 10.0

        # Check for wavelngth overlap
        blue_only = np.where(wl <= redSpec.min_wl + BUFFER)
        overlap = np.where(
            (wl < blueSpec.max_wl - BUFFER) & (wl > redSpec.min_wl + BUFFER)
        )
        red_only = np.where(wl >= blueSpec.max_wl - BUFFER)

        # Average the two taking into account the buffer region and weighting

        flux = np.zeros(len(flux_B), float)
        fluxVar = np.zeros(len(flux_B), float)

        flux[blue_only] = flux_B[blue_only]
        fluxVar[blue_only] = fluxvar_B[blue_only]

        fluxVar[overlap] = 1.0 / (1.0 / fluxvar_B[overlap] + 1.0 / fluxvar_R[overlap])
        flux[overlap] = (
            flux_B[overlap] / fluxvar_B[overlap] + flux_R[overlap] / fluxvar_R[overlap]
        ) * fluxVar[overlap]

        flux[red_only] = flux_R[red_only]
        fluxVar[red_only] = fluxvar_R[red_only]

        return flux, fluxVar


def splice_spectra(blue_spec_path, red_spec_path, output_path):
    """
    The main routine
    """
    # Create to instances of SingleSpec
    blueSpec = SingleSpec(blue_spec_path)
    redSpec = SingleSpec(red_spec_path)

    # Join the spectra
    print("Splicing spectra.")
    flux, fluxVar = join_spectra(blueSpec, redSpec)

    if flux is None:
        print("No spectral overlap.")

    else:
        # Write out the results
        # Use the header in red arm to start with
        # Add additional blue CCD keywords as required
        hdulist = fits.HDUList(fits.PrimaryHDU())
        hdulist[0].data = flux
        hdulist[0].header = redSpec.header.copy()
        hdulist[0].header["CRPIX1"] = blueSpec.header["CRPIX1"]
        hdulist[0].header["CRVAL1"] = blueSpec.header["CRVAL1"]
        hdulist[0].header["CDELT1"] = blueSpec.header["CDELT1"]
        hdulist[0].header["CTYPE1"] = "WAVE"
        hdulist[0].header["CUNIT1"] = "Angstrom"

        hdr_fluxvar = fits.Header()
        hdr_fluxvar["EXTNAME"] = "VARIANCE"
        hdr_fluxvar["CRPIX1"] = blueSpec.header["CRPIX1"]
        hdr_fluxvar["CRVAL1"] = blueSpec.header["CRVAL1"]
        hdr_fluxvar["CDELT1"] = blueSpec.header["CDELT1"]
        hdr_fluxvar["CTYPE1"] = "WAVE"
        hdr_fluxvar["CUNIT1"] = "Angstrom"
        hdr_fluxvar["BUNIT"] = "(count / Angstrom)^2"

        hdu_fluxvar = fits.ImageHDU(data=fluxVar, header=hdr_fluxvar)
        hdulist.append(hdu_fluxvar)

        print("Saving spliced spectra.")
        hdulist.writeto(output_path, overwrite=True)
        hdulist.close()


def join_cubes(blue_path, red_path):

    # Read red data and metadata
    red_flux_cube, red_header = fits.getdata(red_path, 0, header=True)
    red_fluxvar_cube = fits.getdata(red_path, 1, header=False)
    red_CRVAL1 = red_header["CRVAL3"]
    red_CDELT1 = red_header["CDELT3"]
    red_CRPIX1 = red_header["CRPIX3"]
    red_NAXIS1 = red_header["NAXIS3"]
    red_wavelength = (np.arange(red_NAXIS1) - red_CRPIX1 + 1) * red_CDELT1 + red_CRVAL1
    red_min_wavelength = min(red_wavelength)

    # Read blue data and metadata
    blue_flux_cube, blue_header = fits.getdata(blue_path, 0, header=True)
    blue_fluxvar_cube = fits.getdata(blue_path, 1, header=False)
    blue_CRVAL1 = blue_header["CRVAL3"]
    blue_CDELT1 = blue_header["CDELT3"]
    blue_CRPIX1 = blue_header["CRPIX3"]
    blue_NAXIS1 = blue_header["NAXIS3"]
    blue_wavelength = (
        np.arange(blue_NAXIS1) - blue_CRPIX1 + 1
    ) * blue_CDELT1 + blue_CRVAL1
    blue_max_wavelength = max(blue_wavelength)

    # No wavelength overlap
    if red_min_wavelength > blue_max_wavelength:
        return None, None

    # Cubes overlap in wavelength
    wl = calculate_wavelength_array(
        blue_CRVAL1, blue_CDELT1, blue_CRPIX1, red_CRVAL1, red_CDELT1, red_NAXIS1
    )

    # One does not need to interplolate the blue spectra if the waveelength is
    # set to the blue spectrum. However the code is robust to this.

    AB = a_lanczos(blue_wavelength, wl, 3).tocsr()
    AR = a_lanczos(red_wavelength, wl, 3).tocsr()

    wave_dim, y_dim, x_dim = np.shape(red_flux_cube)
    wave_dim = len(wl)

    flux_cube = np.zeros((wave_dim, y_dim, x_dim))
    fluxvar_cube = np.zeros((wave_dim, y_dim, x_dim))

    # Run over every point (i,j) in the spatial plane
    for i in range(y_dim):
        for j in range(x_dim):
            # Red
            red_flux = red_flux_cube[:, i, j]
            flux_R = np.array(AR * red_flux).ravel()

            red_fluxvar = red_fluxvar_cube[:, i, j]
            diag_R = sp.dia_matrix(([red_fluxvar], [0]), shape=[red_NAXIS1, red_NAXIS1])
            fluxvar_R = np.array((AR * diag_R * AR.T).sum(axis=1)).ravel()

            # Blue
            blue_flux = blue_flux_cube[:, i, j]
            flux_B = np.array(AB * blue_flux).ravel()

            blue_fluxvar = blue_fluxvar_cube[:, i, j]
            diag_B = sp.dia_matrix(
                ([blue_fluxvar], [0]), shape=[blue_NAXIS1, blue_NAXIS1]
            )
            fluxvar_B = np.array((AB * diag_B * AB.T).sum(axis=1)).ravel()

            BUFFER = 10.0

            # Check for wavelength overlap
            blue_only = np.where(wl <= red_min_wavelength + BUFFER)
            overlap = np.where(
                (wl < blue_max_wavelength - BUFFER) & (wl > red_min_wavelength + BUFFER)
            )
            red_only = np.where(wl >= blue_max_wavelength - BUFFER)

            # Average the two taking into account the buffer region and weighting
            flux = np.zeros(len(flux_B), float)
            fluxVar = np.zeros(len(flux_B), float)

            flux[blue_only] = flux_B[blue_only]
            fluxVar[blue_only] = fluxvar_B[blue_only]

            fluxVar[overlap] = 1.0 / (
                1.0 / fluxvar_B[overlap] + 1.0 / fluxvar_R[overlap]
            )
            flux[overlap] = (
                flux_B[overlap] / fluxvar_B[overlap]
                + flux_R[overlap] / fluxvar_R[overlap]
            ) * fluxVar[overlap]

            flux[red_only] = flux_R[red_only]
            fluxVar[red_only] = fluxvar_R[red_only]
            flux_cube[:, i, j] = flux
            fluxvar_cube[:, i, j] = fluxVar

    # Reshape the cube to the corresponding shape
    wave_dim = len(flux)
    flux_cube = flux_cube.reshape(wave_dim, y_dim, x_dim)
    fluxvar_cube = fluxvar_cube.reshape(wave_dim, y_dim, x_dim)

    return flux_cube, fluxvar_cube


def splice_cubes(blue_path, red_path, output_path):
    """
    Main routine
    """

    red_header = fits.getheader(red_path, 0)
    blue_header = fits.getheader(blue_path, 0)

    # Join the cubes
    print("Splicing cubes.")
    flux, fluxVar = join_cubes(blue_path, red_path)

    if flux is None:
        print("No spectral overlap between the cubes.")

    else:

        # Write out the results. Use the header in red arm to start with then add additional blue CCD keywords as required
        hdulist = fits.HDUList(fits.PrimaryHDU())
        hdulist[0].data = flux
        hdulist[0].header = red_header
        hdulist[0].header["CRPIX3"] = blue_header["CRPIX3"]
        hdulist[0].header["CRVAL3"] = blue_header["CRVAL3"]
        hdulist[0].header["CDELT3"] = blue_header["CDELT3"]
        hdulist[0].header["CTYPE3"] = "wavelength"
        hdulist[0].header["CUNIT3"] = "angstrom"

        hdr_fluxvar = fits.Header()
        hdr_fluxvar["EXTNAME"] = "VARIANCE"
        hdr_fluxvar["CRPIX3"] = blue_header["CRPIX3"]
        hdr_fluxvar["CRVAL3"] = blue_header["CRVAL3"]
        hdr_fluxvar["CDELT3"] = blue_header["CDELT3"]
        hdr_fluxvar["CTYPE3"] = "wavelength"
        hdr_fluxvar["CUNIT3"] = "angstrom"
        hdr_fluxvar["BUNIT"] = "(count / Angstrom)^2"

        hdu_fluxvar = fits.ImageHDU(data=fluxVar, header=hdr_fluxvar)
        hdulist.append(hdu_fluxvar)

        print("Saving spliced cube.")
        hdulist.writeto(output_path, overwrite=True)
        hdulist.close()
