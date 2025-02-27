import astropy.io.fits as fits
import numpy
import scipy.sparse as sp


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
    - numpy.ndarray: Wavelength array calculated from the given parameters.
    """
    # Calculate the number of wavelength points needed
    nwl = int((red_CRVAL1 + red_CDELT1 * red_NAXIS1 - blue_CRVAL1) / blue_CDELT1)
    # Generate the wavelength array
    wl = (numpy.arange(nwl) - blue_CRPIX1 + 1) * blue_CDELT1 + blue_CRVAL1
    return wl


def is_evenly_sampled(x_original):
    dx = x_original[1:] - x_original[:-1]
    return abs((dx.max() - dx.min()) / numpy.median(dx)) < 1e-3


def a_lanczos(x_original, x, a=2, missing=numpy.nan):
    """
    Returns transformation matrix for Lanczos (sinc) interpolation

    Inputs
        x_original:  x-values at which the original function is sampled
        x:  x-values to which we want to interpolate
        a:  integer order of the kernel
        missing:  value to use for missing data (default:  numpy.nan)
    Outputs
        A:  sparse representation of the transformation matrix
    """
    if not is_evenly_sampled(x_original):
        print("a_lanczos:  warning -- x_original not evenly sampled!")
        print("Lanczos interpolation may not give good results.")
    # Figure out which x values are in range
    x = numpy.atleast_1d(x)
    ok = numpy.all([x >= x_original[0], x < x_original[-1]], axis=0)
    # Find bins and fractional positions within bins, and interpolate
    new_fractional_pos = numpy.array(
        (len(x_original) - 1) * (x - x_original[0]) / (x_original[-1] - x_original[0])
    )
    new_pos = numpy.array([int(f) for f in new_fractional_pos])
    # Write down the matrix for Lanczos interpolation
    epsilon = numpy.finfo(numpy.double).eps

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
            xk = numpy.pi * (new_fractional_pos[i] - numpy.array(ik) + epsilon)
            w = a * numpy.sin(xk) * numpy.sin(xk / a) / xk**2
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
        self.flux, self.header = fits.getdata(fits_path, extname='SCI', header=True)
        self.wl = (
            numpy.arange(self.header["NAXIS1"]) - self.header["CRPIX1"] + 1
        ) * self.header["CDELT1"] + self.header["CRVAL1"]
        self.fluxvar = fits.getdata(fits_path, extname='VAR', header=False)
        try:
            self.dq = fits.getdata(fits_path, extname='DQ', header=False)
        except KeyError:
            self.dq = None
        try:
            self.sky = fits.getdata(fits_path, extname='SKY', header=False)
        except KeyError:
            self.sky = None

        self.min_wl = numpy.min(self.wl)
        self.max_wl = numpy.max(self.wl)
        return


def join_spectra(blueSpec, redSpec, get_dq=False):

    if redSpec.min_wl > blueSpec.max_wl:
        return None, None, None

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
        flux_B = numpy.array(AB * blueSpec.flux).ravel()
        diag_B = sp.dia_matrix(
            ([blueSpec.fluxvar], [0]), shape=[blue_NAXIS1, blue_NAXIS1]
        ).astype(numpy.float64)
        fluxvar_B = numpy.array((AB * diag_B * AB.T).sum(axis=1)).ravel()

        # Red
        flux_R = numpy.array(AR * redSpec.flux).ravel()
        diag_R = sp.dia_matrix(
            ([redSpec.fluxvar], [0]), shape=[red_NAXIS1, red_NAXIS1]
        ).astype(numpy.float64)
        fluxvar_R = numpy.array((AR * diag_R * AR.T).sum(axis=1)).ravel()

        BUFFER = 10.0

        # Check for wavelngth overlap
        blue_only = numpy.where(wl <= redSpec.min_wl + BUFFER)
        overlap = numpy.where(
            (wl < blueSpec.max_wl - BUFFER) & (wl > redSpec.min_wl + BUFFER)
        )
        red_only = numpy.where(wl >= blueSpec.max_wl - BUFFER)

        # Average the two taking into account the buffer region and weighting

        flux = numpy.zeros(len(flux_B), dtype=numpy.float32)
        fluxVar = numpy.zeros(len(flux_B), dtype=numpy.float32)

        flux[blue_only] = flux_B[blue_only]
        fluxVar[blue_only] = fluxvar_B[blue_only]

        fluxVar[overlap] = 1.0 / (1.0 / fluxvar_B[overlap] + 1.0 / fluxvar_R[overlap])
        flux[overlap] = (
            flux_B[overlap] / fluxvar_B[overlap] + flux_R[overlap] / fluxvar_R[overlap]
        ) * fluxVar[overlap]

        flux[red_only] = flux_R[red_only]
        fluxVar[red_only] = fluxvar_R[red_only]

        # Ensure no bad effects from Lanczos interpolation.
        fluxVar = numpy.abs(fluxVar)

        if get_dq:
            dq = numpy.zeros(len(flux_B), dtype=numpy.int16)
            dq_B = numpy.array(AB * blueSpec.dq).ravel()
            dq_R = numpy.array(AR * redSpec.dq).ravel()
            dq[blue_only] = dq_B[blue_only]
            dq[overlap] = numpy.amin((dq_B[overlap], dq_R[overlap]), axis=0)
            dq[red_only] = dq_R[red_only]
        else:
            dq = None

        if blueSpec.sky is not None and redSpec.sky is not None:
            sky = numpy.zeros(len(flux_B), dtype=numpy.float32)
            sky_B = numpy.array(AB * blueSpec.sky).ravel()
            sky_R = numpy.array(AR * redSpec.sky).ravel()
            sky[blue_only] = sky_B[blue_only]
            sky[overlap] = numpy.nanmean((sky_B[overlap], sky_R[overlap]), axis=0)
            sky[red_only] = sky_R[red_only]
        else:
            sky = None

        return flux, fluxVar, dq, sky


def splice_spectra(blue_spec_path, red_spec_path, output_path, get_dq=False):
    """
    Main routine to splice two spectra together.

    Parameters
    ----------
    blue_spec_path : str
        Path to the blue spectrum FITS file.
    red_spec_path : str
        Path to the red spectrum FITS file.
    output_path : str
        Path to save the spliced spectrum FITS file.
    get_dq : bool, optional
        Whether to include the DQ extension in the output.
        Default: False.

    Returns
    -------
    None

    Notes
    -----
    This function reads two spectra from FITS files, joins them together, and saves the
    spliced spectrum to a new FITS file.

    The splicing process involves combining the flux data from the blue and red spectra,
    and updating the header information of the red spectrum to match the blue spectrum.

    If there is no spectral overlap between the blue and red spectra, a message will be
    printed indicating the lack of overlap.

    The spliced spectrum is saved as a FITS file with the provided output path.

    Examples
    --------
    >>> splice_spectra('blue_spectrum.fits', 'red_spectrum.fits', 'spliced_spectrum.fits')
    """
    # Create to instances of SingleSpec
    blueSpec = SingleSpec(blue_spec_path)
    redSpec = SingleSpec(red_spec_path)

    # Join the spectra
    print("Splicing spectra.")
    flux, fluxVar, dq, sky = join_spectra(blueSpec, redSpec, get_dq=get_dq)

    if flux is None:
        print("No spectral overlap.")

    else:
        # Write out the results
        # Use the header in red arm to start with
        # Add additional blue CCD keywords as required
        hdulist = fits.HDUList(fits.PrimaryHDU())
        hdulist[0].data = flux.astype("float32", casting="same_kind")
        hdulist[0].scale("float32")

        hdulist[0].header = redSpec.header.copy()
        hdulist[0].header["CRPIX1"] = blueSpec.header["CRPIX1"]
        hdulist[0].header["CRVAL1"] = blueSpec.header["CRVAL1"]
        hdulist[0].header["CDELT1"] = blueSpec.header["CDELT1"]
        hdulist[0].header["CTYPE1"] = "WAVE"
        hdulist[0].header["CUNIT1"] = "Angstrom"

        hdr_fluxvar = fits.Header()
        hdr_fluxvar["EXTNAME"] = "VAR"
        hdr_fluxvar["CRPIX1"] = blueSpec.header["CRPIX1"]
        hdr_fluxvar["CRVAL1"] = blueSpec.header["CRVAL1"]
        hdr_fluxvar["CDELT1"] = blueSpec.header["CDELT1"]
        hdr_fluxvar["CTYPE1"] = "WAVE"
        hdr_fluxvar["CUNIT1"] = "Angstrom"
        hdr_fluxvar["BUNIT"] = "(count / Angstrom)^2"

        hdu_fluxvar = fits.ImageHDU(data=fluxVar.astype("float32", casting="same_kind"),
                                    header=hdr_fluxvar)
        hdu_fluxvar.scale("float32")
        hdulist.append(hdu_fluxvar)

        if get_dq:
            hdr_dq = fits.Header()
            hdr_dq["EXTNAME"] = "DQ"
            hdr_dq["CRPIX1"] = blueSpec.header["CRPIX1"]
            hdr_dq["CRVAL1"] = blueSpec.header["CRVAL1"]
            hdr_dq["CDELT1"] = blueSpec.header["CDELT1"]
            hdr_dq["CTYPE1"] = "WAVE"
            hdr_dq["CUNIT1"] = "Angstrom"

            hdu_dq = fits.ImageHDU(data=dq.astype("int16", casting="unsafe"),
                                   header=hdr_dq)
            hdu_dq.scale('int16')
            hdulist.append(hdu_dq)

        if sky is not None:
            hdr_sky = fits.Header()
            hdr_sky["EXTNAME"] = "SKY"
            hdr_sky["CRPIX1"] = blueSpec.header["CRPIX1"]
            hdr_sky["CRVAL1"] = blueSpec.header["CRVAL1"]
            hdr_sky["CDELT1"] = blueSpec.header["CDELT1"]
            hdr_sky["CTYPE1"] = "WAVE"
            hdr_sky["CUNIT1"] = "Angstrom"

            hdu_sky = fits.ImageHDU(data=sky.astype("float32", casting="same_kind"),
                                    header=hdr_sky)
            hdu_sky.scale("float32")
            hdulist.append(hdu_sky)

        print("Saving spliced spectra.")
        hdulist.writeto(output_path, overwrite=True)
        hdulist.close()


def join_cubes(blue_path, red_path, get_dq=False):

    # Read red data and metadata
    red_flux_cube, red_header = fits.getdata(red_path, extname='SCI', header=True)
    red_fluxvar_cube = fits.getdata(red_path, extname='VAR', header=False)
    if get_dq:
        try:
            red_dq_cube = fits.getdata(red_path, extname='DQ', header=False)
        except IndexError:
            red_dq_cube = None
    red_CRVAL1 = red_header["CRVAL3"]
    red_CDELT1 = red_header["CDELT3"]
    red_CRPIX1 = red_header["CRPIX3"]
    red_NAXIS1 = red_header["NAXIS3"]
    red_wavelength = (numpy.arange(red_NAXIS1) - red_CRPIX1 + 1) * red_CDELT1 + red_CRVAL1
    red_min_wavelength = min(red_wavelength)

    # Read blue data and metadata
    blue_flux_cube, blue_header = fits.getdata(blue_path, extname='SCI', header=True)
    blue_fluxvar_cube = fits.getdata(blue_path, extname='VAR', header=False)
    if get_dq:
        try:
            blue_dq_cube = fits.getdata(blue_path, extname='DQ', header=False)
        except IndexError:
            blue_dq_cube = None
    blue_CRVAL1 = blue_header["CRVAL3"]
    blue_CDELT1 = blue_header["CDELT3"]
    blue_CRPIX1 = blue_header["CRPIX3"]
    blue_NAXIS1 = blue_header["NAXIS3"]
    blue_wavelength = (
        numpy.arange(blue_NAXIS1) - blue_CRPIX1 + 1
    ) * blue_CDELT1 + blue_CRVAL1
    blue_max_wavelength = max(blue_wavelength)

    # No wavelength overlap
    if red_min_wavelength > blue_max_wavelength:
        return None, None, None

    # Cubes overlap in wavelength
    wl = calculate_wavelength_array(
        blue_CRVAL1, blue_CDELT1, blue_CRPIX1, red_CRVAL1, red_CDELT1, red_NAXIS1
    )

    # One does not need to interplolate the blue spectra if the waveelength is
    # set to the blue spectrum. However the code is robust to this.

    AB = a_lanczos(blue_wavelength, wl, 3).tocsr()
    AR = a_lanczos(red_wavelength, wl, 3).tocsr()

    wave_dim, y_dim, x_dim = numpy.shape(red_flux_cube)
    wave_dim = len(wl)

    flux_cube = numpy.zeros((wave_dim, y_dim, x_dim))
    fluxvar_cube = numpy.zeros((wave_dim, y_dim, x_dim))
    if get_dq:
        dq_cube = numpy.zeros((wave_dim, y_dim, x_dim))
    else:
        dq_cube = None

    # Check for wavelength overlap
    BUFFER = 10.0
    blue_only = numpy.where(wl <= red_min_wavelength + BUFFER)
    overlap = numpy.where(
        (wl < blue_max_wavelength - BUFFER) & (wl > red_min_wavelength + BUFFER)
    )
    red_only = numpy.where(wl >= blue_max_wavelength - BUFFER)

    # Run over every point (i,j) in the spatial plane
    for i in range(y_dim):
        for j in range(x_dim):
            # Red
            red_flux = red_flux_cube[:, i, j]
            flux_R = numpy.array(AR * red_flux).ravel()

            red_fluxvar = red_fluxvar_cube[:, i, j]
            diag_R = sp.dia_matrix(
                ([red_fluxvar], [0]), shape=[red_NAXIS1, red_NAXIS1]
            ).astype(numpy.float64)
            fluxvar_R = numpy.array((AR * diag_R * AR.T).sum(axis=1)).ravel()

            # Blue
            blue_flux = blue_flux_cube[:, i, j]
            flux_B = numpy.array(AB * blue_flux).ravel()

            blue_fluxvar = blue_fluxvar_cube[:, i, j]
            diag_B = sp.dia_matrix(
                ([blue_fluxvar], [0]), shape=[blue_NAXIS1, blue_NAXIS1]
            ).astype(numpy.float64)
            fluxvar_B = numpy.array((AB * diag_B * AB.T).sum(axis=1)).ravel()

            # Average the two taking into account the buffer region and weighting
            flux = numpy.zeros(len(flux_B), float)
            fluxVar = numpy.zeros(len(flux_B), float)

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
    if get_dq:
        for i in range(y_dim):
            for j in range(x_dim):
                dq_B = numpy.array(AB * blue_dq_cube[:, i, j]).ravel()
                dq_R = numpy.array(AR * red_dq_cube[:, i, j]).ravel()
                dq_cube[blue_only, i, j] = dq_B[blue_only]
                dq_cube[overlap, i, j] = numpy.amax((dq_B[overlap],
                                                    dq_R[overlap]), axis=0)
                dq_cube[red_only, i, j] = dq_R[red_only]
        dq_cube[dq_cube > 0.1] = 1

    return flux_cube, fluxvar_cube, dq_cube


def splice_cubes(blue_path, red_path, output_path, get_dq=False):
    """
    Main function to splice two previously generated cubes from both WiFeS arms
    together.

    Parameters
    ----------
    blue_path : str
        Path to the blue cube FITS file.
    red_path : str
        Path to the red cube FITS file.
    output_path : str
        Path to save the spliced cube FITS file.
    get_dq : bool, optional
        Whether to include the DQ extension in the output.
        Default: False.

    Returns
    -------
    None

    Notes
    -----
    This function joins the blue and red cubes together, and saves the spliced cube
    to the specified output path. It uses the header information from the red cube
    as a starting point and adds additional blue CCD keywords as required.

    If there is no spectral overlap between the cubes, a message will be printed.

    The spliced cube is saved in FITS format with the following extensions:
    - Primary HDU: Contains the spliced flux data.
    - Extension HDU: Contains the variance data.

    The wavelength information is added to the header of both HDUs.

    Examples
    --------
    >>> splice_cubes('blue_cube.fits', 'red_cube.fits', 'spliced_cube.fits')
    """

    red_header = fits.getheader(red_path, 0)
    blue_header = fits.getheader(blue_path, 0)

    # Join the cubes
    print("Splicing cubes.")
    flux, fluxVar, dq = join_cubes(blue_path, red_path, get_dq=get_dq)

    if flux is None:
        print("No spectral overlap between the cubes.")
    else:
        # Write out the results. Use the header in red arm to start with then add additional blue CCD keywords as required
        hdulist = fits.HDUList(fits.PrimaryHDU())
        hdulist[0].data = flux.astype("float32", casting="same_kind")
        hdulist[0].scale("float32")

        hdulist[0].header = red_header
        hdulist[0].header["CRPIX3"] = blue_header["CRPIX3"]
        hdulist[0].header["CRVAL3"] = blue_header["CRVAL3"]
        hdulist[0].header["CDELT3"] = blue_header["CDELT3"]
        hdulist[0].header["CTYPE3"] = "Wavelength"
        hdulist[0].header["CUNIT3"] = "Angstrom"

        hdr_fluxvar = fits.Header()
        hdr_fluxvar["EXTNAME"] = "VAR"
        hdr_fluxvar["CRPIX3"] = blue_header["CRPIX3"]
        hdr_fluxvar["CRVAL3"] = blue_header["CRVAL3"]
        hdr_fluxvar["CDELT3"] = blue_header["CDELT3"]
        hdr_fluxvar["CTYPE3"] = "Wavelength"
        hdr_fluxvar["CUNIT3"] = "Angstrom"
        hdr_fluxvar["BUNIT"] = "(count / Angstrom)^2"

        hdu_fluxvar = fits.ImageHDU(data=fluxVar.astype("float32", casting="same_kind"), header=hdr_fluxvar)
        hdu_fluxvar.scale("float32")
        hdulist.append(hdu_fluxvar)

        if get_dq:
            hdr_dq = fits.Header()
            hdr_dq["EXTNAME"] = "DQ"
            hdr_dq["CRPIX3"] = blue_header["CRPIX3"]
            hdr_dq["CRVAL3"] = blue_header["CRVAL3"]
            hdr_dq["CDELT3"] = blue_header["CDELT3"]
            hdr_dq["CTYPE3"] = "Wavelength"
            hdr_dq["CUNIT3"] = "Angstrom"

            hdu_dq = fits.ImageHDU(data=dq.astype("int16", casting="unsafe"), header=hdr_dq)
            hdu_dq.scale('int16')
            hdulist.append(hdu_dq)

        print("Saving spliced cube.")
        hdulist.writeto(output_path, overwrite=True)
        hdulist.close()
