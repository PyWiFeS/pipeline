import astropy.io.fits as fits
import numpy as np
import scipy.sparse as sp


def is_evenly_sampled(_x):
    dx = _x[1:] - _x[:-1]
    return abs((dx.max() - dx.min()) / np.median(dx)) < 1e-3


def _A_lanczos(_x, x, a=2, missing=np.nan):
    """
    Returns transformation matrix for Lanczos (sinc) interpolation

    Inputs
       _x:  x-values at which the original function is sampled
        x:  x-values to which we want to interpolate
        a:  integer order of the kernel
        missing:  value to use for missing data (default:  np.nan)
    Outputs
        A:  sparse representation of the transformation matrix
    """
    if not is_evenly_sampled(_x):
        print("_A_lanczos:  warning -- _x not evenly sampled!")
        print("Lanczos interpolation may not give good results.")
    # Figure out which x values are in range
    x = np.atleast_1d(x)
    ok = np.all([x >= _x[0], x < _x[-1]], axis=0)
    # Find bins and fractional positions within bins, and interpolate
    fi = np.array((len(_x) - 1) * (x - _x[0]) / (_x[-1] - _x[0]))
    i = np.array([int(f) for f in fi])
    # Write down the matrix for Lanczos interpolation
    EPSILON = np.finfo(np.double).eps
    # The following line of code results in a deprecation waring
    # A sparce matrix in link list format
    A = sp.lil_matrix((len(x), len(_x)))
    for m in np.arange(len(x)):
        # If x[m] is in interpolation range...
        if ok[m]:
            # Generate a range of i within the original signal
            ik = range(max(0, i[m] - a + 1), min(len(_x) - 1, i[m] + a + 1))
            # Find L(x-i) for each of these
            xk = np.pi * (fi[m] - np.array(ik) + EPSILON)
            w = a * np.sin(xk) * np.sin(xk / a) / xk**2
            A[m, ik] = w / w.sum()
        elif x[m] == _x[-1]:
            A[m, -1] = 1
        else:
            A[m, 0] = missing
    return A


# Unclear if needed
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


def joinSpectra(blueSpec, redSpec):
    # Generate the resampling matrices
    blueCRVAL1 = blueSpec.header["CRVAL1"]
    blueCDELT1 = blueSpec.header["CDELT1"]
    blueCRPIX1 = blueSpec.header["CRPIX1"]
    blueNAXIS1 = blueSpec.header["NAXIS1"]

    redCRVAL1 = redSpec.header["CRVAL1"]
    redCDELT1 = redSpec.header["CDELT1"]
    redCRPIX1 = redSpec.header["CRPIX1"]
    redNAXIS1 = redSpec.header["NAXIS1"]

    nwl = int((redCRVAL1 + redCDELT1 * redNAXIS1 - blueCRVAL1) / blueCDELT1)
    wl = (np.arange(nwl) - blueCRPIX1 + 1) * blueCDELT1 + blueCRVAL1

    # One does not need to interplolate the blue spectra if the waveelength is
    # set to the blue spectrum. However the code is robust to this.

    AB = _A_lanczos(blueSpec.wl, wl, 3).tocsr()
    AR = _A_lanczos(redSpec.wl, wl, 3).tocsr()

    flux_B = np.array(AB * blueSpec.flux).ravel()
    flux_R = np.array(AR * redSpec.flux).ravel()

    # Two sparse diagonal matrices
    diag_B = sp.dia_matrix(([blueSpec.fluxvar], [0]), shape=[blueNAXIS1, blueNAXIS1])
    diag_R = sp.dia_matrix(([redSpec.fluxvar], [0]), shape=[redNAXIS1, redNAXIS1])

    fluxvar_B = np.array((AB * diag_B * AB.T).sum(axis=1)).ravel()
    fluxvar_R = np.array((AR * diag_R * AR.T).sum(axis=1)).ravel()

    buffer = 10.0

    # Should check for wavelngth overlap

    blueOnly = np.where(wl <= redSpec.minWL + buffer)
    overlap = np.where((wl < blueSpec.maxWL - buffer) & (wl > redSpec.minWL + buffer))
    redOnly = np.where(wl >= blueSpec.maxWL - buffer)

    if redSpec.minWL > blueSpec.maxWL:
        print("WARNING: No spectral overlap")

    # Average the two taking into account the buffer region and weighting

    flux = np.zeros(len(flux_B), float)
    fluxVar = np.zeros(len(flux_B), float)

    flux[blueOnly] = flux_B[blueOnly]
    fluxVar[blueOnly] = fluxvar_B[blueOnly]

    fluxVar[overlap] = 1.0 / (1.0 / fluxvar_B[overlap] + 1.0 / fluxvar_R[overlap])
    flux[overlap] = (
        flux_B[overlap] / fluxvar_B[overlap] + flux_R[overlap] / fluxvar_R[overlap]
    ) * fluxVar[overlap]

    flux[redOnly] = flux_R[redOnly]
    fluxVar[redOnly] = fluxvar_R[redOnly]

    return flux, fluxVar


def splice_spectra(blue_spec_path, red_spec_path, output):
    """
    The main routine
    """
    # Create to instances of SingleSpec
    blueSpec = SingleSpec(blue_spec_path)
    redSpec = SingleSpec(red_spec_path)

    # Join the spectra
    flux, fluxVar = joinSpectra(blueSpec, redSpec)

    # Write out the results
    # Use the header in red arm to start with
    # Add additional blue CCD keywords as required - pending
    hdulist = fits.HDUList(fits.PrimaryHDU())
    hdulist[0].data = flux
    hdulist[0].header = redSpec.header
    hdulist[0].header["CRPIX1"] = blueSpec.header["CRPIX1"]
    hdulist[0].header["CRVAL1"] = blueSpec.header["CRVAL1"]
    hdulist[0].header["CDELT1"] = blueSpec.header["CDELT1"]
    hdulist[0].header["CTYPE1"] = "wavelength"
    hdulist[0].header["CUNIT1"] = "angstrom"

    hdr_fluxvar = fits.Header()
    hdr_fluxvar["EXTNAME"] = "VARIANCE"
    hdr_fluxvar["CRPIX1"] = blueSpec.header["CRPIX1"]
    hdr_fluxvar["CRVAL1"] = blueSpec.header["CRVAL1"]
    hdr_fluxvar["CDELT1"] = blueSpec.header["CDELT1"]
    hdr_fluxvar["CTYPE1"] = "wavelength"
    hdr_fluxvar["CUNIT1"] = "angstrom"
    hdr_fluxvar["BUNIT"] = "(count / Angstrom)^2"

    hdu_fluxvar = fits.ImageHDU(data=fluxVar, header=hdr_fluxvar)
    hdulist.append(hdu_fluxvar)

    hdulist.writeto(output, overwrite=True)
    hdulist.close()

    return



def join_cubes(blue_path, red_path):

    # Read red data and metadata
    red_flux_cube = fits.getdata(red_path, 0, header=False)
    red_fluxvar_cube = fits.getdata(red_path, 1, header=False)

    redCRVAL1 = fits.getheader(red_path)["CRVAL3"]
    redCDELT1 = fits.getheader(red_path)["CDELT3"]
    redCRPIX1 = fits.getheader(red_path)["CRPIX3"]
    redNAXIS1 = fits.getheader(red_path)["NAXIS3"]
    red_wl = (np.arange(redNAXIS1) - redCRPIX1 + 1) * redCDELT1 + redCRVAL1
    red_minWL = min(red_wl)
    red_maxWL = max(red_wl)

    # Read blue data and metadata
    blue_flux_cube = fits.getdata(blue_path, 0, header=False)
    blue_fluxvar_cube = fits.getdata(blue_path, 1, header=False)
    blueCRVAL1 = fits.getheader(blue_path)["CRVAL3"]
    blueCDELT1 = fits.getheader(blue_path)["CDELT3"]
    blueCRPIX1 = fits.getheader(blue_path)["CRPIX3"]
    blueNAXIS1 = fits.getheader(blue_path)["NAXIS3"]
    blue_wl = (np.arange(blueNAXIS1) - blueCRPIX1 + 1) * blueCDELT1 + blueCRVAL1
    blue_minWL = min(blue_wl)
    blue_maxWL = max(blue_wl)



    nwl = int((redCRVAL1 + redCDELT1 * redNAXIS1 - blueCRVAL1) / blueCDELT1)
    wl = (np.arange(nwl) - blueCRPIX1 + 1) * blueCDELT1 + blueCRVAL1

    # One does not need to interplolate the blue spectra if the waveelength is
    # set to the blue spectrum. However the code is robust to this.

    AB = _A_lanczos(blue_wl, wl, 3).tocsr()
    AR = _A_lanczos(red_wl, wl, 3).tocsr()
    

    flux_spliced_list = []
    fluxvar_spliced_list = []

    wave_dim, y_dim, x_dim  = np.shape(red_flux_cube)
    wave_dim = len(wl)

    flux_cube = np.zeros((wave_dim, y_dim, x_dim))
    fluxvar_cube = np.zeros((wave_dim, y_dim, x_dim))



    # Run over every point (i,j) in the spatial plane
    for i in range(y_dim):
        for j in range(x_dim):

            red_flux =  red_flux_cube[:, i, j]
            red_fluxvar =  red_fluxvar_cube[:, i, j]

            blue_flux =  blue_flux_cube[:, i, j]
            blue_fluxvar =  blue_fluxvar_cube[:, i, j]

            flux_R = np.array(AR * red_flux).ravel()
            flux_B = np.array(AB * blue_flux).ravel()

            # Two sparse diagonal matrices
            diag_R = sp.dia_matrix(([red_fluxvar], [0]), shape=[redNAXIS1, redNAXIS1])
            diag_B = sp.dia_matrix(([blue_fluxvar], [0]), shape=[blueNAXIS1, blueNAXIS1])

            fluxvar_B = np.array((AB * diag_B * AB.T).sum(axis=1)).ravel()
            fluxvar_R = np.array((AR * diag_R * AR.T).sum(axis=1)).ravel()

            buffer = 10.0

            # Should check for wavelength overlap

            blueOnly = np.where(wl <= red_minWL + buffer)
            overlap = np.where((wl < blue_maxWL - buffer) & (wl > red_minWL + buffer))
            redOnly = np.where(wl >= blue_maxWL - buffer)

            if red_minWL > blue_maxWL:
                print("WARNING: No spectral overlap")


            # Average the two taking into account the buffer region and weighting

            flux = np.zeros(len(flux_B), float)
            fluxVar = np.zeros(len(flux_B), float)

            flux[blueOnly] = flux_B[blueOnly]
            fluxVar[blueOnly] = fluxvar_B[blueOnly]

            fluxVar[overlap] = 1.0 / (1.0 / fluxvar_B[overlap] + 1.0 / fluxvar_R[overlap])
            flux[overlap] = (
                flux_B[overlap] / fluxvar_B[overlap] + flux_R[overlap] / fluxvar_R[overlap]
            ) * fluxVar[overlap]

            flux[redOnly] = flux_R[redOnly]
            fluxVar[redOnly] = fluxvar_R[redOnly]
            flux_cube[:,i,j] = flux
            fluxvar_cube[:,i,j] = fluxVar


    #         flux_spliced_list.append(flux)
    #         fluxvar_spliced_list.append(fluxVar)            

    # # Stack the flux vectors to form a cube
    # flux_cube = np.stack(flux_spliced_list)
    # fluxvar_cube = np.stack(fluxvar_spliced_list)


    

    # Reshape the cube to the corresponding shape
    wave_dim = len(flux)
    flux_cube = flux_cube.reshape(wave_dim, y_dim, x_dim)
    fluxvar_cube = fluxvar_cube.reshape(wave_dim, y_dim, x_dim)

    return flux_cube, fluxvar_cube



def splice_cubes(blue_path, red_path, output):
    """
    The main routine
    """

    red_header = fits.getheader(red_path, 0)
    blue_header = fits.getheader(blue_path, 0)


    # Join the spectra
    flux, fluxVar = join_cubes(blue_path, red_path)

    # Write out the results
    # Use the header in red arm to start with
    # Add additional blue CCD keywords as required - pending
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

    hdulist.writeto(output, overwrite=True)
    hdulist.close()

    return
