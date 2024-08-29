from astropy.io import fits as pyfits
from inspect import getargvalues, stack
import numpy
import re


# Diagnose inputs to any function
def arguments():
    """
    Report current functions source file and as-called arguments.

    Adapted from http://kbyanc.blogspot.com/2007/07/python-aggregating-function-arguments.html
    """
    this_stack = stack()
    posname, kwname, args = getargvalues(this_stack[1][0])[-3:]
    posargs = args.pop(posname, [])
    args.update(args.pop(kwname, []))
    args_string = ', '.join([str(a).strip('\'')
                             + (f"='{args[a]}'" if isinstance(args[a], str) else f"={args[a]}")
                             for a in args]).rstrip(', ')
    if posargs:
        posargs_string = ', ' + ', '.join([str(a).strip('\'')
                                           + (f"='{args[a]}'" if isinstance(args[a], str) else f"={args[a]}")
                                           for a in posargs]).rstrip(', ')
    else:
        posargs_string = ''
    return f"In file {this_stack[1][1]}, calling function {this_stack[1][3]}({args_string}{posargs_string})"


def fits_scale_from_bitpix(bitpix):
    """
    Map common BITPIX header value to astropy.io.fits 'scale' argument.
    """
    if bitpix == -64:
        return 'float64'
    elif bitpix == -32:
        return 'float32'
    elif bitpix == 32:
        return 'int32'
    elif bitpix == 16:
        return 'int16'
    elif bitpix == 8:
        return 'uint8'
    else:
        print(f"BITPIX value {bitpix} has no mapping, using current data type")
        return None


# ------------------------------------------------------------------------
# high-level functions to check if an observation is half-frame or N+S
def is_halfframe(inimg, data_hdu=0):
    """
    Report whether this exposure (filename or HDUList) is a half-frame (a.k.a. Stellar mode) image.
    """
    if isinstance(inimg, str):
        extnum = data_hdu + 1 if re.search('.fz', inimg) else data_hdu
        header = pyfits.getheader(inimg, ext=extnum)
        detsec = header["DETSEC"]
    elif isinstance(inimg, pyfits.hdu.hdulist.HDUList):
        f = inimg
        detsec = f[data_hdu].header["DETSEC"]
    else:
        raise ValueError(f"is_halfframe takes filepath or HDUList as inputs, not type {type(inimg)}")
    ystart, ystop = [int(pix) for pix in detsec.split(",")[1].rstrip(']').split(":")]
    return ystop - ystart + 1 == 2056


def is_nodshuffle(inimg, data_hdu=0):
    """
    Report whether this exposure is a Nod & Shuffle frame.
    """
    f = pyfits.open(inimg)
    ns = f[data_hdu].header["WIFESOBS"]
    f.close()
    return ns == "NodAndShuffle"


def is_subnodshuffle(inimg, data_hdu=0):
    """
    Report whether this exposure is a SubNodAndShuffle frame.
    """
    f = pyfits.open(inimg)
    sns = f[data_hdu].header["WIFESOBS"]
    f.close()
    return sns == "SubNodAndShuffle"


def is_taros(inimg):
    """
    Report whether this exposure was taken with TAROS.
    Useful for determining whether last slit of half-frame readouts is truncated.
    """
    if isinstance(inimg, pyfits.header.Header):
        header = inimg
    elif isinstance(inimg, str):
        extnum = 1 if re.search('.fz', inimg) else 0
        header = pyfits.getheader(inimg, ext=extnum)
    elif isinstance(inimg, pyfits.hdu.hdulist.HDUList):
        extnum = 1 if re.search('.fz', inimg.filename()) else 0
        header = inimg[extnum].header
    else:
        raise TypeError(f"The is_taros function takes headers, filenames, or HDULists; given {type(inimg)}")
    return 'TAROSREQ' in header


def nan_helper(y):
    """
    Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs

    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices

    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= numpy.interp(x(nans), x(~nans), y[~nans])

    https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array
    """

    return numpy.isnan(y), lambda z: z.nonzero()[0]
