from astropy.io import fits as pyfits
import numpy
import scipy.signal
import scipy.ndimage
import scipy.interpolate

from pywifes.logger_config import custom_print
import logging

from pywifes.multiprocessing_utils import get_task, map_tasks
from pywifes.wifes_imtrans import blkrep, blkavg, transform_data, detransform_data
from pywifes.wifes_utils import arguments, is_halfframe, is_taros

# Redirect print statements to logger
logger = logging.getLogger("PyWiFeS")
print = custom_print(logger)


# -----------------------------------------------------------------------------
laplace_kernel = numpy.array([[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]])

growth_kernel = numpy.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])


def lacos_spec_data(
    data,
    input_dq=None,
    gain=1.0,
    rdnoise=0.0,
    wave=None,
    sig_clip=4.0,
    sig_frac=0.5,
    obj_lim=1.0,
    niter=4,
    n_nx=1,
    n_ny=5,
    verbose=True,
):
    """
    Perform cosmic ray rejection on spectroscopic data using the L.A.Cosmic algorithm. The L.A.Cosmic algorithm is a Laplacian cosmic ray detection algorithm that uses a noise model to identify cosmic rays. The algorithm is described in van Dokkum (2001, PASP, 113, 1420).

    Parameters
    ----------
    data : ndarray
        The input spectroscopic data.
    input_dq : ndarray
        An optional input pixel mask array.
        Default: None.
    gain : float, optional
        The gain of the detector. Default is 1.0.
    rdnoise : float, optional
        The read noise of the detector. Default is 0.0.
    wave : ndarray, optional
        The wavelength array. Default is None.
    sig_clip : float, optional
        The sigma clip threshold for identifying cosmic rays. Default is 4.0.
    sig_frac : float, optional
        The fraction of the sigma clip threshold used for identifying neighboring pixels. Default is 0.5.
    obj_lim : float, optional
        The object flux limit used for identifying cosmic rays. Default is 1.0.
    niter : int, optional
        The number of iterations for the cosmic ray rejection. Default is 4.
    n_nx : int, optional
        The number of neighboring pixels to consider in the x-direction. Default is 1.
    n_ny : int, optional
        The number of neighboring pixels to consider in the y-direction. Default is 5.
    verbose : bool, optional
        Whether to print verbose output. Default is True.

    Returns
    -------
    clean_data : ndarray
        The cosmic ray cleaned spectroscopic data.
    global_bpm : ndarray
        The bad pixel mask indicating the locations of cosmic rays.
    """
    ny, nx = numpy.shape(data)
    x_mg, y_mg = numpy.meshgrid(numpy.arange(nx), numpy.arange(ny))
    # set up bad pix mask that will persist through multiple iterations
    # retain previously flagged bad pixels if mask if supplied
    global_bpm = numpy.zeros(numpy.shape(data)) if input_dq is None else input_dq
    clean_data = 1.0 * data

    # ------------------------------------------------------------------------
    # MULTIPLE ITERATIONS
    for i in range(niter):
        # ------------------------------------
        # step 1 - subtract sky lines
        if wave is None:
            sky_model = scipy.ndimage.median_filter(clean_data, size=[7, 1])
            m5_model = scipy.ndimage.median_filter(clean_data, size=[5, 5])
        else:
            rect_data = transform_data(clean_data, wave)
            med_rect = scipy.ndimage.median_filter(rect_data, size=[7, 1])
            m5_rect = scipy.ndimage.median_filter(rect_data, size=[5, 5])
            sky_model = detransform_data(med_rect, clean_data, wave)
            m5_model = detransform_data(m5_rect, clean_data, wave)
        subbed_data = clean_data - sky_model

        # ------------------------------------
        # step 2 - subtract object spectra... let's skip this

        # ------------------------------------
        # step 3 - take 2nd order derivative (laplacian) of input image
        blkdata = blkrep(subbed_data, 2, 2)
        init_conv_data = scipy.signal.convolve2d(blkdata, laplace_kernel)[1:-1, 1:-1]
        init_conv_data[numpy.nonzero(init_conv_data <= 0.0)] = 0.0
        conv_data = blkavg(init_conv_data, 2, 2)
        # ------------------------------------
        # step 4 - make noise model, create sigma map
        noise = (1.0 / gain) * ((gain * m5_model + rdnoise**2) ** 0.5)
        noise_min = 0.00001
        noise[numpy.nonzero(noise <= noise_min)] = noise_min
        # div by 2 to correct convolution counting
        sigmap = conv_data / (2.0 * noise)
        # subtract large structure in sigma map
        sig_smooth = scipy.ndimage.median_filter(sigmap, size=[5, 5])
        sig_detrend = sigmap - sig_smooth

        # ------------------------------------
        # step 5 - identify potential cosmic rays!!!!
        # must be greater than sig_clip*noise and obj_lim*obj_flux
        obj_flux = scipy.ndimage.median_filter(
            subbed_data, size=[3, 3]
        ) - scipy.ndimage.median_filter(subbed_data, size=[7, 7])
        obj_sig = obj_flux / noise
        obj_sig[numpy.nonzero(obj_sig <= 0.01)] = 0.01
        bad_pix = numpy.nonzero(
            (sig_detrend >= sig_clip) * (sigmap / obj_sig > obj_lim)
        )
        new_bpm = numpy.zeros(numpy.shape(subbed_data))
        new_bpm[bad_pix] = 1

        # ------------------------------------
        # step 6 - identify neighboring pixels that might be CRs
        neighbor_sc = sig_frac * sig_clip
        neighbor_test_bpm = scipy.signal.convolve2d(new_bpm, growth_kernel)[1:-1, 1:-1]
        neighbor_test_bpm[numpy.nonzero(neighbor_test_bpm >= 0.5)] = 1
        neighbor_sig_map = sigmap * neighbor_test_bpm
        # first check against original pixel map
        # and excise the object flux requirement
        init_neighbor_bpix = numpy.nonzero(neighbor_sig_map >= sig_clip)
        init_neighbor_bpm = numpy.zeros(numpy.shape(subbed_data))
        init_neighbor_bpm[init_neighbor_bpix] = 1
        # now check neighboring pixels with a lower threshold
        new_neighbor_bpm = scipy.signal.convolve2d(init_neighbor_bpm, growth_kernel)[
            1:-1, 1:-1
        ]
        new_neighbor_bpm[numpy.nonzero(new_neighbor_bpm >= 0.5)] = 1
        new_neighbor_sig = sigmap * new_neighbor_bpm
        final_neighbor_bpix = numpy.nonzero(new_neighbor_sig >= neighbor_sc)
        final_neighbor_bpm = numpy.zeros(numpy.shape(subbed_data))
        final_neighbor_bpm[final_neighbor_bpix] = 1
        # make final map accounting for all confirmed CRs
        new_bpm = numpy.zeros(numpy.shape(data))
        new_bpm[bad_pix] = 1
        new_bpm[init_neighbor_bpix] = 1
        new_bpm[final_neighbor_bpix] = 1
        global_bpm[bad_pix] = 1
        global_bpm[init_neighbor_bpix] = 1
        global_bpm[final_neighbor_bpix] = 1
        new_bpix = numpy.nonzero(new_bpm)
        if verbose:
            print("%d CR pixels found in iteration %d" % (len(new_bpix[0]), i + 1))
        if len(new_bpix) > 1000:
            import pdb

            pdb.set_trace()
        # if no new CRs found, exit loop with no new pixels to interpolate over
        if len(new_bpix[0]) == 0:
            break

        # ------------------------------------
        # step 7 - interpolate over neighboring pixels to fill the bad CR pix
        all_bpix = numpy.nonzero(new_bpm)  # tends to spread original bad pixels if looping over global_bpm here
        n_bpix = len(all_bpix[0])
        for i in range(n_bpix):
            bpy = all_bpix[0][i]
            bpx = all_bpix[1][i]
            n_inds = numpy.nonzero(
                (numpy.abs(y_mg - bpy) <= n_ny)
                * (numpy.abs(x_mg - bpx) <= n_nx)
                * (global_bpm == 0)
            )
            clean_data[bpy, bpx] = numpy.nanmedian(data[n_inds])
        clean_data[numpy.isnan(clean_data)] = data[numpy.isnan(clean_data)]

    # ------------------------------------------------------
    # RETURN FINAL RESULT
    return clean_data, global_bpm


# -----------------------------------------------------------------------------
# function for doing LA Cosmic on a wifes MEF file
def lacos_wifes(
    inimg,
    outimg,
    gain=1.0,  # assume data has been scaled by its gain
    rdnoise=5.0,
    wsol_fn=None,
    sig_clip=4.0,
    sig_frac=0.5,
    obj_lim=1.0,
    niter=4,
    n_nx=1,
    n_ny=5,
    is_multithread=False,
    max_processes=-1,
    debug=False
):
    """
    Performs cosmic ray rejection on spectroscopic data using the L.A.Cosmic algorithm. The L.A.Cosmic algorithm is a Laplacian cosmic ray detection algorithm that uses a noise model to identify cosmic rays. The algorithm is described in van Dokkum (2001, PASP, 113, 1420).
    If the user aims to run in multi-threaded mode, the number of processes to use can be specified. If the user does not specify the number of processes, the function will use all available processes.


    Parameters
    ----------
    inimg : str
        Input image file path.
    outimg : str
        Output image file path.
    gain : float, optional
        Gain of the data (default is 1.0).
    rdnoise : float, optional
        Read noise of the data (default is 5.0).
    wsol_fn : str, optional
        File path to the wavelength solution file (default is None).
    sig_clip : float, optional
        Sigma clip threshold for cosmic ray detection (default is 4.0).
    sig_frac : float, optional
        Fractional detection threshold for cosmic ray detection (default is 0.5).
    obj_lim : float, optional
        Object detection limit for cosmic ray detection (default is 1.0).
    niter : int, optional
        Number of iterations for cosmic ray detection (default is 4).
    n_nx : int, optional
        Number of pixels to use in the x-direction for cosmic ray detection (default is 1).
    n_ny : int, optional
        Number of pixels to use in the y-direction for cosmic ray detection (default is 5).
    is_multithread : bool, optional
        Flag indicating whether to use multi-threading (default is False).
    max_processes : int, optional
        Maximum number of processes to use for multi-threading (default is -1, which uses all available processes).
    debug : bool
        Whether to report the parameters used in this function call.
        Default: False.

    Returns
    -------
    None
    """

    if debug:
        arguments()
    if is_multithread:
        lacos_wifes_multithread(
            inimg,
            outimg,
            gain=gain,
            rdnoise=rdnoise,
            wsol_filepath=wsol_fn,
            sig_clip=sig_clip,
            sig_frac=sig_frac,
            obj_lim=obj_lim,
            niter=niter,
            n_nx=n_nx,
            n_ny=n_ny,
            max_processes=max_processes,
        )
    else:
        lacos_wifes_oneproc(
            inimg,
            outimg,
            gain=gain,
            rdnoise=rdnoise,
            wsol_filepath=wsol_fn,
            sig_clip=sig_clip,
            sig_frac=sig_frac,
            obj_lim=obj_lim,
            niter=niter,
            n_nx=n_nx,
            n_ny=n_ny,
        )
    return


def lacos_wifes_oneproc(
    in_img_filepath,
    out_filepath,
    gain=1.0,  # assume data has been scaled by its gain
    rdnoise=5.0,
    wsol_filepath=None,
    sig_clip=4.0,
    sig_frac=0.5,
    obj_lim=1.0,
    niter=4,
    n_nx=1,
    n_ny=5,
):
    """
    Apply the L.A.Cosmic cosmic ray removal algorithm to a single WiFeS image in a single process. This function is intended to be used when the user does not want to use multi-threading.

    Parameters:
    ----------
    in_img_filepath : str
        Filepath of the input image.
    out_filepath : str
        Filepath to save the output image.
    gain : float, optional
        Gain of the detector. Default is 1.0.
    rdnoise : float, optional
        Read noise of the detector. Default is 5.0.
    wsol_filepath : str, optional
        Filepath of the wavelength solution file. Default is None.
    sig_clip : float, optional
        Sigma clipping threshold for cosmic ray detection. Default is 4.0.
    sig_frac : float, optional
        Fractional threshold for cosmic ray detection. Default is 0.5.
    obj_lim : float, optional
        Object limit for cosmic ray detection. Default is 1.0.
    niter : int, optional
        Number of iterations for cosmic ray detection. Default is 4.
    n_nx : int, optional
        Number of pixels to use in the x-direction for cosmic ray detection. Default is 1.
    n_ny : int, optional
        Number of pixels to use in the y-direction for cosmic ray detection. Default is 5.

    Returns:
    -------
    None
    """

    hdus = pyfits.open(in_img_filepath)
    halfframe = is_halfframe(hdus)

    if halfframe:
        if is_taros(hdus):
            nslits = 12
        else:
            nslits = 13
    else:
        nslits = 25

    outfits = pyfits.HDUList(hdus)
    if wsol_filepath:
        wsol_hdus = pyfits.open(wsol_filepath)

    for i in range(nslits):
        curr_hdu = i + 1
        orig_data = hdus[curr_hdu].data
        orig_dq = hdus[curr_hdu + 2 * nslits].data

        if wsol_filepath:
            wave = wsol_hdus[curr_hdu].data
        else:
            wave = None
        clean_data, global_bpm = lacos_spec_data(
            orig_data,
            input_dq=orig_dq,
            gain=gain,
            rdnoise=rdnoise,
            wave=wave,
            sig_clip=sig_clip,
            sig_frac=sig_frac,
            obj_lim=obj_lim,
            niter=niter,
            n_nx=n_nx,
            n_ny=n_ny,
            verbose=False,
        )
        # update the data hdu
        outfits[curr_hdu].data = clean_data.astype("float32", casting="same_kind")
        outfits[curr_hdu].scale("float32")
        # save the bad pixel mask in the DQ extentions
        # trim data beyond range
        global_bpm[global_bpm > 32767] = 32767
        global_bpm[global_bpm < -32768] = -32768
        outfits[curr_hdu + 2 * nslits].data = global_bpm.astype("int16", casting="unsafe")
        outfits[curr_hdu + 2 * nslits].scale("int16")
    if wsol_filepath:
        wsol_hdus.close()

    outfits.writeto(out_filepath, overwrite=True)
    hdus.close()
    return


def lacos_wifes_multithread(
    in_img_filepath,
    out_filepath,
    gain=1.0,  # assume data has been scaled by its gain
    rdnoise=5.0,
    wsol_filepath=None,
    sig_clip=4.0,
    sig_frac=0.5,
    obj_lim=1.0,
    niter=4,
    n_nx=1,
    n_ny=5,
    max_processes=-1,
):
    """
    Apply the L.A.Cosmic cosmic ray removal algorithm to WiFeS data using multiple threads.

    Parameters:
    ----------
    in_img_filepath : str
        Filepath of the input image.
    out_filepath : str
        Filepath to save the output image.
    gain : float, optional
        Gain of the data (default is 1.0).
    rdnoise : float, optional
        Read noise of the data (default is 5.0).
    wsol_filepath : str, optional
        Filepath of the wavelength solution data (default is None).
    sig_clip : float, optional
        Sigma clipping threshold for cosmic ray detection (default is 4.0).
    sig_frac : float, optional
        Fractional threshold for cosmic ray detection (default is 0.5).
    obj_lim : float, optional
        Object limit for cosmic ray detection (default is 1.0).
    niter : int, optional
        Number of iterations for cosmic ray detection (default is 4).
    n_nx : int, optional
        Number of pixels in the x-direction for cosmic ray detection (default is 1).
    n_ny : int, optional
        Number of pixels in the y-direction for cosmic ray detection (default is 5).
    max_processes : int, optional
        Maximum number of processes to use for parallelization (default is -1, which uses all available processes).

    Returns:
    -------
    None
    """
    hdus = pyfits.open(in_img_filepath)
    halfframe = is_halfframe(hdus)

    if halfframe:
        if is_taros(hdus):
            nslits = 12
        else:
            nslits = 13
    else:
        nslits = 25

    tasks = []
    if wsol_filepath:
        wsol_hdus = pyfits.open(wsol_filepath)
    for i in range(nslits):
        curr_hdu = i + 1
        orig_data = hdus[curr_hdu].data
        orig_dq = hdus[curr_hdu + 2 * nslits].data
        if wsol_filepath:
            wave = wsol_hdus[curr_hdu].data
        else:
            wave = None
        task = get_task(
            lacos_spec_data,
            orig_data,
            input_dq=orig_dq,
            gain=gain,
            rdnoise=rdnoise,
            wave=wave,
            sig_clip=sig_clip,
            sig_frac=sig_frac,
            obj_lim=obj_lim,
            niter=niter,
            n_nx=n_nx,
            n_ny=n_ny,
            verbose=False,
        )
        tasks.append(task)
    if wsol_filepath:
        wsol_hdus.close()

    results = map_tasks(tasks, max_processes=max_processes)

    outfits = pyfits.HDUList(hdus)
    for i, (clean_data, global_bpm) in enumerate(results):
        i_slit = i + 1
        i_dq_slit = i_slit + 2 * nslits

        # update the data hdu
        outfits[i_slit].data = clean_data.astype("float32", casting="same_kind")
        outfits[i_slit].scale("float32")
        # save the bad pixel mask in the DQ extentions
        # trim data beyond range
        global_bpm[global_bpm > 32767] = 32767
        global_bpm[global_bpm < -32768] = -32768
        outfits[i_dq_slit].data = global_bpm.astype("int16", casting="unsafe")
        outfits[i_dq_slit].scale("int16")

    outfits.writeto(out_filepath, overwrite=True)
    hdus.close()
    return
