from astropy.io import fits as pyfits
import numpy
import scipy.signal
import scipy.ndimage
import scipy.interpolate
from pywifes.logger_config import custom_print
import logging

# Redirect print statements to logger
logger = logging.getLogger("PyWiFeS")
print = custom_print(logger)

from .multiprocessing_utils import get_task, map_tasks
from .wifes_imtrans import blkrep, blkavg, transform_data, detransform_data
import logging


# ------------------------------------------------------------------------
# high-level functions to check if an observation is half-frame or N+S
def _is_halfframe(hdus, data_hdu=0):
    detsec = hdus[data_hdu].header["DETSEC"]
    ystart = int(float(detsec.split(",")[1].split(":")[0]))
    return ystart == 1029


# -----------------------------------------------------------------------------
laplace_kernel = numpy.array([[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]])

growth_kernel = numpy.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])


def lacos_spec_data(
    data,
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
    ny, nx = numpy.shape(data)
    x_mg, y_mg = numpy.meshgrid(numpy.arange(nx), numpy.arange(ny))
    # set up bad pix mask that will persist through multiple iterations
    global_bpm = numpy.zeros(numpy.shape(data))

    # ------------------------------------------------------------------------
    # MULTIPLE ITERATIONS
    clean_data = 1.0 * data
    # if verbose: print 'Beginning cosmic ray search'
    for i in range(niter):
        # ------------------------------------
        # step 1 - subtract sky lines
        # if verbose: print 'Generating model of sky background for iteration %d' % (i+1)
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
        # plot_sig = sig_detrend[numpy.nonzero(sig_detrend < 1000)]
        # print numpy.mean(plot_sig), numpy.std(plot_sig)
        # import pylab
        # pylab.hist(plot_sig.flatten(), bins=1000)
        # pylab.show()
        # print squirrel

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
        # if no new CRs found, exit loop
        if len(new_bpix[0]) == 0:
            break

        # ------------------------------------
        # step 7 - interpolate over neighboring pixels to fill the bad CR pix
        all_bpix = numpy.nonzero(global_bpm)
        n_bpix = len(all_bpix[0])
        for i in range(n_bpix):
            bpy = all_bpix[0][i]
            bpx = all_bpix[1][i]
            n_inds = numpy.nonzero(
                (numpy.abs(y_mg - bpy) <= n_ny)
                * (numpy.abs(x_mg - bpx) <= n_nx)
                * (global_bpm == 0)
            )
            clean_data[bpy, bpx] = numpy.median(data[n_inds])

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
):
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


# def _get_slit_indexes(slit_index):
#     slit_hdu_index = slit_index + 1
#     dq_hdu_index = slit_hdu_index + 50
#     return slit_hdu_index, dq_hdu_index


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

    hdus = pyfits.open(in_img_filepath)
    halfframe = _is_halfframe(hdus)

    if halfframe:
        nslits = 12
        first = 7
        last = 19
    else:
        nslits = 25
        first = 1
        last = 26

    outfits = pyfits.HDUList(hdus)
    if wsol_filepath:
        wsol_hdus = pyfits.open(wsol_filepath)

    for i in range(first, last):

        orig_data = hdus[i].data

        if wsol_filepath:
            wave = wsol_hdus[i - first + 1].data
        else:
            wave = None
        clean_data, global_bpm = lacos_spec_data(
            orig_data,
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
        outfits[i].data = clean_data
        # save the bad pixel mask in the DQ extentions
        outfits[i + 50].data = global_bpm
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
    hdus = pyfits.open(in_img_filepath)
    halfframe = _is_halfframe(hdus)
    if halfframe:
        nslits = 12
        first = 7
        last = 19
    else:
        nslits = 25
        first = 1
        last = 26

    tasks = []
    if wsol_filepath:
        wsol_hdus = pyfits.open(wsol_filepath)
    for i in range(first, last):
        orig_data = hdus[i].data
        if wsol_filepath:
            wave = wsol_hdus[i - first + 1].data
        else:
            wave = None
        task = get_task(
            lacos_spec_data,
            orig_data,
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
        i_slit = first + i
        i_dq_slit = first + 50 + i
        # update the data hdu
        outfits[i_slit].data = clean_data
        # save the bad pixel mask in the DQ extention
        outfits[i_dq_slit].data = global_bpm

    outfits.writeto(out_filepath, overwrite=True)
    hdus.close()
    return
