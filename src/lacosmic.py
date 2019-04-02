from astropy.io import fits as pyfits
import numpy
import subprocess
import multiprocessing
import scipy.signal
import scipy.ndimage
import scipy.interpolate
from wifes_imtrans import blkrep, blkavg, transform_data, detransform_data

import warnings
warnings.simplefilter('ignore', scipy.ComplexWarning)

numpy.seterr(invalid='ignore')

#------------------------------------------------------------------------
# high-level functions to check if an observation is half-frame or N+S
def is_halfframe(inimg, data_hdu=0):
    f = pyfits.open(inimg)
    ccdsec = f[data_hdu].header['CCDSEC']
    f.close()
    ystart = int(float(ccdsec.split(',')[1].split(':')[0]))
    if ystart == 2049:
        return True
    else:
        return False

#-----------------------------------------------------------------------------
laplace_kernel = numpy.array(
    [[ 0.0, -1.0,  0.0],
     [-1.0,  4.0, -1.0],
     [ 0.0, -1.0,  0.0]])

growth_kernel = numpy.array(
    [[ 1.0,  1.0,  1.0],
     [ 1.0,  1.0,  1.0],
     [ 1.0,  1.0,  1.0]])

def lacos_spec_data(data,
                    gain=1.0,
                    rdnoise=0.0,
                    wave=None,
                    sig_clip = 4.0,
                    sig_frac = 0.5,
                    obj_lim  = 1.0,
                    niter = 4,
                    n_nx = 1,
                    n_ny = 5,
                    verbose=True):
    ny, nx = numpy.shape(data)
    x_mg, y_mg = numpy.meshgrid(numpy.arange(nx),
                                numpy.arange(ny))
    # set up bad pix mask that will persist through multiple iterations
    global_bpm = numpy.zeros(numpy.shape(data))
    
    #------------------------------------------------------------------------
    # MULTIPLE ITERATIONS
    clean_data = 1.0*data
    #if verbose: print 'Beginning cosmic ray search'
    for i in range(niter):
        #------------------------------------
        # step 1 - subtract sky lines
        #if verbose: print 'Generating model of sky background for iteration %d' % (i+1)
        if wave is None:
            sky_model = scipy.ndimage.median_filter(clean_data, size=[7,1])
            m5_model = scipy.ndimage.median_filter(clean_data, size=[5,5])
        else:
            rect_data = transform_data(clean_data, wave)
            med_rect = scipy.ndimage.median_filter(rect_data, size=[7,1])
            m5_rect = scipy.ndimage.median_filter(rect_data, size=[5,5])
            sky_model = detransform_data(med_rect, clean_data, wave)
            m5_model  = detransform_data(m5_rect,  clean_data, wave)
        subbed_data = clean_data - sky_model
        
        #------------------------------------
        # step 2 - subtract object spectra... let's skip this
        
        #------------------------------------
        # step 3 - take 2nd order derivative (laplacian) of input image
        blkdata = blkrep(subbed_data,2,2)
        init_conv_data = scipy.signal.convolve2d(blkdata, laplace_kernel)[1:-1,1:-1]
        init_conv_data[numpy.nonzero(init_conv_data <= 0.0)] = 0.0
        conv_data = blkavg(init_conv_data,2,2)
        
        #------------------------------------
        # step 4 - make noise model, create sigma map
        noise = (1.0/gain) * ((gain*m5_model + rdnoise**2)**0.5)
        noise_min = 0.00001
        noise[numpy.nonzero(noise <= noise_min)] = noise_min
        # div by 2 to correct convolution counting
        sigmap = conv_data / (2.0*noise)
        # subtract large structure in sigma map
        sig_smooth = scipy.ndimage.median_filter(sigmap, size=[5,5])
        sig_detrend = sigmap-sig_smooth
        #plot_sig = sig_detrend[numpy.nonzero(sig_detrend < 1000)]
        #print numpy.mean(plot_sig), numpy.std(plot_sig)
        #import pylab
        #pylab.hist(plot_sig.flatten(), bins=1000)
        #pylab.show()
        #print squirrel
        
        #------------------------------------
        # step 5 - identify potential cosmic rays!!!!
        # must be greater than sig_clip*noise and obj_lim*obj_flux
        obj_flux = (
            scipy.ndimage.median_filter(subbed_data, size=[3,3])-
            scipy.ndimage.median_filter(subbed_data, size=[7,7]))
        obj_sig = obj_flux/noise
        obj_sig[numpy.nonzero(obj_sig <= 0.01)] = 0.01
        bad_pix = numpy.nonzero(
            (sig_detrend >= sig_clip)*
            (sigmap/obj_sig > obj_lim))
        new_bpm = numpy.zeros(numpy.shape(subbed_data))
        new_bpm[bad_pix] = 1
        
        #------------------------------------
        # step 6 - identify neighboring pixels that might be CRs
        neighbor_sc = sig_frac*sig_clip
        neighbor_test_bpm = scipy.signal.convolve2d(
            new_bpm, growth_kernel)[1:-1,1:-1]
        neighbor_test_bpm[numpy.nonzero(neighbor_test_bpm>=0.5)] = 1
        neighbor_sig_map = sigmap*neighbor_test_bpm
        # first check against original pixel map
        # and excise the object flux requirement
        init_neighbor_bpix = numpy.nonzero(neighbor_sig_map >= sig_clip)
        init_neighbor_bpm = numpy.zeros(numpy.shape(subbed_data))
        init_neighbor_bpm[init_neighbor_bpix] = 1
        # now check neighboring pixels with a lower threshold
        new_neighbor_bpm = scipy.signal.convolve2d(
            init_neighbor_bpm, growth_kernel)[1:-1,1:-1]
        new_neighbor_bpm[numpy.nonzero(new_neighbor_bpm>=0.5)] = 1
        new_neighbor_sig = sigmap*new_neighbor_bpm
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
            print('%d CR pixels found in iteration %d' % (
                len(new_bpix[0]), i+1))
        # if no new CRs found, exit loop
        if len(new_bpix[0]) == 0:
            break
        
        #------------------------------------
        # step 7 - interpolate over neighboring pixels to fill the bad CR pix
        all_bpix = numpy.nonzero(global_bpm)
        n_bpix = len(all_bpix[0])
        for i in range(n_bpix):
            bpy = all_bpix[0][i]
            bpx = all_bpix[1][i]
            n_inds = numpy.nonzero(
                (numpy.abs(y_mg-bpy) <= n_ny)*
                (numpy.abs(x_mg-bpx) <= n_nx)*
                (global_bpm == 0))
            clean_data[bpy, bpx] = numpy.median(data[n_inds])

    #------------------------------------------------------
    # RETURN FINAL RESULT
    return clean_data, global_bpm

def lacos_data_savefits(data,
                        output_fn,
                        gain=1.0,
                        rdnoise=0.0,
                        wave=None,
                        sig_clip = 4.0,
                        sig_frac = 0.5,
                        obj_lim  = 1.0,
                        niter = 4,
                        n_nx = 1,
                        n_ny = 5,
                        verbose=False):
    # run lacosmic
    clean_data, global_bpm = lacos_spec_data(
        data,
        gain=gain,
        rdnoise=rdnoise,
        wave=wave,
        sig_clip = sig_clip,
        sig_frac = sig_frac,
        obj_lim  = obj_lim,
        niter = niter,
        n_nx = n_nx,
        n_ny = n_ny,
        verbose=verbose)
    # save it to a fits file
    outfits = pyfits.HDUList([pyfits.PrimaryHDU(data=clean_data)])
    new_hdu = pyfits.ImageHDU(global_bpm)
    outfits.append(new_hdu)
    outfits.writeto(output_fn, overwrite=True)
    # exit
    return

#-----------------------------------------------------------------------------
# function for doing LA Cosmic on a wifes MEF file
def lacos_wifes(inimg, outimg,
                gain=1.0,       # assume data has been scaled by its gain 
                rdnoise=5.0,
                wsol_fn=None,
                sig_clip = 4.0,
                sig_frac = 0.5,
                obj_lim  = 1.0,
                niter = 4,
                n_nx = 1,
                n_ny = 5,
                multithread=False):
    if multithread:
        lacos_wifes_multithread(
            inimg, outimg,
            gain=gain,
            rdnoise=rdnoise,
            wsol_fn=wsol_fn,
            sig_clip = sig_clip,
            sig_frac = sig_frac,
            obj_lim  = obj_lim,
            niter = niter,
            n_nx = n_nx,
            n_ny = n_ny)
    else:
        lacos_wifes_oneproc(
            inimg, outimg,
            gain=gain,
            rdnoise=rdnoise,
            wsol_fn=wsol_fn,
            sig_clip = sig_clip,
            sig_frac = sig_frac,
            obj_lim  = obj_lim,
            niter = niter,
            n_nx = n_nx,
            n_ny = n_ny)
    return

def lacos_wifes_oneproc(inimg, outimg,
                        gain=1.0,       # assume data has been scaled by its gain 
                        rdnoise=5.0,
                        wsol_fn=None,
                        sig_clip = 4.0,
                        sig_frac = 0.5,
                        obj_lim  = 1.0,
                        niter = 4,
                        n_nx = 1,
                        n_ny = 5):
    # check if halfframe
    halfframe = is_halfframe(inimg)
    if halfframe:
        nslits = 12
    else:
        nslits = 25
    # open and operate on data
    f = pyfits.open(inimg)
    outfits = pyfits.HDUList(f)
    for i in range(nslits):
        curr_hdu = i+1
        print('Starting slitlet %d of image %s' % (curr_hdu,
                                                   inimg.split('/')[-1]))
        curr_dq_hdu = 50+curr_hdu
        orig_data = f[curr_hdu].data
        if wsol_fn != None:
            f2 = pyfits.open(wsol_fn)
            wave = f2[i+1].data
            f2.close()
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
            n_ny=n_ny)
        # update the data hdu
        outfits[curr_hdu].data = clean_data
        # save the bad pixel mask in the DQ extention
        outfits[curr_dq_hdu].data = global_bpm
    outfits.writeto(outimg, overwrite=True)
    f.close()
    return

def lacos_wifes_multithread(
    inimg, outimg,
    gain=1.0,       # assume data has been scaled by its gain 
    rdnoise=5.0,
    wsol_fn=None,
    sig_clip = 4.0,
    sig_frac = 0.5,
    obj_lim  = 1.0,
    niter = 4,
    n_nx = 1,
    n_ny = 5):
    f = pyfits.open(inimg)
    # check if halfframe
    halfframe = is_halfframe(inimg)
    if halfframe:
        nslits = 12
    else:
        nslits = 25
    # open and operate on data
    # 1 - spool off thread for each set of data!
    threads = []
    for i in range(nslits):
        curr_hdu = i+1
        curr_dq_hdu = 50+curr_hdu
        orig_data = f[curr_hdu].data
        if wsol_fn != None:
            f2 = pyfits.open(wsol_fn)
            wave = f2[i+1].data
            f2.close()
        else:
            wave = None
        # HERE DO THE THREADING!!
        temp_fn = 'tmp_lacos_wifes_s%02d.fits' % (i+1)
        new_thread = multiprocessing.Process(
            target=lacos_data_savefits,
            args=(orig_data, temp_fn),
            kwargs={'gain':gain,
                    'rdnoise':rdnoise,
                    'wave':wave,
                    'sig_clip':sig_clip,
                    'sig_frac':sig_frac,
                    'obj_lim':obj_lim,
                    'niter':niter,
                    'n_nx':n_nx,
                    'n_ny':n_ny})
        new_thread.start()
        threads.append(new_thread)
    for t in threads:
        t.join()
    # 2 - gather fitted data
    outfits = pyfits.HDUList(f)
    for i in range(nslits):
        curr_hdu = i+1
        curr_dq_hdu = 50+curr_hdu
        # read in the temporary output!
        temp_fn = 'tmp_lacos_wifes_s%02d.fits' % (i+1)
        f2 = pyfits.open(temp_fn)
        clean_data = f2[0].data
        global_bpm = f2[1].data
        # update the data hdu
        outfits[curr_hdu].data = clean_data
        # save the bad pixel mask in the DQ extention
        outfits[curr_dq_hdu].data = global_bpm
        # delete the temporary file
        f2.close()
        subprocess.call(['rm', '-f', temp_fn])
    # 3 - save it!
    outfits.writeto(outimg, overwrite=True)
    f.close()
    return
