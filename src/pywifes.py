from __future__ import division, print_function
from astropy.io import fits as pyfits
import numpy
import pickle
import scipy.interpolate
import multiprocessing
import subprocess
from wifes_metadata import metadata_dir
from wifes_imtrans import blkrep, blkavg, transform_data, detransform_data
import gc
import mpfit
import sys
import os
import scipy.ndimage as ndimage
import scipy.interpolate as interp
import pylab
import pdb
import matplotlib.pyplot as plt

# CODE VERSION
from wifes_metadata import __version__

#------------------------------------------------------------------------
# NEED TO OPEN / ACCESS WIFES METADATA FILE!!
f0 = open(os.path.join(metadata_dir,'basic_wifes_metadata.pkl'), 'rb')
try:
    wifes_metadata = pickle.load(f0, fix_imports=True, encoding='latin')
except:
    wifes_metadata = pickle.load(f0) # fix_imports doesn't work in python 2.7.
f0.close()
blue_slitlet_defs = wifes_metadata['blue_slitlet_defs']
red_slitlet_defs = wifes_metadata['red_slitlet_defs']
nslits = len(blue_slitlet_defs.keys())

#------------------------------------------------------------------------
# functions that operate on data rather than images
def single_centroid_prof_fit(y,
                             x=None,
                             ctr_guess=None,
                             width_guess=None,
                             return_width=False):
    N = len(y)
    if x is None:
        x = numpy.arange(N,dtype='d')
    # choose x,y subregions for this line
    if ctr_guess != None and width_guess != None:
        ifit_lo = ctr_guess-5*width_guess
        ifit_hi = ctr_guess+5*width_guess
    else:
        ifit_lo = 0
        ifit_hi = N
    xfit = x[ifit_lo:ifit_hi]
    yfit = y[ifit_lo:ifit_hi]
    new_xctr = numpy.sum(xfit*(yfit**2))/numpy.sum((yfit**2))
    new_x2 = numpy.sum(((xfit-new_xctr)**2)*(yfit**2))/numpy.sum((yfit**2))
    new_rms = new_x2**0.5
    new_sig = new_rms/2.235
    if return_width:
        return new_xctr, new_sig
    else:
        return new_xctr

#------------------------------------------------------------------------
# high-level functions to check if an observation is half-frame or N+S
def is_halfframe(inimg, data_hdu=0):
    f = pyfits.open(inimg)
    ccdsec = f[data_hdu].header['CCDSEC']
    f.close()
    ystart = int(float(ccdsec.split(',')[1].split(':')[0]))
    if ystart > 2000:
        return True
    else:
        return False

def is_nodshuffle(inimg, data_hdu=0):
    f = pyfits.open(inimg)
    ns = f[data_hdu].header['WIFESOBS']
    f.close()
    if ns == 'NodAndShuffle':
        return True
    else:
        return False

#------------------------------------------------------------------------
def imcombine(inimg_list, outimg,
              method='median',
              scale=None,
              data_hdu=0):
    # read in data from inimg_list[0] to get image size
    f = pyfits.open(inimg_list[0])
    outfits = pyfits.HDUList(f)
    orig_data = f[data_hdu].data
    orig_hdr  = f[data_hdu].header
    f.close()
    nimg = len(inimg_list)
    ny, nx = numpy.shape(orig_data)
    coadd_arr = numpy.zeros([ny, nx, nimg],dtype='d')
    coadd_arr[:,:,0] = orig_data
    # gather data for all
    exptime_list = []
    airmass_list = []
    for i in range(nimg):
        f = pyfits.open(inimg_list[i])
        new_data = f[data_hdu].data
        exptime_list.append(f[data_hdu].header['EXPTIME'])
        try:
            airmass_list.append(f[data_hdu].header['AIRMASS'])
        except:
            airmass_list.append(1.0)
        f.close()
        if scale == None:
            scale_factor = 1.0
        elif scale == 'median':
            scale_factor = numpy.median(new_data)
        elif scale == 'median_nonzero':
            nonzero_inds = numpy.nonzero(new_data > 100.0)
            scale_factor = numpy.median(new_data[nonzero_inds])
        elif scale == 'exptime':
            scale_factor = exptime_list[-1]
        else:
            raise ValueError('scaling method not yet supported')
        coadd_arr[:,:,i] = new_data/scale_factor
        # later - scale data by some value, e.g. median or exptime
        gc.collect()
    # now do median (or something else - tbd later)
    if method == 'median':
        coadd_data = numpy.median(coadd_arr,axis=2)
    elif method == 'sum':
        coadd_data = numpy.sum(coadd_arr,axis=2)
    else:
        raise ValueError('combine method not yet supported')
    outfits[data_hdu].data = coadd_data
    # fix ephemeris data if images are co-added!!!
    if method == 'sum' and scale == None:
        f1 = pyfits.open(inimg_list[0])
        first_hdr = f1[data_hdu].header
        f1.close()
        f2 = pyfits.open(inimg_list[-1])
        last_hdr = f2[data_hdu].header
        f2.close()
        # HAEND, ZDEND, EXPTIME
        outfits[data_hdu].header.set(
            'EXPTIME', sum(exptime_list))
        outfits[data_hdu].header.set(
            'LSTEND', last_hdr['LSTEND'])
        outfits[data_hdu].header.set(
            'UTCEND', last_hdr['UTCEND'])
        outfits[data_hdu].header.set(
            'HAEND', last_hdr['HAEND'])
        outfits[data_hdu].header.set(
            'ZDEND', last_hdr['ZDEND'])
        outfits[data_hdu].header.set(
            'AIRMASS',
            numpy.mean(numpy.array(airmass_list)))
    # (5) write to outfile!
    outfits[data_hdu].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    gc.collect()
    return

def imcombine_mef(inimg_list, outimg,
                  data_hdu_list = list(range(1,26)),
                  var_hdu_list = list(range(26,51)),
                  dq_hdu_list = list(range(51,76)),
                  scale=None,
                  method='median'):
    nimg = len(inimg_list)
    # read in data from inimg_list[0] to get image size
    f = pyfits.open(inimg_list[0])
    outfits = pyfits.HDUList(f)
    # eventually will have different methods for each...
    for data_hdu in data_hdu_list+var_hdu_list+dq_hdu_list:
        orig_data = f[data_hdu].data
        ny, nx = numpy.shape(orig_data)
        coadd_arr = numpy.zeros([ny, nx, nimg])
        # gather data for all
        for i in range(nimg):
            f2 = pyfits.open(inimg_list[i])
            new_data = f2[data_hdu].data
            exptime = f2[data_hdu].header['EXPTIME']
            f2.close()
            if scale == None:
                scale_factor = 1.0
            elif scale == 'median':
                scale_factor = numpy.median(new_data)
            elif scale == 'exptime':
                scale_factor = exptime
            else:
                raise ValueError('scaling method not yet supported')
            coadd_arr[:,:,i] = new_data/scale_factor
            # later - scale data by some value, e.g. median or exptime
            gc.collect()
        # now do median (or something else - tbd later)
        if method == 'median':
            coadd_data = numpy.median(coadd_arr,axis=2)
        elif method == 'sum':
            coadd_data = numpy.sum(coadd_arr,axis=2)
        else:
            raise ValueError('combine method not yet supported')
        outfits[data_hdu].data = coadd_data
        gc.collect()
    # fix ephemeris data if images are co-added!!!
    if method == 'sum' and scale == None:
        airmass_list = []
        exptime_list = []
        for i in range(nimg):
            f2 = pyfits.open(inimg_list[i])
            try:
                airmass_list.append(f2[1].header['AIRMASS'])
            except:
                pass
            exptime_list.append(f2[1].header['EXPTIME'])
            f2.close()
        f1 = pyfits.open(inimg_list[0])
        first_hdr = f1[1].header
        f1.close()
        f2 = pyfits.open(inimg_list[-1])
        last_hdr = f2[1].header
        f2.close()
        # HAEND, ZDEND, EXPTIME
        outfits[1].header.set(
            'EXPTIME', sum(exptime_list))
        try:
            outfits[1].header.set(
                'LSTEND', last_hdr['LSTEND'])
        except:
            pass
        try:
            outfits[1].header.set(
                'UTCEND', last_hdr['UTCEND'])
        except:
            pass
        try:
            outfits[1].header.set(
                'HAEND', last_hdr['HAEND'])
        except:
            pass
        try:
            outfits[1].header.set(
                'ZDEND', last_hdr['ZDEND'])
        except:
            pass
        if len(airmass_list) > 0:
            outfits[1].header.set(
                'AIRMASS',
                numpy.mean(numpy.array(airmass_list)))
    # (5) write to outfile!
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits[1].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f.close()
    gc.collect()
    return

#------------------------------------------------------------------------
def imarith_mef(inimg1, operator, inimg2, outimg):
    # check if halfframe
    halfframe = is_halfframe(inimg1)
    if halfframe:
        nslits = 12
        data_hdu_list = range(1,13)
        var_hdu_list = range(26,38)
        dq_hdu_list = range(51,63)
    else:
        nslits = 25
        data_hdu_list = range(1,26)
        var_hdu_list = range(26,51)
        dq_hdu_list = range(51,76)
    # read in data from the two images, set up the output hdu
    f1 = pyfits.open(inimg1)
    f2 = pyfits.open(inimg2)
    outfits = pyfits.HDUList(f1)
    # PART 1 - data HDUs
    for data_hdu in data_hdu_list:
        data1 = f1[data_hdu].data
        data2 = f2[data_hdu].data
        # do the desired operation
        if operator == '+':
            op_data = data1+data2
        elif operator == '-':
            op_data = data1-data2
        elif operator == '*':
            op_data = data1*data2
        elif operator == '/':
            op_data = data1/data2
        else:
            raise ValueError
        outfits[data_hdu].data = op_data
        gc.collect()
    # PART 2 - var HDUs
    # NOTE: var_hdu_list must correspond directly with data_hdu_list!!
    for i in range(len(var_hdu_list)):
        try:
            var_hdu = var_hdu_list[i]
            data_hdu = data_hdu_list[i]
            var1 = f1[var_hdu].data
            var2 = f2[var_hdu].data
            data1 = f1[data_hdu].data
            data2 = f2[data_hdu].data
        except:
            continue
        # do the desired operation
        if (operator == '+') or (operator == '-'):
            op_var = var1+var2
        elif operator == '*':
            op_var = var1*(data2**2)+var2*(data1**2)
        elif operator == '/':
            op_var = var1/(data2**2) +var2*((data1/(data2**2))**2)
        else:
            raise ValueError
        outfits[var_hdu].data = op_var
        gc.collect()
    # PART 3 - dq HDUs
    for dq_hdu in dq_hdu_list:
        try:
            dq1 = f1[dq_hdu].data
            dq2 = f2[dq_hdu].data
        except:
            continue
        # always add the DQ images!!
        op_dq = dq1+dq2
        outfits[dq_hdu].data = op_dq
        gc.collect()
    # (5) write to outfile!
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f1.close()
    f2.close()
    gc.collect()
    return

def scaled_imarith_mef(inimg1, operator, inimg2, outimg,
                       scale=None):
    # check if halfframe
    halfframe = is_halfframe(inimg1)
    if halfframe:
        nslits = 12
        data_hdu_list = range(1,13)
        var_hdu_list = range(26,38)
        dq_hdu_list = range(51,63)
    else:
        nslits = 25
        data_hdu_list = range(1,26)
        var_hdu_list = range(26,51)
        dq_hdu_list = range(51,76)
    # read in data from the two images, set up the output hdu
    f1 = pyfits.open(inimg1)
    f2 = pyfits.open(inimg2)
    outfits = pyfits.HDUList(f1)
    # calculate the scale factor!
    if type(scale) == type(1.0):
        scale_factor = scale
    elif scale == 'exptime':
        exptime1 = f1[1].header['EXPTIME']
        exptime2 = f2[1].header['EXPTIME']
        scale_factor = exptime1/exptime2
    else:
        scale_factor = 1.0
    # PART 1 - data HDUs
    for data_hdu in data_hdu_list:
        data1 = f1[data_hdu].data
        data2 = scale_factor*(f2[data_hdu].data)
        # do the desired operation
        if operator == '+':
            op_data = data1+data2
        elif operator == '-':
            op_data = data1-data2
        elif operator == '*':
            op_data = data1*data2
        elif operator == '/':
            op_data = data1/data2
        else:
            raise ValueError
        outfits[data_hdu].data = op_data
    # PART 2 - var HDUs
    # NOTE: var_hdu_list must correspond directly with data_hdu_list!!
    for i in range(len(var_hdu_list)):
        var_hdu = var_hdu_list[i]
        data_hdu = data_hdu_list[i]
        var1 = f1[var_hdu].data
        var2 = (scale_factor**2)*(f2[var_hdu].data)
        data1 = f1[data_hdu].data
        data2 = scale_factor*(f2[data_hdu].data)
        # do the desired operation
        if (operator == '+') or (operator == '-'):
            op_var = var1+var2
        elif operator == '*':
            op_var = var1*(data2**2)+var2*(data1**2)
        elif operator == '/':
            op_var = var1/(data2**2) +var2*((data1/(data2**2))**2)
        else:
            raise ValueError
        outfits[var_hdu].data = op_var
    # PART 3 - dq HDUs
    for dq_hdu in dq_hdu_list:
        dq1 = f1[dq_hdu].data
        dq2 = f2[dq_hdu].data
        # always add the DQ images!!
        op_dq = dq1+dq2
        outfits[dq_hdu].data = op_dq
    # (5) write to outfile!
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f1.close()
    f2.close()
    return

def imarith(inimg1, operator, inimg2, outimg, data_hdu=0):
    f1 = pyfits.open(inimg1)
    f2 = pyfits.open(inimg2)
    outfits = pyfits.HDUList(f1)
    data1 = f1[data_hdu].data
    data2 = f2[data_hdu].data
    # do the desired operation
    if operator == '+':
        op_data = data1+data2
    elif operator == '-':
        op_data = data1-data2
    elif operator == '*':
        op_data = data1*data2
    elif operator == '/':
        op_data = data1/data2
    else:
        raise ValueError
    outfits[data_hdu].data = op_data
    outfits[data_hdu].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f1.close()
    f2.close()
    return

def imarith_float_mef(inimg1, operator, scale, outimg):
    # check if halfframe
    halfframe = is_halfframe(inimg1)
    if halfframe:
        nslits = 12
        data_hdu_list = range(1,13)
        var_hdu_list = range(26,38)
        dq_hdu_list = range(51,63)
    else:
        nslits = 25
        data_hdu_list = range(1,26)
        var_hdu_list = range(26,51)
        dq_hdu_list = range(51,76)
    # read in data from the two images, set up the output hdu
    f1 = pyfits.open(inimg1)
    outfits = pyfits.HDUList(f1)
    # PART 1 - data HDUs
    for data_hdu in data_hdu_list:
        data1 = f1[data_hdu].data
        # do the desired operation
        if operator == '+':
            op_data = data1+scale
        elif operator == '-':
            op_data = data1-scale
        elif operator == '*':
            op_data = data1*scale
        elif operator == '/':
            op_data = data1/scale
        else:
            raise ValueError
        outfits[data_hdu].data = op_data
    # PART 2 - var HDUs
    # NOTE: var_hdu_list must correspond directly with data_hdu_list!!
    for i in range(len(var_hdu_list)):
        var_hdu = var_hdu_list[i]
        data_hdu = data_hdu_list[i]
        var1 = f1[var_hdu].data
        data1 = f1[data_hdu].data
        # do the desired operation
        if (operator == '+') or (operator == '-'):
            op_var = var1
        elif operator == '*':
            op_var = var1*(scale**2)
        elif operator == '/':
            op_var = var1/(scale**2)
        else:
            raise ValueError
        outfits[var_hdu].data = op_var
    # PART 3 - dq HDUs
    for dq_hdu in dq_hdu_list:
        dq1 = f1[dq_hdu].data
        # always add the DQ images!!
        op_dq = dq1
        outfits[dq_hdu].data = op_dq
    # (5) write to outfile!
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f1.close()
    return

def imarith_float(inimg1, operator, scale, outimg, data_hdu=0):
    return imarith_float_mef(inimg1, operator, scale, outimg,
                             data_hdu_list=[data_hdu],
                             var_hdu_list=[],
                             dq_hdu_list=[])

def imcopy(inimg, outimg):
    f = pyfits.open(inimg)
    outfits = pyfits.HDUList(f)
    outfits.writeto(outimg, overwrite=True)
    f.close()
    return

#------------------------------------------------------------------------
# OVERSCAN SUBTRACTION
#### Default values for 1x1 binning
###default_det_reg_bin1x1 = [[1,4096,62,2109],
###                          [1,4096,2115,4162]]
####default_ovs_reg_bin1x1 = [1,4096,48,61]
###default_ovs_reg_bin1x1 = [[1,4096,48,61],
###                          [1,4096,48,61]]
#### Default values for 2x binning in Y
###default_det_reg_bin1x2 = [[1,2048,62,2109],
###                          [1,2048,2115,4162]]
####default_ovs_reg_bin1x2 = [1,2048,48,61]
###default_ovs_reg_bin1x2 = [[1,2048,48,61],
###                          [1,2048,48,61]]

default_detector_values = {
    'B1' : {'det_regs' : [[1,2048,22,2069],
                          [1,2048,2070,4117],
                          [2049,4096,2070,4117],
                          [2049,4096,22,2069]],
            'ovs_regs' : [[1,2048,8,17],
                          [1,2048,4122,4131],
                          [2049,4096,4122,4131],
                          [2049,4096,8,17]],
            'sci_regs' : [[1,2048,1,2048],
                          [1,2048,2049,4096],
                          [2049,4096,2049,4096],
                          [2049,4096,1,2048]],
            'gain'     : [0.82, 0.92, 0.92, 0.96],
            'rdnoise'  : [7.9, 8.95, 4.4, 5.7]},
    'B2' : {'det_regs' : [[1,2048,19,2066],
                          [1,2048,2115,4162],
                          [2049,4096,2115,4162],
                          [2049,4096,19,2066]],
            'ovs_regs' : [[1,2048,7,14],
                          [1,2048,4166,4173],
                          [2049,4096,4166,4173],
                          [2049,4096,7,14]],
            'sci_regs' : [[1,2048,1,2048],
                          [1,2048,2049,4096],
                          [2049,4096,2049,4096],
                          [2049,4096,1,2048]],
            'gain'     : [0.82, 0.92, 0.92, 0.96],
            'rdnoise'  : [7.9, 8.95, 4.4, 5.7]},
    'B3' : {'det_regs' : [[1,4096,62,2109],
                          [1,4096,2115,4162]],
            'ovs_regs' : [[1,4096,14,40],
                          [1,4096,14,40]],
            'sci_regs' : [[1,4096,1,2048],
                          [1,4096,2049,4096]],
            'gain'     : [1.0, 1.0],
            'rdnoise'  : [7.0, 7.0]},
    'B4' : {'det_regs' : [[1,4096,54,4149],],
            'ovs_regs' : [[1,4096,6,53],],
            'sci_regs' : [[1,4096,1,4096],],
            'gain'     : [1.0],
             'rdnoise'  : [7.0]},
    'R1' : {'det_regs' : [[1,2048,22,2069],
                          [1,2048,2070,4117],
                          [2049,4096,2070,4117],
                          [2049,4096,22,2069]],
            'ovs_regs' : [[1,2048,8,17],
                          [1,2048,4122,4131],
                          [2049,4096,4122,4131],
                          [2049,4096,8,17]],
            'sci_regs' : [[1,2048,1,2048],
                          [1,2048,2049,4096],
                          [2049,4096,2049,4096],
                          [2049,4096,1,2048]],
            'gain'     : [1.01, 0.98, 1.02, 0.84],
            'rdnoise'  : [4.4, 7.0, 5.2, 5.0]},
    'R2' : {'det_regs' : [[1,2048,19,2066],
                          [1,2048,2115,4162],
                          [2049,4096,2115,4162],
                          [2049,4096,19,2066]],
            'ovs_regs' : [[1,2048,7,14],
                          [1,2048,4166,4173],
                          [2049,4096,4166,4173],
                          [2049,4096,7,14]],
            'sci_regs' : [[1,2048,1,2048],
                          [1,2048,2049,4096],
                          [2049,4096,2049,4096],
                          [2049,4096,1,2048]],
            'gain'     : [1.01, 0.98, 1.02, 0.84],
            'rdnoise'  : [4.4, 7.0, 5.2, 5.0]},
    'R3' : {'det_regs' : [[1,4096,62,2109],
                          [1,4096,2115,4162]],
            'ovs_regs' : [[1,4096,14,40],
                          [1,4096,14,40]],
            'sci_regs' : [[1,4096,1,2048],
                          [1,4096,2049,4096]],
            'gain'     : [1.0, 1.0],
            'rdnoise'  : [7.0, 7.0]},
    'R4a' : {'det_regs' : [[1,4096,64,4159],],
             'ovs_regs' : [[1,4096,1,63],],
             'sci_regs' : [[1,4096,1,4096],],
             'gain'     : [1.0],
             'rdnoise'  : [7.0]},
    'R4' : {'det_regs' : [[1,4096,54,4149],],
            'ovs_regs' : [[1,4096,6,53],],
            'sci_regs' : [[1,4096,1,4096],],
            'gain'     : [1.0],
             'rdnoise'  : [7.0]},
    }

def determine_detector_epoch(inimg, data_hdu=0):
    f = pyfits.open(inimg)
    orig_hdr = f[data_hdu].header
    f.close()
    utc_str = orig_hdr['DATE-OBS'].split('T')[0]
    utc_date = int(float(utc_str.replace('-','')))
    camera = orig_hdr['CAMERA']
    if camera == 'WiFeSRed':
        if utc_date < 20100511:
            epoch = 'R1'
        elif utc_date < 20100801:
            epoch = 'R2'
        elif utc_date < 20130225:
            epoch = 'R3'
        elif utc_date < 20130322:
            epoch = 'R4a'
        else:
            epoch = 'R4'
    else:
        if utc_date < 20100511:
            epoch = 'B1'
        elif utc_date < 20110616:
            epoch = 'B2'
        elif utc_date < 20130522:
            epoch = 'B3'
        else:
            epoch = 'B4'
    return epoch

def convert_ccd_to_bindata_pix(pix_defs, bin_x, bin_y):
    mod_pix_defs = [(pix_defs[0]-1)//bin_y,
                    (pix_defs[1]-1)//bin_y,
                    (pix_defs[2]-1)//bin_x,
                    (pix_defs[3]-1)//bin_x]
    return mod_pix_defs

def subtract_overscan(inimg, outimg,
                      data_hdu = 0,
                      detector_regions=None,
                      overscan_regions=None,
                      gain=None,
                      rdnoise=None):
    # (0) open file, initialize HDUList that will be saved at output
    f = pyfits.open(inimg)
    outfits = pyfits.HDUList(f)
    orig_data = f[data_hdu].data
    orig_hdr = f[data_hdu].header
    f.close()
    # (1) format detector and overscan regions
    #     - if 'None' grab default values
    ccdsum = orig_hdr['CCDSUM']
    bins = ccdsum.split()
    bin_x = int(float(bins[0]))
    bin_y = int(float(bins[1]))
    # default values if all of the values are not specified
    if ((detector_regions == None) or
        (overscan_regions == None) or
        (gain == None) or
        (rdnoise == None)):
        # determine the epoch
        epoch = determine_detector_epoch(inimg)
        # get detector characteristics
        init_det_reg = default_detector_values[epoch]['det_regs']
        detector_regions = [convert_ccd_to_bindata_pix(x, bin_x, bin_y)
                            for x in init_det_reg]
        init_ovs_reg = default_detector_values[epoch]['ovs_regs']
        overscan_regions = [convert_ccd_to_bindata_pix(x, bin_x, bin_y)
                            for x in init_ovs_reg]
        init_sci_reg = default_detector_values[epoch]['sci_regs']
        science_regions = [convert_ccd_to_bindata_pix(x, bin_x, bin_y)
                           for x in init_sci_reg]
        gain = default_detector_values[epoch]['gain']
        rdnoise = default_detector_values[epoch]['rdnoise']
    # FORMAT REGION DEFINITIONS
    fmt_det_reg = [[x[0],(x[1]+1),x[2],(x[3]+1)]
                   for x in detector_regions]
    fmt_ovs_reg = [[x[0],(x[1]+1),x[2],(x[3]+1)]
                   for x in overscan_regions]
    fmt_sci_reg = [[x[0],(x[1]+1),x[2],(x[3]+1)]
                   for x in science_regions]
    # (2) create data array - MUST QUERY FOR HALF-FRAME
    halfframe = is_halfframe(inimg, data_hdu=data_hdu)
    if halfframe:
        ny = 2048//bin_y
        nx = 4096//bin_x
        # note: offset will be equal to ny
        subbed_data = numpy.zeros([ny,nx], dtype=numpy.float)
        for i in range(len(fmt_det_reg)):
            init_det = fmt_det_reg[i]
            det = [init_det[0], init_det[1]-ny, init_det[2], init_det[3]]
            init_ovs = fmt_ovs_reg[i]
            ovs = [init_ovs[0], init_ovs[1]-ny, init_ovs[2], init_ovs[3]]
            init_sci = fmt_sci_reg[i]
            sci = [init_sci[0], init_sci[1]-ny, init_sci[2], init_sci[3]]
            curr_data = (
                orig_data[det[0]:det[1],det[2]:det[3]])
            # (3) create mean/median overscan, subtract from data
            curr_ovs_data = orig_data[ovs[0]:ovs[1],ovs[2]:ovs[3]]
            curr_ovs_val = numpy.median(curr_ovs_data)
            subbed_data[sci[0]:sci[1],sci[2]:sci[3]] = (
                gain[i]*(curr_data-curr_ovs_val))
    else:
        ny = 4096//bin_y
        nx = 4096//bin_x
        subbed_data = numpy.zeros([ny,nx], dtype=numpy.float)
        for i in range(len(fmt_det_reg)):
            det = fmt_det_reg[i]
            ovs = fmt_ovs_reg[i]
            sci = fmt_sci_reg[i]
            curr_data = (
                orig_data[det[0]:det[1],det[2]:det[3]])
            # (3) create mean/median overscan, subtract from data
            curr_ovs_data = orig_data[ovs[0]:ovs[1],ovs[2]:ovs[3]]
            curr_ovs_val = numpy.median(curr_ovs_data)
            subbed_data[sci[0]:sci[1],sci[2]:sci[3]] = (
                gain[i]*(curr_data-curr_ovs_val))
    # (3b) - if epoch R4a, flip data!
    if epoch == 'R4a':
        temp_data = subbed_data[:,::-1]
        next_data = temp_data[::-1,:]
        subbed_data = next_data
    # (4) only modify data part of data_hdu for output
    detsize_str = '[%d:%d,%d:%d]' % (1, ny, 1, nx)
    outfits[data_hdu].header.set('DETSIZE', detsize_str)
    outfits[data_hdu].header.set('CCDSIZE', detsize_str)
    #outfits[data_hdu].header.set('CCDSEC',  detsize_str)
    outfits[data_hdu].header.set('DATASEC', detsize_str)
    outfits[data_hdu].header.set('TRIMSEC', detsize_str)
    outfits[data_hdu].header.set('RDNOISE', max(rdnoise))
    outfits[data_hdu].header.set('GAIN', 1.0)
    outfits[data_hdu].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits[data_hdu].data = subbed_data
    # (5) write to outfile!
    outfits.writeto(outimg, overwrite=True)
    return

#------------------------------------------------------------------------
# NEW 2013-06-06
def repair_blue_bad_pix(inimg, outimg,
                        data_hdu=0,
                        bin_x=1, bin_y=1):
    # first check it is the new blue detector, otherwise skip
    epoch = determine_detector_epoch(inimg)
    if epoch[0] != 'B' or float(epoch[1]) < 4:
        imcopy(inimg, outimg)
        return
    # get data and header
    f = pyfits.open(inimg)
    outfits = pyfits.HDUList(f)
    orig_data = f[data_hdu].data
    orig_hdr = f[data_hdu].header
    f.close()
    # figure out binning
    ccdsum = orig_hdr['CCDSUM']
    bins = ccdsum.split()
    bin_x = int(float(bins[0]))
    bin_y = int(float(bins[1]))
    # now interpolate
    halfframe = is_halfframe(inimg)
    if halfframe:
        y_bad_min = 0
        y_bad_max = 2048//bin_y
    else:
        y_bad_min = 752//bin_y-1
        y_bad_max = 4096//bin_y
    x_bad = 1529//bin_x-1
    interp_data = 1.0*orig_data
    for i in range(y_bad_min, y_bad_max):
        interp_data[i,x_bad] = 0.5*(interp_data[i,x_bad-1]+
                                    interp_data[i,x_bad+1])
    # save it!
    outfits[data_hdu].data = interp_data
    outfits.writeto(outimg, overwrite=True)
    return

def repair_red_bad_pix(inimg, outimg,
                       data_hdu=0,
                       bin_x=1, bin_y=1):
    # first check it is the new blue detector, otherwise skip
    epoch = determine_detector_epoch(inimg)
    if epoch[0] != 'R' or float(epoch[1]) < 4:
        imcopy(inimg, outimg)
        return
    # get data and header
    f = pyfits.open(inimg)
    outfits = pyfits.HDUList(f)
    orig_data = f[data_hdu].data
    orig_hdr = f[data_hdu].header
    f.close()
    # figure out binning
    ccdsum = orig_hdr['CCDSUM']
    bins = ccdsum.split()
    bin_x = int(float(bins[0]))
    bin_y = int(float(bins[1]))
    # now interpolate
    halfframe = is_halfframe(inimg)
    if halfframe:
        y_bad_min = (3242-2048)//bin_y-1
        y_bad_max = 2048//bin_y
    else:
        y_bad_min = 3282//bin_y-1
        y_bad_max = 4096//bin_y
    x_bad_1 = 774//bin_x-1
    x_bad_2 = 775//bin_x-1
    interp_data = 1.0*orig_data
    for i in range(y_bad_min, y_bad_max):
        interp_data[i,x_bad_1] = (2.0*interp_data[i,x_bad_1-1]+
                                  1.0*interp_data[i,x_bad_2+1])/3.0
        interp_data[i,x_bad_2] = (1.0*interp_data[i,x_bad_1-1]+
                                  2.0*interp_data[i,x_bad_2+1])/3.0
    # save it!
    outfits[data_hdu].data = interp_data
    outfits.writeto(outimg, overwrite=True)
    return

#------------------------------------------------------------------------
# NEW 2012-11-08
# inter-slitlet bias subtraction method
def fit_wifes_interslit_bias(inimg,
                             data_hdu=0,
                             slitlet_def_file=None,
                             method='row_med',
                             x_polydeg=1,
                             y_polydeg=1,
                             plot=False):
    #---------------------------
    # (0) open file, initialize HDUList that will be saved at output
    f = pyfits.open(inimg)
    outfits = pyfits.HDUList(f)
    orig_data = f[data_hdu].data
    orig_hdr = f[data_hdu].header
    f.close()
    # check if halfframe
    halfframe = is_halfframe(inimg)
    # grab necessary info from header
    ccdsum = orig_hdr['CCDSUM']
    bins = ccdsum.split()
    bin_x = int(float(bins[0]))
    bin_y = int(float(bins[1]))
    utc_str = orig_hdr['DATE-OBS'].split('T')[0]
    utc_date = int(float(utc_str.replace('-','')))
    camera = orig_hdr['CAMERA']
    #---------------------------
    # (1) make a mask of inter-slitlet zones
    interstice_map = numpy.ones(numpy.shape(orig_data))
    interstice_mask = numpy.ones(numpy.shape(orig_data)[0])
    if slitlet_def_file != None:
        f2 = open(slitlet_def_file, 'r')
        slitlet_defs = pickle.load(f2)
        f2.close()
    elif camera == 'WiFeSRed':
        slitlet_defs = red_slitlet_defs
    else:
        slitlet_defs = blue_slitlet_defs
    if halfframe:
        nslits = 12
        offset = 4096//bin_y
    else:
        nslits = 25
        offset = 0
    for i in range(nslits):
        init_curr_defs = slitlet_defs[str(i+1)]
        # convert to binned values
        curr_defs = [
            ((init_curr_defs[0]-1)//bin_x)+1,
            ((init_curr_defs[1]-1)//bin_x)+1,
            ((init_curr_defs[2]-1)//bin_y)+1,
            ((init_curr_defs[3]-1)//bin_y)+1]
        # horrible kluge to make sure everything has the same dimensions!
        if (curr_defs[1]-curr_defs[0]+1) != (4096//bin_x):
            curr_defs[1] -= 1
        if (curr_defs[3]-curr_defs[2]+1) != (86//bin_y):
            curr_defs[3] -= 1
        # define a buffer zone
        ybuff=6//bin_y
        mod_defs = [curr_defs[0]-1, curr_defs[1],
                    curr_defs[2]-1-ybuff-offset,
                    curr_defs[3]+ybuff-offset]
        # everything in slitlet zone gets set to zero
        interstice_map[mod_defs[2]:mod_defs[3],
                       mod_defs[0]:mod_defs[1]] = 0
        interstice_mask[mod_defs[2]:mod_defs[3]] = 0
    #---------------------------
    # (2) from its epoch, determine if is 1-amp or 4-amp readout
    if camera == 'WiFeSRed':
        if utc_date < 20110616:
            # namps=4
            sci_regs = [[0,2047//bin_y,0,2047//bin_x],
                        [0,2047//bin_y,2048//bin_x,4095//bin_x],
                        [2048//bin_y,4095//bin_y,0,2047//bin_x],
                        [2048//bin_y,4095//bin_y,2048//bin_x,4095//bin_x]]
        else:
            # namps=1
            sci_regs = [[0,4095//bin_y,0,4095//bin_x]]
    else:
        if utc_date < 20100801:
            # namps=4
            sci_regs = [[0,2047//bin_y,0,2047//bin_x],
                        [0,2047//bin_y,2048//bin_x,4095//bin_x],
                        [2048//bin_y,4095//bin_y,0,2047//bin_x],
                        [2048//bin_y,4095//bin_y,2048//bin_x,4095//bin_x]]
        else:
            # namps=1
            sci_regs = [[0,4095//bin_y,0,4095//bin_x]]
    #---------------------------
    # (3) for each sci region, get interstitial bias level
    out_data = numpy.zeros(numpy.shape(orig_data), dtype='d')
    for i in range(len(sci_regs)):
        init_reg = sci_regs[i]
        if halfframe:
            ny = 2048//bin_y
            reg = [init_reg[0], init_reg[1]-ny, init_reg[2], init_reg[3]]
        else:
            reg = init_reg
        curr_data = orig_data[reg[0]:reg[1]+1,reg[2]:reg[3]+1]
        curr_map  = interstice_map[reg[0]:reg[1]+1,reg[2]:reg[3]+1]
        curr_mask = interstice_mask[reg[0]:reg[1]+1]
        ny, nx = numpy.shape(curr_data)
        linx = numpy.arange(nx, dtype='d')
        liny = numpy.arange(ny, dtype='d')
        full_x, full_y = numpy.meshgrid(linx, liny)        
        if method == 'row_med':
            #ny, nx = numpy.shape(curr_data)
            row_med = numpy.zeros(nx, dtype='d')
            # iterate over columns to excise CRs... bleh
            for i in range(nx):
                curr_col = curr_data[:,i]
                curr_med = numpy.median(curr_col[
                    numpy.nonzero(curr_mask)[0]])
                good_inds = numpy.nonzero(
                    (numpy.abs(curr_col - curr_med) < 20.0)*
                    (curr_mask))[0]
                pylab.figure()
                pylab.hist(curr_col[good_inds], bins=50)
                pylab.show()
                row_med[i] = numpy.mean(curr_col[good_inds])
            #row_med = numpy.median(curr_data, axis=0)
            bias_sub = row_med**numpy.ones(
                numpy.shape(curr_data), dtype='d')
            out_data[reg[0]:reg[1]+1,reg[2]:reg[3]+1] = bias_sub
        elif method == 'surface':
            curr_inds = numpy.nonzero(curr_map)
            fit_x = full_x[curr_inds].flatten()
            fit_y = full_y[curr_inds].flatten()
            fit_b = curr_data[curr_inds].flatten()
            from wifes_wsol import fit_wsol_poly, evaluate_wsol_poly
            init_x_poly, init_y_poly = fit_wsol_poly(fit_x, fit_y, fit_b,
                                                     x_polydeg, y_polydeg)
            init_bias_fvals = evaluate_wsol_poly(
                fit_x, fit_y, init_x_poly, init_y_poly)
            resids = fit_b - init_bias_fvals
            resids_rms = numpy.median(resids**2)**0.5
            good_inds = numpy.nonzero(numpy.abs(resids) <
                                      3.0*(resids_rms))
            x_poly, y_poly = fit_wsol_poly(fit_x[good_inds],
                                           fit_y[good_inds],
                                           fit_b[good_inds],
                                           x_polydeg, y_polydeg)
            bias_fvals = evaluate_wsol_poly(
                full_x, full_y, x_poly, y_poly)
            out_data[reg[0]:reg[1]+1,reg[2]:reg[3]+1] = bias_fvals
            # PLOT IT!!!
            if plot:
                randplot_data = numpy.random.permutation(len(fit_x))[:1000]
                randplot_inds = numpy.random.permutation(len(full_x.flatten()))[:1000]
                randplot_lin_mask = numpy.zeros(len(full_x.flatten()))
                randplot_lin_mask[randplot_inds] = 1
                randplot_mask = numpy.reshape(randplot_lin_mask, numpy.shape(full_x))
                randplot_fit  = numpy.nonzero(randplot_mask)
                import pylab
                from mpl_toolkits.mplot3d import axes3d
                fig = pylab.figure()
                sp = fig.gca(projection='3d')
                sp.plot(fit_x[randplot_data],
                        fit_y[randplot_data],
                        fit_b[randplot_data], 'b.')
                sp.plot_surface(full_x,
                                full_y,
                                bias_fvals,
                                alpha=0.3, rstride=200, cstride=200, color='r')
                fval_med = numpy.median(bias_fvals)
                sp.set_zlim([fval_med-5.0*resids_rms,
                             fval_med+5.0*resids_rms])
                pylab.show()
        elif method == 'median':
            bias_val = numpy.median(curr_data[numpy.nonzero(curr_map)])
            out_data[reg[0]:reg[1]+1,reg[2]:reg[3]+1] = bias_val
        else:
            print('only row_med and median methods currently supported')
    #---------------------------
    # (4) return it!
    return out_data
    return

def save_wifes_interslit_bias(inimg, outimg,
                              data_hdu=0,
                              slitlet_def_file=None,
                              method='surface',
                              x_polydeg=1,
                              y_polydeg=1,
                              plot=False):
    # (0) open file, initialize HDUList that will be saved at output
    f = pyfits.open(inimg)
    outfits = pyfits.HDUList(f)
    # (1) fit bias with requested settings
    out_data = fit_wifes_interslit_bias(inimg,
                                        data_hdu=data_hdu,
                                        slitlet_def_file=slitlet_def_file,
                                        method=method,
                                        x_polydeg=x_polydeg,
                                        y_polydeg=y_polydeg,
                                        plot=plot)
    # (2) save it!
    outfits[data_hdu].data = out_data
    outfits[data_hdu].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f.close()
    return

def subtract_wifes_interslit_bias(inimg, outimg,
                                  data_hdu=0,
                                  slitlet_def_file=None,
                                  method='surface',
                                  x_polydeg=1,
                                  y_polydeg=1,
                                  plot=False):
    # (0) open file, initialize HDUList that will be saved at output
    f = pyfits.open(inimg)
    outfits = pyfits.HDUList(f)
    orig_data = f[data_hdu].data
    # (1) fit bias with requested settings
    out_data = fit_wifes_interslit_bias(inimg,
                                        data_hdu=data_hdu,
                                        slitlet_def_file=slitlet_def_file,
                                        method=method,
                                        x_polydeg=x_polydeg,
                                        y_polydeg=y_polydeg,
                                        plot=plot)
    # (2) save it!
    outfits[data_hdu].data = orig_data - out_data
    outfits[data_hdu].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    return

#------------------------------------------------------------------------
# NEW 2012-10-02
# bias subtraction method

#--------------------------- Fred's update ------------------------------
def wifes_bias_model (p,x,camera) :
    if camera == 'WiFeSRed' :
        model = p[0]
        model += p[1] * numpy.exp(p[2]/numpy.abs(x-p[3]) )
        model += p[4] * numpy.exp(p[5]/numpy.abs(x-p[6]) )
        model += p[7] * x
    else :
        model = p[0]
        model += p[1] * numpy.exp(p[2]/numpy.abs(x-p[3]) )
        model += p[4] * numpy.exp(p[5]/numpy.abs(x-p[6]) )
        model += p[7] * x
        model += p[8]* numpy.exp(-(p[9]*(x-p[10]))**2)
    return model

#--------------------------- Fred's update ------------------------------
def error_wifes_bias_model (p,x,z,err,camera, fjac=None) :
    status = 0
    residual = (wifes_bias_model(p,x,camera) - z)/err
    return ([status,residual])

#--------------------------- Fred's update ------------------------------
def generate_wifes_bias_fit(bias_img, outimg, data_hdu=0,
                            plot=False, verbose=False, method = 'row_med'):
    # get object and bias data
    f1 = pyfits.open(bias_img)
    orig_data = f1[data_hdu].data
    orig_hdr  = f1[data_hdu].header
    # from its epoch, determine if is 1-amp or 4-amp readout
    ccdsum = orig_hdr['CCDSUM']
    bins = ccdsum.split()
    bin_x = int(float(bins[0]))
    bin_y = int(float(bins[1]))
    utc_str = orig_hdr['DATE-OBS'].split('T')[0]
    utc_date = int(float(utc_str.replace('-','')))
    camera = orig_hdr['CAMERA']
    if camera == 'WiFeSRed':
        if utc_date < 20110616:
            # namps=4
            sci_regs = [[0,2047//bin_y,0,2047//bin_x],
                        [0,2047//bin_y,2048//bin_x,4095//bin_x],
                        [2048//bin_y,4095//bin_y,0,2047//bin_x],
                        [2048//bin_y,4095//bin_y,2048//bin_x,4095//bin_x]]
            # Currently not supported !
            if method == 'fit' :
                print(' ---')
                print('Bias frame not compatible with surface fitting (4 amps) !')
                print(' --- ')
        else:
            # namps=1
            sci_regs = [[0,4095//bin_y,0,4095//bin_x]]
    else:
        if utc_date < 20100801:
            # namps=4
            sci_regs = [[0,2047//bin_y,0,2047//bin_x],
                        [0,2047//bin_y,2048//bin_x,4095//bin_x],
                        [2048//bin_y,4095//bin_y,0,2047//bin_x],
                        [2048//bin_y,4095//bin_y,2048//bin_x,4095//bin_x]]
            # Currently not supported !
            if method == 'fit' :
                print(' ---')
                print('Bias frame not compatible with surface fitting (4 amps) !')
                print(' --- ')

        else:
            # namps=1
            sci_regs = [[0,4095//bin_y,0,4095//bin_x]]

    # method == fit :
    # for each region, fit the bias structure
    # It is time consuming to fit the whole 2D structure of the bias. 
    # Instead, I collapse it along the y-direction and fit it.
    # It's (much) faster, and accurate enough. 
    # method == row_med :
    # This is faster still, and as accurate. No fitting required, we just collapsed
    # the bias, and take the mean after rejecting outliers. Should be the default,
    # unless you know what you are doing ...

    out_data = numpy.zeros(numpy.shape(orig_data), dtype='d')
    for i in range(len(sci_regs)):
        reg = sci_regs[i]
        curr_data = orig_data[reg[0]:reg[1]+1,reg[2]:reg[3]+1]
    
        ny, nx = numpy.shape(curr_data)
        linx = numpy.arange(nx, dtype='d')
        liny = numpy.arange(ny, dtype='d')
        full_x, full_y = numpy.meshgrid(linx, liny)

        if method == 'row_med':
           #ny, nx = numpy.shape(curr_data)
           row_med = numpy.zeros(nx, dtype='d')
           # iterate over columns to excise CRs... bleh
           for i in range(nx):
               curr_col = curr_data[:,i]
               curr_med = numpy.median(curr_col)
               good_inds = numpy.nonzero(
                   numpy.abs(curr_col - curr_med) < 20.0)[0]
               row_med[i] = numpy.mean(curr_col[good_inds])
           #row_med = numpy.median(curr_data, axis=0)
           bias_sub = row_med**numpy.ones(
               numpy.shape(curr_data), dtype='d')
           # Fred's update (bias fit) ------
           # To remove the variations (some at least) along the 
           # y-direction, let's smooth the residuals !
           residual = curr_data - bias_sub
           residual[numpy.where(residual>150)]=0.0 # remove 'CR'
           residual_blur = ndimage.gaussian_filter(residual,sigma=[50,50])
           bias_sub += residual_blur
           # -----
           out_data[reg[0]:reg[1]+1,reg[2]:reg[3]+1] = bias_sub

        if method == 'fit':
        # Fit the bias with an appropriate model
        # Initial conditions differ depending on the camera

            if camera == 'WiFeSRed' :
                p0 = [0.,1., -200.,4100., 1., 20.,-800., 0.0001]
                fa = {'x': linx, 'z':curr_data.mean(axis=0), 
                      'err':curr_data.std(axis=0), 'camera':camera }
                constraints = [{'limited':[0,0]},
                               {'limited':[0,0]},
                               {'limited':[0,0]},            
                               {'limited':[1,0], 'limits':[numpy.max(full_x[0,:]),0]},
                               {'limited':[0,0]},
                               {'limited':[0,0]},
                               {'limited':[0,1], 'limits':[0,numpy.min(full_x[0,:])]},
                               {'limited':[0,0]}
                               ]

            else :
                p0 = [-3.,5., -500.,4100., 0.0001, 200.,-800., 0.0001, 1, 0.01,4000.]
                fa = {'x': linx, 'z':curr_data.mean(axis=0), 
                      'err':curr_data.std(axis=0), 'camera':camera }
                constraints = [{'limited':[0,0]},
                               {'limited':[0,0]},
                               {'limited':[0,0]},            
                               {'limited':[1,0], 'limits':[numpy.max(full_x[0,:]),0]},
                               {'limited':[0,0]},
                               {'limited':[0,0]},
                               {'limited':[0,1], 'limits':[0,numpy.min(full_x[0,:])]},
                               {'limited':[0,0]},
                               {'limited':[1,0], 'limits':[0,0]},
                               {'limited':[0,0]},
                               {'limited':[0,0]}
                               ]

            print(' Fitting bias frame %s' % bias_img.split('/')[-1])
            fit_result = mpfit.mpfit(error_wifes_bias_model, p0, functkw = fa, 
                                     parinfo = constraints, quiet=not verbose) 
            p1 = fit_result.params
            #print p1
            if fit_result.status <= 0 or fit_result.status == 5  :
                print(' Fit may have failed : mpfit status:',fit_result.status)
                print("I'll plot this one for sanity check...")
                plot = True
    
            out_data[reg[0]:reg[1]+1,reg[2]:reg[3]+1] = wifes_bias_model(p1,full_x,
                                                                         camera)
        # Plot for test purposes ...
        if plot:
            plt.figure()
            plt.plot(linx, curr_data.mean(axis=0), 'k-', label='raw bias', lw=2)
            if method == 'row_med':
                plt.plot(linx, out_data.mean(axis=0),'r-', label = 'row_med bias')
                plt.plot(linx, curr_data.mean(axis=0) - out_data.mean(axis=0), 'g-', label='residual')       
            if method == 'fit' :
                plt.plot(linx, curr_data.mean(axis=0) - wifes_bias_model(p1, linx, camera), 'g', label='residual')
                plt.plot(linx, wifes_bias_model(p1, linx, camera), 'r', label='model fit')

            plt.axhline(0, numpy.min(linx), numpy.max(linx), color='k')
            plt.xlabel('x [pixels]')
            plt.ylabel(' bias signal collapsed along y')
            plt.legend(loc='lower left', fancybox=True, shadow=True)
            plt.xlim([numpy.min(linx), numpy.max(linx)])
            plt.ylim([-10, 10])
            plt.title('Fitting bias frame %s' % bias_img.split('/')[-1])
            plt.show()
    
        # ----------------------------------------------------
    # save it!
    outfits = pyfits.HDUList(f1)
    outfits[data_hdu].data = out_data
    outfits[data_hdu].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f1.close()
    return

#------------------------------------------------------------------------
# WIFES specific tasks
def derive_slitlet_profiles(flatfield_fn,
                            output_fn,
                            data_hdu=0,
                            verbose=False,
                            shift_global=True,
                            plot=False,
                            bin_x=None,
                            bin_y=None):
    f = pyfits.open(flatfield_fn)
    flat_data = f[data_hdu].data
    orig_hdr = f[data_hdu].header
    f.close()
    # check if halfframe
    halfframe = is_halfframe(flatfield_fn)
    # check for binning, if not specified read from header
    try:
        bin_temp = orig_hdr['CCDSUM'].split()
        default_bin_x = int(float(bin_temp[0]))
        default_bin_y = int(float(bin_temp[1]))
    except:
        default_bin_x = 1
        default_bin_y = 1
    if bin_x == None:
        bin_x = default_bin_x
    if bin_y == None:
        bin_y = default_bin_y
    # first check which camera it is
    if orig_hdr['CAMERA'] == 'WiFeSRed':
        baseline_defs = red_slitlet_defs
    else:
        baseline_defs = blue_slitlet_defs
    # now fit slitlet profiles!!
    new_slitlet_defs = {}
    y_shift_vals = []
    if halfframe:
        nslits = 12
        offset = 4096//bin_y
    else:
        nslits = 25
        offset = 0
    for i in range(nslits):
        init_curr_defs = baseline_defs[str(i+1)]
        # convert to binned values
        curr_defs = [
            ((init_curr_defs[0]-1)//bin_x)+1,
            ((init_curr_defs[1]-1)//bin_x)+1,
            ((init_curr_defs[2]-1)//bin_y)+1,
            ((init_curr_defs[3]-1)//bin_y)+1]
        # check data outside these regions by 8/bin_y pixels
        y_buff = 20//bin_y
        mod_defs = [curr_defs[0]-1, curr_defs[1],
                    curr_defs[2]-1-y_buff-offset,
                    curr_defs[3]+y_buff-offset]
        expanded_data = flat_data[mod_defs[2]:mod_defs[3],
                                  mod_defs[0]:mod_defs[1]]
        #------------------
        # fit for best new center!
        init_yprof = numpy.sum(expanded_data, axis=1)
        y_prof = (init_yprof - init_yprof.min())/(
            init_yprof.max()-init_yprof.min())
        # center = halfway between edges where it drops below 0.001
        bright_inds = numpy.nonzero(y_prof>0.1)[0]
        new_ymin = bright_inds[0]-1
        new_ymax = bright_inds[-1]+1
        orig_ctr = 0.5*float(len(y_prof))
        new_ctr = 0.5*(new_ymin+new_ymax)
        if plot:
            import pylab
            pylab.figure()
            pylab.plot(y_prof, color='b')
            pylab.axvline(orig_ctr, color='r')
            pylab.axvline(new_ctr, color='g')
            pylab.show()
        # now adjust the slitlet definitions!
        y_shift = bin_y*int(new_ctr-orig_ctr)
        if verbose:
            print('Fitted shift of %d (unbinned) pixels for slitlet %d' % (
                y_shift, i+1))
        y_shift_vals.append(y_shift)
        final_defs = [init_curr_defs[0],
                      init_curr_defs[1],
                      init_curr_defs[2]+y_shift,
                      init_curr_defs[3]+y_shift]
        new_slitlet_defs[str(i+1)] = final_defs
    # finally, use a single global shift if requested
    if shift_global:
        best_shift = int(numpy.mean(numpy.array(y_shift_vals)))
        if verbose:
            print('Best global shift is %d (unbinned) pixels' % best_shift)
        final_slitlet_defs = {}
        for i in range(25):
            init_curr_defs = baseline_defs[str(i+1)]
            final_defs = [init_curr_defs[0],
                          init_curr_defs[1],
                          init_curr_defs[2]+best_shift,
                          init_curr_defs[3]+best_shift]
            final_slitlet_defs[str(i+1)] = final_defs
    else:
        final_slitlet_defs = new_slitlet_defs
    # save it!
    f3 = open(output_fn, 'wb')
    pickle.dump(final_slitlet_defs, f3)
    f3.close()
    return

# Fred's update (sag)
def interslice_cleanup(input_fn, output_fn,
                       slitlet_def_file=None,
                       bin_x=None, bin_y=None,
                       data_hdu=0,
                       offset = 0.4,
                       buffer = 0,
                       radius=10,
                       nsig_lim=3.0,
                       verbose=False,
                       plot=False,
                       savefigs=False,
                       save_prefix='cleanup_',
                       method='2D'):
    #------------------------------------
    # 1) Open the flat field
    f = pyfits.open(input_fn)
    header = f[data_hdu].header
    data = f[data_hdu].data
    f.close()
    # Check if half-frame
    halfframe = is_halfframe(input_fn)
    # check which channel (blue / red) it is!!
    camera = header['CAMERA']
    #------------------------------------
    # 2) Get the slitlets boundaries
    if slitlet_def_file != None:
        f2 = open(slitlet_def_file, 'rb')
        init_slitlet_defs = pickle.load(f2)
        f2.close()
    elif camera == 'WiFeSRed':
        init_slitlet_defs = red_slitlet_defs
    else:
        init_slitlet_defs = blue_slitlet_defs
    # NEW MJC CODE: ADD IN BUFFER
    slitlet_defs = {}
    for i in range(1,26):
        slit_num = str(i)
        init_slitdefs = init_slitlet_defs[slit_num]
        slitlet_defs[slit_num] = [
            init_slitdefs[0],
            init_slitdefs[1],
            max(init_slitdefs[2]-buffer,1),
            min(init_slitdefs[3]+buffer,4096)]
    # check for binning, if no specified read from header
    try:
        bin_temp = header['CCDSUM'].split()
        default_bin_x = int(float(bin_temp[0]))
        default_bin_y = int(float(bin_temp[1]))
    except:
        default_bin_x = 1
        default_bin_y = 1
    if bin_x == None:
        bin_x = default_bin_x
    if bin_y == None:
        bin_y = default_bin_y
    #------------------------------------
    # 3) Create temporary storage structure
    inter_smooth = numpy.zeros_like(data)
    fitted = numpy.zeros_like(data)
    if halfframe :
        nslits = 12
    else :
        nslits = 25
    slitlets_n = numpy.arange(26-nslits,26,1)
    #------------------------------------
    # 4) Get rid of the 'science-full region' and leave only the interslice
    # (and smooth it as well)
    for slit in slitlets_n:
        # Get the slit boundaries
        [xmin,xmax,ymin,ymax] = slitlet_defs[numpy.str(slit)]
        # Account for binning
        xmin = numpy.round(xmin//bin_x)
        xmax = numpy.round(xmax//bin_x)
        ymin = numpy.round(ymin//bin_y)
        ymax = numpy.round(ymax//bin_y)
        # Clear the appropriate region in temporary structure
        inter_smooth[ymin:ymax,xmin:xmax] = numpy.nan
        # Now, smooth the remaining interslice region
        if slit == 1:
            symin = ymax
            symax = numpy.shape(data)[0]
        else :
            symin = ymax
            symax = numpy.round(slitlet_defs[numpy.str(slit-1)][2]//bin_y)
        # Need to get rid of cosmic rays
        # here's a quick and dirty way of doing it ...
        # better idea anyone ?
        tmp = numpy.zeros_like(data[symin:symax,xmin:xmax])
        median = numpy.median(data[symin:symax,xmin:xmax])
        std = numpy.std(data[symin:symax,xmin:xmax])        
        tmp[data[symin:symax,xmin:xmax]< (median+nsig_lim*std)] = \
            data[symin:symax,xmin:xmax][data[symin:symax,xmin:xmax]< 
                                        (median+nsig_lim*std)]
        if method=='2D':
            inter_smooth[symin:symax,xmin:xmax] = \
                ndimage.gaussian_filter(tmp,
                                        sigma=[radius,radius])
        elif method=='1D':
            inter_smooth[symin:symax,xmin:xmax] += numpy.median(tmp)
        # And don't forget the extra slice
        if slit == 25 :
            symin=1
            symax=ymin
            if method=='2D':
                inter_smooth[symin:symax,xmin:xmax] = \
                    ndimage.gaussian_filter(data[symin:symax,xmin:xmax],
                                            sigma=[radius,radius])
            elif method=='1D':
                inter_smooth[symin:symax,xmin:xmax] += numpy.median(tmp)                    
    #------------------------------------
    # 5) Great, now we can interpolate this and reconstruct the contamination
    # Do each slitlet individually to avoid overloading the memory
    # Sampling
    dx = 10
    dy = 3
    for slit in slitlets_n:
        [xmin,xmax,y2,y3] = slitlet_defs[numpy.str(slit)]
        # Account for binning
        xmin = numpy.round(xmin//bin_x)
        xmax = numpy.round(xmax//bin_x)
        y2 = numpy.round(y2//bin_y)
        y3 = numpy.round(y3//bin_y)
        if slit == 1:
            y4 = numpy.shape(data)[0]
            y1 = numpy.round(slitlet_defs['2'][3]//bin_y)
        elif slit == 25 :
            y4 = numpy.round(slitlet_defs['24'][2]//bin_y)
            y1 = 1
        else:
            y4 = numpy.round(slitlet_defs[numpy.str(slit-1)][2]//bin_y)
            y1 = numpy.round(slitlet_defs[numpy.str(slit+1)][3]//bin_y)         
        # Select a subsample of point to do the integration
        x = numpy.arange(xmin+3,xmax-3,dx)
        y = numpy.append(numpy.arange(y1+1,y2-1,dy),
                         numpy.arange(y3+1,y4-1,dy))
        grid = numpy.zeros((len(y),len(x)))
        for i in range(len(x)):
            for j in range(len(y)):
                grid[j,i] = inter_smooth[y[j],x[i]]
        #------------------------------------
        # 6) Actually perform the interpolation
        func = interp.RectBivariateSpline(y,x,grid, kx=1,ky=1)
        # Note : because of the large gap to fill, kx and ky have little effect
        # reconstruct missing slice
        xall = numpy.arange(xmin,xmax,1)
        yall = numpy.arange(y1,y4,1)
        fitted[y1:y4,xmin:xmax] = func(yall,xall)
    #------------------------------------
    # 7) All done ! Let's save it all ...
    f = pyfits.open(input_fn)
    f[0].data = fitted
    f.writeto(input_fn[:-5]+'_corr.fits',overwrite=True) # Save the corretion image separately just in case
    f[0].data = data - fitted + offset*numpy.mean(fitted)
    f[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    f.writeto(output_fn,overwrite=True)
    f.close()
    if verbose:    
        print(' Additive offset:',offset*numpy.mean(fitted))
    # 8) Plot anything ?
    if plot or savefigs:
        myvmax = numpy.max(fitted)
        #------------------
        fig = pylab.figure()
        pylab.imshow(data, vmin=0,vmax=myvmax,cmap='nipy_spectral', origin='lower')
        pylab.title('Pre-corrected '+input_fn.split('/')[-1])
        if savefigs:
            save_fn = save_prefix+'flat_orig_data.png'
            pylab.savefig(save_fn)
        #------------------
        pylab.figure()
        pylab.imshow(fitted,vmin=0,vmax=myvmax,cmap='nipy_spectral', origin='lower')
        pylab.title('Fitted contamination for '+input_fn.split('/')[-1])
        if savefigs:
            save_fn = save_prefix+'flat_glow_data.png'
            pylab.savefig(save_fn)
        #------------------            
        pylab.figure()
        pylab.imshow(data-fitted+offset*numpy.mean(fitted),vmin=0,
                     vmax=myvmax,cmap='nipy_spectral', origin='lower')
        pylab.title('Corrected '+output_fn.split('/')[-1])
        if savefigs:
            save_fn = save_prefix+'flat_sub_data.png'
            pylab.savefig(save_fn)
        if plot:
            pylab.show()
    return

def wifes_slitlet_mef(inimg, outimg, data_hdu=0,
                      bin_x=None, bin_y=None,
                      slitlet_def_file=None):
    f = pyfits.open(inimg)
    #outfits = pyfits.HDUList([f[0]])
    outfits = pyfits.HDUList([pyfits.PrimaryHDU(header=f[0].header)])
    old_hdr = f[data_hdu].header
    full_data = f[data_hdu].data
    f.close()
    # check if halfframe
    halfframe = is_halfframe(inimg)
    # check which channel (blue / red) it is!!
    camera = old_hdr['CAMERA']
    rdnoise = old_hdr['RDNOISE']
    # get slitlet definitions!
    # new ones if defined, otherwise use baseline values!
    if slitlet_def_file != None:
        f2 = open(slitlet_def_file, 'rb')
        try:
            slitlet_defs = pickle.load(f2, fix_imports=True, encoding='latin')
        except:
            slitlet_defs = pickle.load(f2) # for python 2.7
        f2.close()
    elif camera == 'WiFeSRed':
        slitlet_defs = red_slitlet_defs
    else:
        slitlet_defs = blue_slitlet_defs
    # check for binning, if no specified read from header
    try:
        bin_temp = old_hdr['CCDSUM'].split()
        default_bin_x = int(float(bin_temp[0]))
        default_bin_y = int(float(bin_temp[1]))
    except:
        default_bin_x = 1
        default_bin_y = 1
    if bin_x == None:
        bin_x = default_bin_x
    if bin_y == None:
        bin_y = default_bin_y
    #---------------------------
    # ** NEW VERSION 0.5.8 - CREATE DQ IMAGE **
    dq_img = numpy.zeros(numpy.shape(full_data))
    # flag bad pixels on new detectors
    epoch = determine_detector_epoch(inimg)
    if int(float(epoch[1])) > 3:
        if camera == 'WiFeSRed':
            if halfframe:
                y_bad_min = (3242-2048)//bin_y-1
                y_bad_max = 2048//bin_y
            else:
                y_bad_min = 3282//bin_y-1
                y_bad_max = 4096//bin_y
            x_bad_1 = 774//bin_x-1
            x_bad_2 = 775//bin_x-1
            dq_img[y_bad_min:y_bad_max,x_bad_1] = 1
            dq_img[y_bad_min:y_bad_max,x_bad_2] = 1
        else:
            x_bad = 1529//bin_x-1
            if halfframe:
                y_bad_min = 0
                y_bad_max = 2048//bin_y
            else:
                y_bad_min = 752//bin_y-1
                y_bad_max = 4096//bin_y
            dq_img[y_bad_min:y_bad_max,x_bad] = 1
    #---------------------------
    # for each slitlet, save it to a single header extension
    if halfframe:
        nslits = 12
        #offset = 4096//bin_y
        old_ymin = int(float(
            old_hdr['CCDSEC'].split(',')[1].split(':')[0]))-1
        offset = old_ymin//bin_y
    else:
        nslits = 25
        offset = 0
    for i in range(25):
        init_curr_defs = slitlet_defs[str(i+1)]
        # convert to binned values
        curr_defs = [
            ((init_curr_defs[0]-1)//bin_x)+1,
            ((init_curr_defs[1]-1)//bin_x)+1,
            ((init_curr_defs[2]-1)//bin_y)+1,
            ((init_curr_defs[3]-1)//bin_y)+1]
        # horrible kluge to make sure everything has the same dimensions!
        if (curr_defs[1]-curr_defs[0]+1) != (4096//bin_x):
            curr_defs[1] -= 1
        if (curr_defs[3]-curr_defs[2]+1) != (86//bin_y):
            curr_defs[3] -= 1
        # get the data
        dim_str = '[%d:%d,%d:%d]' % (
            curr_defs[0], curr_defs[1],
            curr_defs[2], curr_defs[3])
        mod_defs = [curr_defs[0]-1, curr_defs[1],
                    curr_defs[2]-1-offset, curr_defs[3]-offset]
        #print curr_defs, (curr_defs[3]-curr_defs[2]), init_curr_defs, (
        #    init_curr_defs[3]-init_curr_defs[2])
        if halfframe and ((i+1) > nslits):
            new_data = numpy.zeros([86//bin_y, 4096//bin_x], dtype='d')
        else:
            new_data = full_data[mod_defs[2]:mod_defs[3],
                                 mod_defs[0]:mod_defs[1]]
        # create fits hdu
        hdu_name = 'SCI%d' % (i+1)
        new_hdu = pyfits.ImageHDU(new_data, old_hdr, name=hdu_name)
        new_hdu.header.set('CCDSEC',  dim_str)
        new_hdu.header.set('DATASEC', dim_str)
        new_hdu.header.set('TRIMSEC', dim_str)
        outfits.append(new_hdu)
        gc.collect()
    #print squirrel
    # VARIANCE EXTENSIONS
    for i in range(25):
        init_curr_defs = slitlet_defs[str(i+1)]
        # convert to binned values
        curr_defs = [
            ((init_curr_defs[0]-1)//bin_x)+1,
            ((init_curr_defs[1]-1)//bin_x)+1,
            ((init_curr_defs[2]-1)//bin_y)+1,
            ((init_curr_defs[3]-1)//bin_y)+1]
        # horrible kluge to make sure everything has the same dimensions!
        if (curr_defs[1]-curr_defs[0]+1) != (4096//bin_x):
            curr_defs[1] -= 1
        if (curr_defs[3]-curr_defs[2]+1) != (86//bin_y):
            curr_defs[3] -= 1
        #print curr_defs
        #print squirrel
        # get the data
        dim_str = '[%d:%d,%d:%d]' % (
            curr_defs[0], curr_defs[1],
            curr_defs[2], curr_defs[3])
        mod_defs = [curr_defs[0]-1, curr_defs[1],
                    curr_defs[2]-1-offset, curr_defs[3]-offset]
        if halfframe and ((i+1) > nslits):
            new_data = numpy.zeros([86//bin_y, 4096//bin_x], dtype='d')
        else:
            new_data = full_data[mod_defs[2]:mod_defs[3],
                                 mod_defs[0]:mod_defs[1]]
        var_data = new_data+rdnoise**2
        # create fits hdu
        hdu_name = 'VAR%d' % (i+1)
        new_hdu = pyfits.ImageHDU(var_data, old_hdr, name=hdu_name)
        new_hdu.header.set('CCDSEC',  dim_str)
        new_hdu.header.set('DATASEC', dim_str)
        new_hdu.header.set('TRIMSEC', dim_str)
        outfits.append(new_hdu)
        gc.collect()
    # DATA QUALITY EXTENSIONS
    for i in range(25):
        init_curr_defs = slitlet_defs[str(i+1)]
        # convert to binned values
        curr_defs = [
            ((init_curr_defs[0]-1)//bin_x)+1,
            ((init_curr_defs[1]-1)//bin_x)+1,
            ((init_curr_defs[2]-1)//bin_y)+1,
            ((init_curr_defs[3]-1)//bin_y)+1]
        # horrible kluge to make sure everything has the same dimensions!
        if (curr_defs[1]-curr_defs[0]+1) != (4096//bin_x):
            curr_defs[1] -= 1
        if (curr_defs[3]-curr_defs[2]+1) != (86//bin_y):
            curr_defs[3] -= 1
        # get the data
        dim_str = '[%d:%d,%d:%d]' % (
            curr_defs[0], curr_defs[1],
            curr_defs[2], curr_defs[3])
        mod_defs = [curr_defs[0]-1, curr_defs[1],
                    curr_defs[2]-1-offset, curr_defs[3]-offset]
        ny = curr_defs[3]-curr_defs[2]+1
        nx = curr_defs[1]-curr_defs[0]+1
        if halfframe and ((i+1) > nslits):
            dq_data = numpy.ones([ny,nx])
        else:
            dq_data = dq_img[mod_defs[2]:mod_defs[3],
                             mod_defs[0]:mod_defs[1]]
        # OLD METHOD
        #if halfframe and ((i+1) > nslits):
        #    dq_data = numpy.ones([ny,nx])
        #else:
        #    dq_data = numpy.zeros([ny,nx])
        # create fits hdu
        hdu_name = 'DQ%d' % (i+1)
        new_hdu = pyfits.ImageHDU(dq_data, old_hdr, name=hdu_name)
        new_hdu.header.set('CCDSEC',  dim_str)
        new_hdu.header.set('DATASEC', dim_str)
        new_hdu.header.set('TRIMSEC', dim_str)
        outfits.append(new_hdu)
        gc.collect()
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    return

def wifes_slitlet_mef_ns(inimg, outimg_obj, outimg_sky,
                         data_hdu=0, bin_x=None, bin_y=None, nod_dy=80,
                         slitlet_def_file=None):
    f = pyfits.open(inimg)
    outfits_obj = pyfits.HDUList([pyfits.PrimaryHDU(header=f[0].header)])
    outfits_sky = pyfits.HDUList([pyfits.PrimaryHDU(header=f[0].header)])
    old_hdr = f[data_hdu].header
    full_data = f[data_hdu].data
    f.close()
    rdnoise = old_hdr['RDNOISE']
    # check which channel (blue / red) it is!!
    camera = old_hdr['CAMERA']
    # get slitlet definitions!
    # new ones if defined, otherwise use baseline values!
    if slitlet_def_file != None:
        f2 = open(slitlet_def_file, 'rb')
        slitlet_defs = pickle.load(f2, fix_imports=True, encoding='latin')
        f2.close()
    elif camera == 'WiFeSRed':
        slitlet_defs = red_slitlet_defs
    else:
        slitlet_defs = blue_slitlet_defs
    # check for binning, if no specified read from header
    try:
        bin_temp = old_hdr['CCDSUM'].split()
        default_bin_x = int(float(bin_temp[0]))
        default_bin_y = int(float(bin_temp[1]))
    except:
        default_bin_x = 1
        default_bin_y = 1
    if bin_x == None:
        bin_x = default_bin_x
    if bin_y == None:
        bin_y = default_bin_y
    #------------------------------------
    # for each slitlet, save it to a single header extension
    for i in range(nslits):
        init_curr_defs = slitlet_defs[str(i+1)]
        #------------------
        # convert to binned values
        obj_defs = [
            ((init_curr_defs[0]-1)//bin_x)+1,
            ((init_curr_defs[1]-1)//bin_x)+1,
            ((init_curr_defs[2]-1)//bin_y)+1,
            ((init_curr_defs[3]-1)//bin_y)+1]
        sky_defs = [
            ((init_curr_defs[0]-1)//bin_x)+1,
            ((init_curr_defs[1]-1)//bin_x)+1,
            ((init_curr_defs[2]-1+nod_dy)//bin_y)+1,
            ((init_curr_defs[3]-1+nod_dy)//bin_y)+1]
        #------------------
        # horrible kluge to make sure everything has the same dimensions!
        if (obj_defs[1]-obj_defs[0]+1) != (4096//bin_x):
            obj_defs[1] -= 1
        if (obj_defs[3]-obj_defs[2]+1) != (86//bin_y):
            obj_defs[3] -= 1
        if (sky_defs[1]-sky_defs[0]+1) != (4096//bin_x):
            sky_defs[1] -= 1
        if (sky_defs[3]-sky_defs[2]+1) != (86//bin_y):
            sky_defs[3] -= (sky_defs[3]-sky_defs[2]+1) - (86//bin_y)
            #sky_defs[3] -= 1
        #------------------
        # get the object data
        obj_dim_str = '[%d:%d,%d:%d]' % (
            obj_defs[0], obj_defs[1],
            obj_defs[2], obj_defs[3])
        obj_mod_defs = [obj_defs[0]-1, obj_defs[1],
                        obj_defs[2]-1, obj_defs[3]]
        obj_data = full_data[obj_mod_defs[2]:obj_mod_defs[3],
                             obj_mod_defs[0]:obj_mod_defs[1]]
        # kill outer 3//bin_y pixels!!
        ykill = 4//bin_y
        obj_data[:ykill,:]  *= 0.0
        obj_data[-ykill:,:] *= 0.0
        # create fits hdu for object
        hdu_name = 'SCI%d' % (i+1)
        obj_hdu = pyfits.ImageHDU(obj_data, old_hdr, name=hdu_name)
        obj_hdu.header.set('CCDSEC',  obj_dim_str)
        obj_hdu.header.set('DATASEC', obj_dim_str)
        obj_hdu.header.set('TRIMSEC', obj_dim_str)
        # fix the exposure time!!
        exptime_true = float(old_hdr['SEXP'])*float(old_hdr['NSUBEXPS'])
        obj_hdu.header.set('EXPTIME', exptime_true,
                              comment='Total NS exposure time')
        outfits_obj.append(obj_hdu)
        #------------------
        # and the sky data
        sky_dim_str = '[%d:%d,%d:%d]' % (
            sky_defs[0], sky_defs[1],
            sky_defs[2], sky_defs[3])
        sky_mod_defs = [sky_defs[0]-1, sky_defs[1],
                        sky_defs[2]-1, sky_defs[3]]
        nsx = sky_mod_defs[1]-sky_mod_defs[0]
        nsy = sky_mod_defs[3]-sky_mod_defs[2]
        # horrible fix to include bad NS regions definition for slitlet 1
        sky_data = numpy.zeros([nsy, nsx])
        true_sky_data = full_data[sky_mod_defs[2]:sky_mod_defs[3],
                                  sky_mod_defs[0]:sky_mod_defs[1]]
        tsy, tsx = numpy.shape(true_sky_data)
        sky_data[:tsy,:tsx] = true_sky_data
        # kill outer 3//bin_y pixels!!
        ykill = 4//bin_y
        sky_data[:ykill,:]  *= 0.0
        sky_data[-ykill:,:] *= 0.0
        # create fits hdu for sky
        hdu_name = 'SCI%d' % (i+1)
        sky_hdu = pyfits.ImageHDU(sky_data, old_hdr, name=hdu_name)
        sky_hdu.header.set('CCDSEC',  sky_dim_str)
        sky_hdu.header.set('DATASEC', sky_dim_str)
        sky_hdu.header.set('TRIMSEC', sky_dim_str)
        # fix the exposure time!!
        exptime_true = float(old_hdr['NEXP'])*float(old_hdr['NSUBEXPS'])
        sky_hdu.header.set('EXPTIME', exptime_true,
                              comment='Total NS exposure time')
        outfits_sky.append(sky_hdu)
        gc.collect()
    #------------------------------------
    #print squirrel
    # VARIANCE EXTENSIONS
    for i in range(nslits):
        init_curr_defs = slitlet_defs[str(i+1)]
        #------------------
        # convert to binned values
        obj_defs = [
            ((init_curr_defs[0]-1)//bin_x)+1,
            ((init_curr_defs[1]-1)//bin_x)+1,
            ((init_curr_defs[2]-1)//bin_y)+1,
            ((init_curr_defs[3]-1)//bin_y)+1]
        sky_defs = [
            ((init_curr_defs[0]-1)//bin_x)+1,
            ((init_curr_defs[1]-1)//bin_x)+1,
            ((init_curr_defs[2]-1+nod_dy)//bin_y)+1,
            ((init_curr_defs[3]-1+nod_dy)//bin_y)+1]
        #------------------
        # horrible kluge to make sure everything has the same dimensions!
        if (obj_defs[1]-obj_defs[0]+1) != (4096//bin_x):
            obj_defs[1] -= 1
        if (obj_defs[3]-obj_defs[2]+1) != (86//bin_y):
            obj_defs[3] -= 1
        if (sky_defs[1]-sky_defs[0]+1) != (4096//bin_x):
            sky_defs[1] -= 1
        if (sky_defs[3]-sky_defs[2]+1) != (86//bin_y):
            sky_defs[3] -= (sky_defs[3]-sky_defs[2]+1) - (86//bin_y)
            #sky_defs[3] -= 1
        #------------------
        # get the object data
        obj_dim_str = '[%d:%d,%d:%d]' % (
            obj_defs[0], obj_defs[1],
            obj_defs[2], obj_defs[3])
        obj_mod_defs = [obj_defs[0]-1, obj_defs[1],
                        obj_defs[2]-1, obj_defs[3]]
        obj_data = full_data[obj_mod_defs[2]:obj_mod_defs[3],
                             obj_mod_defs[0]:obj_mod_defs[1]]
        # kill outer 3//bin_y pixels!!
        ykill = 4//bin_y
        obj_data[:ykill,:]  *= 0.0
        obj_data[-ykill:,:] *= 0.0
        obj_var = obj_data+rdnoise**2
        # create fits hdu
        hdu_name = 'VAR%d' % (i+1)
        obj_hdu = pyfits.ImageHDU(obj_var, old_hdr, name=hdu_name)
        obj_hdu.header.set('CCDSEC',  obj_dim_str)
        obj_hdu.header.set('DATASEC', obj_dim_str)
        obj_hdu.header.set('TRIMSEC', obj_dim_str)
        # fix the exposure time!!
        exptime_true = float(old_hdr['SEXP'])*float(old_hdr['NSUBEXPS'])
        obj_hdu.header.set('EXPTIME', exptime_true,
                              comment='Total NS exposure time')
        outfits_obj.append(obj_hdu)
        #------------------
        # get the sky data
        sky_dim_str = '[%d:%d,%d:%d]' % (
            sky_defs[0], sky_defs[1],
            sky_defs[2], sky_defs[3])
        sky_mod_defs = [sky_defs[0]-1, sky_defs[1],
                        sky_defs[2]-1, sky_defs[3]]
        nsx = sky_mod_defs[1]-sky_mod_defs[0]
        nsy = sky_mod_defs[3]-sky_mod_defs[2]
        # horrible fix to include bad NS regions definition for slitlet 1
        sky_data = numpy.zeros([nsy, nsx])
        true_sky_data = full_data[sky_mod_defs[2]:sky_mod_defs[3],
                                  sky_mod_defs[0]:sky_mod_defs[1]]
        tsy, tsx = numpy.shape(true_sky_data)
        sky_data[:tsy,:tsx] = true_sky_data
        # kill outer 3//bin_y pixels!!
        ykill = 4//bin_y
        sky_data[:ykill,:]  *= 0.0
        sky_data[-ykill:,:] *= 0.0
        sky_var = sky_data+rdnoise**2
        # create fits hdu
        hdu_name = 'VAR%d' % (i+1)
        sky_hdu = pyfits.ImageHDU(sky_var, old_hdr, name=hdu_name)
        sky_hdu.header.set('CCDSEC',  sky_dim_str)
        sky_hdu.header.set('DATASEC', sky_dim_str)
        sky_hdu.header.set('TRIMSEC', sky_dim_str)
        # fix the exposure time!!
        exptime_true = float(old_hdr['NEXP'])*float(old_hdr['NSUBEXPS'])
        sky_hdu.header.set('EXPTIME', exptime_true,
                              comment='Total NS exposure time')
        outfits_sky.append(sky_hdu)
        gc.collect()
    #------------------------------------
    # DATA QUALITY EXTENSIONS
    for i in range(nslits):
        init_curr_defs = slitlet_defs[str(i+1)]
        #------------------
        # convert to binned values
        obj_defs = [
            ((init_curr_defs[0]-1)//bin_x)+1,
            ((init_curr_defs[1]-1)//bin_x)+1,
            ((init_curr_defs[2]-1)//bin_y)+1,
            ((init_curr_defs[3]-1)//bin_y)+1]
        sky_defs = [
            ((init_curr_defs[0]-1)//bin_x)+1,
            ((init_curr_defs[1]-1)//bin_x)+1,
            ((init_curr_defs[2]-1+nod_dy)//bin_y)+1,
            ((init_curr_defs[3]-1+nod_dy)//bin_y)+1]
        #------------------
        # horrible kluge to make sure everything has the same dimensions!
        if (obj_defs[1]-obj_defs[0]+1) != (4096//bin_x):
            obj_defs[1] -= 1
        if (obj_defs[3]-obj_defs[2]+1) != (86//bin_y):
            obj_defs[3] -= 1
        if (sky_defs[1]-sky_defs[0]+1) != (4096//bin_x):
            sky_defs[1] -= 1
        if (sky_defs[3]-sky_defs[2]+1) != (86//bin_y):
            sky_defs[3] -= (sky_defs[3]-sky_defs[2]+1) - (86//bin_y)
            #sky_defs[3] -= 1
        #------------------
        # get the data
        obj_dim_str = '[%d:%d,%d:%d]' % (
            obj_defs[0], obj_defs[1],
            obj_defs[2], obj_defs[3])
        sky_dim_str = '[%d:%d,%d:%d]' % (
            sky_defs[0], sky_defs[1],
            sky_defs[2], sky_defs[3])
        nyo = obj_defs[3]-obj_defs[2]+1
        nxo = obj_defs[1]-obj_defs[0]+1
        nys = sky_defs[3]-sky_defs[2]+1
        nxs = sky_defs[1]-sky_defs[0]+1
        obj_dq = numpy.zeros([nyo,nxo])
        sky_dq = numpy.zeros([nys,nxs])
        #------------------
        # create fits hdu for object
        hdu_name = 'DQ%d' % (i+1)
        obj_hdu = pyfits.ImageHDU(obj_dq, old_hdr, name=hdu_name)
        obj_hdu.header.set('CCDSEC',  obj_dim_str)
        obj_hdu.header.set('DATASEC', obj_dim_str)
        obj_hdu.header.set('TRIMSEC', obj_dim_str)
        # fix the exposure time!!
        exptime_true = float(old_hdr['SEXP'])*float(old_hdr['NSUBEXPS'])
        obj_hdu.header.set('EXPTIME', exptime_true,
                              comment='Total NS exposure time')
        outfits_obj.append(obj_hdu)
        #------------------
        # create fits hdu for object
        hdu_name = 'DQ%d' % (i+1)
        sky_hdu = pyfits.ImageHDU(sky_dq, old_hdr, name=hdu_name)
        sky_hdu.header.set('CCDSEC',  sky_dim_str)
        sky_hdu.header.set('DATASEC', sky_dim_str)
        sky_hdu.header.set('TRIMSEC', sky_dim_str)
        # fix the exposure time!!
        exptime_true = float(old_hdr['NEXP'])*float(old_hdr['NSUBEXPS'])
        sky_hdu.header.set('EXPTIME', exptime_true,
                              comment='Total NS exposure time')
        outfits_sky.append(sky_hdu)
        gc.collect()
    #------------------------------------
    outfits_obj[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits_sky[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits_obj.writeto(outimg_obj, overwrite=True)
    outfits_sky.writeto(outimg_sky, overwrite=True)
    return

#------------------------------------------------------------------------
def wifes_response_pixel(inimg, outimg, wsol_fn = None):
    # check if halfframe
    halfframe = is_halfframe(inimg)
    if halfframe:
        nslits = 12
    else:
        nslits = 25
    # now open and operate on data
    f = pyfits.open(inimg)
    outfits = pyfits.HDUList(f)
    for i in range(nslits):
        curr_hdu = i+1
        orig_data = f[curr_hdu].data
        # NEW 2012-04-20
        # rectify data!
        if wsol_fn != None:
            f3 = pyfits.open(wsol_fn)
            wave = f3[i+1].data
            f3.close()
            print('Transforming data for Slitlet %d' % curr_hdu)
            rect_data =  transform_data(orig_data, wave)
            curr_ff_rowwise_ave = numpy.median(rect_data, axis=0)
            curr_ff_illum = numpy.median(rect_data/curr_ff_rowwise_ave,
                                         axis=1)
            curr_model = ((numpy.ones(numpy.shape(rect_data))
                           *curr_ff_rowwise_ave).T*curr_ff_illum).T
            orig_model = detransform_data(
                curr_model, orig_data, wave)
            normed_data = orig_data / orig_model
        else:
            curr_ff_rowwise_ave = numpy.median(orig_data, axis=0)
            curr_ff_illum = numpy.median(orig_data/curr_ff_rowwise_ave,
                                         axis=1)
            curr_model = ((numpy.ones(numpy.shape(orig_data))
                           *curr_ff_rowwise_ave).T*curr_ff_illum).T
            normed_data = orig_data / curr_model
        outfits[curr_hdu].data = normed_data
        # need to fit this for each slitlet
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f.close()
    return

def wifes_response_poly(inimg, outimg,
                        wsol_fn=None,
                        zero_var=True,
                        polydeg=7):
    # check if halfframe
    halfframe = is_halfframe(inimg)
    if halfframe:
        nslits = 12
        mid_slit = 7
    else:
        nslits = 25
        mid_slit = 13
    # now open and operate on data
    f = pyfits.open(inimg)
    outfits = pyfits.HDUList(f)
    #------------------------------------
    # fit a smooth polynomial to the middle slice
    midslice_data = f[mid_slit].data
    if wsol_fn != None:
        f3 = pyfits.open(wsol_fn)
        wave = f3[mid_slit].data
        f3.close()
        rect_data, lam_array =  transform_data(
            midslice_data, wave, return_lambda=True)
    else:
        rect_data = midslice_data
        lam_array = numpy.arange(len(rect_data[0,:]),dtype='d')
    # fit polynomial to median data
    curr_ff_rowwise_ave = numpy.median(rect_data, axis=0)
    curr_y = numpy.log10(curr_ff_rowwise_ave)
    good_inds = numpy.nonzero((curr_y == curr_y)*
                              (curr_ff_rowwise_ave > 0.0))[0]
    # add points far away to force edge derivatives to be preserved
    next_x = lam_array[good_inds][5:-5]
    next_y = curr_y[good_inds][5:-5]
    yderiv = (next_y[1:]-next_y[:-1])/(next_x[1:]-next_x[:-1])
    nave_lo = 50
    nave_hi = 100
    nextend_lo = 150.0
    nextend_hi = 300.0
    ydave_lo = numpy.mean(yderiv[:nave_lo])
    ydave_hi = numpy.mean(yderiv[-nave_hi:])
    dxlo = numpy.arange(1,nextend_lo+1,dtype='d')*(next_x[1]-next_x[0])
    new_xlo = next_x[0] - dxlo
    new_ylo = next_y[0] - dxlo*ydave_lo
    dxhi = numpy.arange(1,nextend_hi+1,dtype='d')*(next_x[-1]-next_x[-2])
    new_xhi = next_x[-1] + dxhi
    new_yhi = next_y[-1] + dxhi*ydave_hi
    fit_x = numpy.concatenate((
        new_xlo,
        next_x,
        new_xhi))
    fit_y = numpy.concatenate((
        new_ylo,
        next_y,
        new_yhi))
    smooth_poly = numpy.polyfit(fit_x, fit_y, polydeg)
    #------------------------------------
    # divide rectified data by the smooth polynomial
    for i in range(nslits):
        curr_hdu = i+1
        orig_data = f[curr_hdu].data
        # NEW 2012-04-20
        # rectify data!
        if wsol_fn != None:
            f3 = pyfits.open(wsol_fn)
            wave = f3[i+1].data
            f3.close()
            print('Transforming data for Slitlet %d' % curr_hdu)
            rect_data, lam_array =  transform_data(
                orig_data, wave, return_lambda=True)
            curr_norm_array = 10.0**(numpy.polyval(smooth_poly,
                                                   lam_array))
            init_normed_data = rect_data / curr_norm_array
            normed_data = detransform_data(
                init_normed_data, orig_data, wave)
        else:
            lam_array = numpy.arange(len(orig_data[0,:]),dtype='d')
            curr_norm_array = 10.0**(numpy.polyval(smooth_poly,
                                                   lam_array))
            normed_data = rect_data / curr_norm_array
        outfits[curr_hdu].data = normed_data
        if zero_var:
            var_hdu = curr_hdu+25
            outfits[var_hdu].data *= 0.0
        # need to fit this for each slitlet
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f.close()
    return

def wifes_2dim_response(spec_inimg,
                        spatial_inimg,
                        outimg,
                        wsol_fn=None,
                        zero_var=True,
                        plot=False,
                        savefigs=False,
                        save_prefix='response_',
                        polydeg=7,
                        resp_min=1.0e-6):
    # check if halfframe
    halfframe = is_halfframe(spec_inimg)
    if halfframe:
        nslits = 12
        mid_slit = 7
    else:
        nslits = 25
        mid_slit = 13
    # open the two files!
    f1 = pyfits.open(spec_inimg)
    f2 = pyfits.open(spatial_inimg)
    ndy, ndx = numpy.shape(f1[1].data)
    xarr = numpy.arange(ndx)
    yarr = numpy.arange(ndy)
    full_x, full_y = numpy.meshgrid(xarr, yarr)
    all_ypos_flat = full_y.flatten()
    outfits = pyfits.HDUList(f1)
    #------------------------------------
    # SPECTRAL FLAT
    # fit a smooth polynomial to the middle slice
    midslice_data = f1[mid_slit].data
    if wsol_fn != None:
        f3 = pyfits.open(wsol_fn)
        wave = f3[mid_slit].data
        f3.close()
        rect_data, mid_lam_array =  transform_data(
            midslice_data, wave, return_lambda=True)
    else:
        rect_data = midslice_data
        mid_lam_array = numpy.arange(len(rect_data[0,:]),dtype='d')
    out_y_full, out_lambda_full = numpy.meshgrid(yarr, mid_lam_array)
    disp_ave = abs(numpy.mean(mid_lam_array[1:]-mid_lam_array[:-1]))
    # fit polynomial to median data
    curr_ff_rowwise_ave = numpy.median(rect_data, axis=0)
    curr_y = numpy.log10(curr_ff_rowwise_ave)
    good_inds = numpy.nonzero((curr_y == curr_y)*
                              (curr_ff_rowwise_ave > 0.0))[0]
    # add points far away to force edge derivatives to be preserved
    if pyfits.getval(spec_inimg, 'CAMERA') == 'WiFeSRed':
        next_x = mid_lam_array[good_inds][500:-10]
        next_y = curr_y[good_inds][500:-10]
    else:
        next_x = mid_lam_array[good_inds][50:-100]
        next_y = curr_y[good_inds][50:-100]
    yderiv = (next_y[1:]-next_y[:-1])/(next_x[1:]-next_x[:-1])
    nave_lo = 50
    nave_hi = 100
    nextend_lo = 150.0
    nextend_hi = 300.0
    ydave_lo = numpy.mean(yderiv[:nave_lo])
    ydave_hi = numpy.mean(yderiv[-nave_hi:])
    dxlo = numpy.arange(1,nextend_lo+1,dtype='d')*(next_x[1]-next_x[0])
    new_xlo = next_x[0] - dxlo
    new_ylo = next_y[0] - dxlo*ydave_lo
    dxhi = numpy.arange(1,nextend_hi+1,dtype='d')*(next_x[-1]-next_x[-2])
    new_xhi = next_x[-1] + dxhi
    new_yhi = next_y[-1] + dxhi*ydave_hi
    fit_x = numpy.concatenate((
        new_xlo,
        next_x,
        new_xhi))
    fit_y = numpy.concatenate((
        new_ylo,
        next_y,
        new_yhi))
    smooth_poly = numpy.polyfit(fit_x, fit_y, polydeg)
    #------------------------------------
    # SPATIAL FLAT
    # get median spatial flat spectrum
    midslice_data = f2[mid_slit].data
    if wsol_fn != None:
        f3 = pyfits.open(wsol_fn)
        wave = f3[mid_slit].data
        f3.close()
        rect_data =  transform_data(
            midslice_data, wave)
    else:
        rect_data = midslice_data
    # fit polynomial to median data
    spatial_flat_spec = numpy.median(rect_data, axis=0)
    spec_norm = curr_ff_rowwise_ave / (
        10.0**(numpy.polyval(smooth_poly, mid_lam_array)))
    #import pylab
    #pylab.figure()
    #pylab.plot(mid_lam_array, spec_norm)
    #pylab.show()
    spat_interp = scipy.interpolate.interp1d(
        mid_lam_array, spatial_flat_spec/spec_norm,
        bounds_error=False,
        fill_value=0.0)
    #------------------------------------
    #------------------------------------
    illum = numpy.zeros([ndy, 25], dtype='d')
    # divide rectified data by the smooth polynomial
    for i in range(nslits):
        curr_hdu = i+1
        orig_spec_data = f1[curr_hdu].data
        orig_spat_data = f2[curr_hdu].data
        # NEW 2012-04-20
        # rectify data!
        if wsol_fn != None:
            f3 = pyfits.open(wsol_fn)
            wave = f3[i+1].data
            f3.close()
            # convert the *x* pixels to lambda
            dw = numpy.abs(wave[:,1:]-wave[:,:-1])
            full_dw = numpy.zeros(numpy.shape(wave))
            full_dw[:,1:] = dw
            full_dw[:,0] = dw[:,0]
            print('Transforming data for Slitlet %d' % curr_hdu)
            # SPECTRAL FLAT
            rect_spec_data, lam_array =  transform_data(
                orig_spec_data, wave, return_lambda=True)
            curr_norm_array = 10.0**(numpy.polyval(smooth_poly,
                                                   lam_array))
            init_normed_data = rect_spec_data / curr_norm_array
            xstart = int(0.25*len(lam_array))
            xstop  = int(0.75*len(lam_array))
            norm_region = init_normed_data[:,xstart:xstop]
            next_normed_data = (init_normed_data.T/
                                numpy.median(norm_region,axis=1)).T
            # SPATIAL FLAT
            rect_spat_data =  transform_data(
                orig_spat_data, wave, return_lambda=False)
            alt_dw = lam_array[1]-lam_array[0]
            alt_flat_spec = spat_interp(lam_array)
            curr_interp_spat = numpy.zeros(numpy.shape(out_y_full),
                                           dtype='d')
            alt_interp_spat = numpy.zeros(numpy.shape(rect_spat_data),
                                          dtype='d')
            #import pylab
            #f1 = pylab.figure()
            #sp1 = f1.add_subplot(111)
            #pylab.title(str(curr_hdu))
            #f2 = pylab.figure()
            #sp2 = f2.add_subplot(111)
            for q in range(ndy):
                curr_x = wave[q,:]
                curr_y = orig_spat_data[q,:]/full_dw[q,:]
                sort_order = curr_x.argsort()
                curr_interp = scipy.interpolate.interp1d(
                    curr_x[sort_order], curr_y[sort_order],
                    bounds_error=False,
                    fill_value=0.0)
                new_data = disp_ave*curr_interp(mid_lam_array)
                curr_interp_spat[:,q] = new_data
                #norm_func = (new_data / spatial_flat_spec)
                alt_y = rect_spat_data[q,:]
                norm_func = (alt_y / (
                    alt_flat_spec*next_normed_data[q,:]))
                alt_interp_spat[q,:] = norm_func
                #print len(new_data)
                #print len(norm_func)
                #norm_func = (alt_y / (
                #    alt_flat_spec))
                #sp1.plot(norm_func[500:-200])
                #sp2.plot(next_normed_data[q,:])
            #pylab.show()
            #spat_ratio = curr_interp_spat.T / spatial_flat_spec
            spat_ratio = alt_interp_spat
            spat_ratio[numpy.nonzero(spat_ratio!=spat_ratio)] = 0.0
            spat_ratio[numpy.nonzero(spat_ratio<0.0)] = 0.0
            spat_flat = numpy.median(spat_ratio[:,xstart:xstop], axis=1)
            # transform back
            final_normed_data = (next_normed_data.T*spat_flat).T
            normed_data = detransform_data(
                final_normed_data, orig_spec_data, wave)
            normed_data[numpy.nonzero(normed_data <= resp_min)] = resp_min
        else:
            lam_array = numpy.arange(len(orig_spec_data[0,:]),dtype='d')
            curr_norm_array = 10.0**(numpy.polyval(smooth_poly,
                                                   lam_array))
            normed_data = orig_spec_data / curr_norm_array
        outfits[curr_hdu].data = normed_data
        illum[:,i] = numpy.sum(normed_data, axis=1)
        if zero_var:
            var_hdu = curr_hdu+25
            outfits[var_hdu].data *= 0.0
        # need to fit this for each slitlet
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f1.close()
    f2.close()
    #---------------------------------------------
    # diagnostic plots!
    if plot or savefigs:
        #------------------
        # (1) spectral fit
        pylab.figure()
        pylab.plot(mid_lam_array, curr_ff_rowwise_ave, 'b')
        pylab.plot(mid_lam_array,
                   10.0**(numpy.polyval(smooth_poly, mid_lam_array)),
                   color='g', lw=3)
        if savefigs:
            save_fn = save_prefix+'flat_spectrum.png'
            pylab.savefig(save_fn)
        #------------------
        # (2) illumination correction
        pylab.figure()
        pylab.imshow(illum,
                     interpolation='nearest',
                     origin='lower',
                     cmap=pylab.cm.Greys_r)
        if savefigs:
            save_fn = save_prefix+'flat_illumination.png'
            pylab.savefig(save_fn)
        #------------------
        if plot:
            pylab.show()
    #---------------------------------------------
    return

def wifes_illumination(spatial_inimg,
                       outimg,
                       wsol_fn=None,
                       zero_var=True,
                       polydeg=7,
                       resp_min=1.0e-6):
    # check if halfframe
    halfframe = is_halfframe(spatial_inimg)
    if halfframe:
        nslits = 12
        mid_slit = 7
    else:
        nslits = 25
        mid_slit = 13
    # open the two files!
    f2 = pyfits.open(spatial_inimg)
    ndy, ndx = numpy.shape(f2[1].data)
    flat_ones = numpy.ones([ndy, ndx], dtype='d')
    xarr = numpy.arange(ndx)
    yarr = numpy.arange(ndy)
    full_x, full_y = numpy.meshgrid(xarr, yarr)
    all_ypos_flat = full_y.flatten()
    outfits = pyfits.HDUList(f2)
    #------------------------------------
    # SPATIAL FLAT
    # get median spatial flat spectrum
    midslice_data = f2[mid_slit].data
    if wsol_fn != None:
        f3 = pyfits.open(wsol_fn)
        wave = f3[mid_slit].data
        f3.close()
        rect_data, mid_lam_array =  transform_data(
            midslice_data, wave, return_lambda=True)
    else:
        rect_data = midslice_data
        mid_lam_array = numpy.arange(len(rect_data[0,:]),dtype='d')
    out_y_full, out_lambda_full = numpy.meshgrid(yarr, mid_lam_array)
    # derive median data
    spatial_flat_spec = numpy.median(rect_data, axis=0)
    #------------------------------------
    #------------------------------------
    # divide rectified data by the smooth polynomial
    for i in range(nslits):
        curr_hdu = i+1
        orig_spat_data = f2[curr_hdu].data
        # NEW 2012-04-20
        # rectify data!
        if wsol_fn != None:
            f3 = pyfits.open(wsol_fn)
            wave = f3[i+1].data
            f3.close()
            print('Transforming data for Slitlet %d' % curr_hdu)
            # SPATIAL FLAT
            wave_flat = wave.flatten()
            dw = wave[:,1:]-wave[:,:-1]
            full_dw = numpy.zeros(numpy.shape(wave))
            full_dw[:,1:] = dw
            full_dw[:,0] = dw[:,0]
            wdisp = numpy.mean(dw)
            new_flux = orig_spat_data
            curr_flux_flat = (new_flux/full_dw).flatten()
            curr_interp_spat = wdisp*scipy.interpolate.griddata(
                (wave_flat, all_ypos_flat),
                curr_flux_flat,
                (out_lambda_full, out_y_full),
                method='linear',
                fill_value=0.0)
            spat_ratio = curr_interp_spat.T / spatial_flat_spec
            spat_ratio[numpy.nonzero(spat_ratio!=spat_ratio)] = 0.0
            spat_ratio[numpy.nonzero(spat_ratio<0.0)] = 0.0
            spat_flat = numpy.median(spat_ratio, axis=1)
        else:
            spat_flat = (numpy.sum(orig_spat_data, axis=1)
                         /numpy.sum(spatial_flat_spec))
        spat_flat[numpy.nonzero(spat_flat < resp_min)[0]] = resp_min    
        normed_data = (flat_ones.T*spat_flat).T
        outfits[curr_hdu].data = normed_data
        if zero_var:
            var_hdu = curr_hdu+25
            outfits[var_hdu].data *= 0.0
        # need to fit this for each slitlet
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f1.close()
    f2.close()
    return

#------------------------------------------------------------------------
# function to fit the wavelength solution!
from wifes_wsol import derive_wifes_wave_solution
from wifes_wsol import derive_wifes_skyline_solution

#------------------------------------------------------------------------
# function to fit the wire solution!
def derive_wifes_wire_solution(inimg, out_file,
                               bin_x=None, bin_y=None,
                               #fit_zones=[10,30,45,65],
                               fit_zones=[16,26,54,70],
                               flux_threshold=20,
                               wire_polydeg=1,
                               xlims='default'):
    # note: later have these parameters as user-input kwargs
    # SOME PARAMTERS FOR THE FITTING
    init_nave = 10
    bg_polydeg = 3
    ctr_polydeg = 15
    xlim_defaults = {'U7000' : [1,2500],
                     'B7000' : [1,4096],
                     'R7000' : [1000,3000],
                     'I7000' : [1,4096],
                     'B3000' : [1,2600],
                     'R3000' : [850,4096]}
    # define region for fitting light profile
    init_fit_pmin_1 = fit_zones[0]
    init_fit_pmax_1 = fit_zones[1]
    init_fit_pmin_2 = fit_zones[2]
    init_fit_pmax_2 = fit_zones[3]
    #------------------------------------
    f = pyfits.open(inimg)
    # check if it is halfframe
    halfframe = is_halfframe(inimg)
    if halfframe:
        nslits = 12
    else:
        nslits = 25    
    # figure out which channel it is
    if f[1].header['CAMERA'] == 'WiFeSRed':
        channel = 'red'
        grating = f[1].header['GRATINGR']
    else:
        channel = 'blue'    
        grating = f[1].header['GRATINGB']
    # figure out the binning!
    try:
        bin_temp = f[1].header['CCDSUM'].split()
        default_bin_x = int(float(bin_temp[0]))
        default_bin_y = int(float(bin_temp[1]))
    except:
        default_bin_x = 1
        default_bin_y = 1
    if bin_x == None:
        bin_x = default_bin_x
    if bin_y == None:
        bin_y = default_bin_y
    # get the xmax for fitting if default selected
    if xlims == 'default':
        xmin = (xlim_defaults[grating][0]-1)//bin_x
        xmax = (xlim_defaults[grating][1]-1)//bin_x
    else:
        xmin = xlims[0]
        xmax = xlims[1]
    # adjust the fit regions accordingly
    nave = init_nave//bin_x
    fit_pmin_1 = init_fit_pmin_1//bin_y
    fit_pmax_1 = init_fit_pmax_1//bin_y
    fit_pmin_2 = init_fit_pmin_2//bin_y
    fit_pmax_2 = init_fit_pmax_2//bin_y
    #------------------------------------
    # get data size and set up output data array
    temp_data = f[1].data
    ny, nx = numpy.shape(temp_data)
    ctr_results = numpy.zeros([25,nx],dtype='f')
    ccd_x = numpy.arange(nx,dtype='f')
    # number of groupings to do
    ng = nx//nave - 1
    for q in range(nslits):
        slit_ind = q+1
        #print slit_ind
        test_data = f[slit_ind].data
        nr, junk = numpy.shape(test_data)
        fit_inds = numpy.nonzero(
            ((numpy.arange(nr) >= fit_pmin_1)*
             (numpy.arange(nr) <= fit_pmax_1))+
            ((numpy.arange(nr) >= fit_pmin_2)*
             (numpy.arange(nr) <= fit_pmax_2)))[0]
        x_full = numpy.arange(nr,dtype='d')
        x_fit = x_full[fit_inds]
        ##initializing at zero on the following line means that the frame_fit_y[-1] works in the first group
        ##but will be excluded by the xlims later, as long as the xmin is always > 0
        frame_fit_x = [0.0]
        frame_fit_y = [0.0]
        for i in range(ng):
            # get median profile
            yprof = numpy.median(
                test_data[:,numpy.max([0,nave*i-nave//2]):numpy.min([nave*i+nave//2,test_data.shape[1]])+1],axis=1)
            # if there is no usable data, use previous value!
            if numpy.max(yprof) < flux_threshold:
                wire_ctr = frame_fit_y[-1]
            # otherwise fit the region outside of ~[30:50]
            else:
                curr_fit = numpy.polyfit(
                    x_fit,
                    numpy.log10(yprof[fit_inds]),
                    bg_polydeg)
                curr_fvals = 10.0**(numpy.polyval(curr_fit, x_full))
                # now take residuals and calculate a simple centroid
                wire_x = numpy.arange(fit_pmin_1,fit_pmax_2)
                wire_y = (curr_fvals-yprof)[fit_pmin_1:fit_pmax_2]
                wire_ctr = single_centroid_prof_fit(wire_y,x=wire_x)
            frame_fit_x.append(float(nave*i))
            frame_fit_y.append(wire_ctr)
        fit_x_arr = numpy.array(frame_fit_x)
        fit_y_arr = numpy.array(frame_fit_y)
        good_inds = numpy.nonzero((fit_x_arr==fit_x_arr)*
                                  (fit_y_arr==fit_y_arr)*
                                  (fit_x_arr>=xmin)*
                                  (fit_x_arr<=xmax))[0]
        wire_trend = numpy.polyfit(fit_x_arr[good_inds],
                                   fit_y_arr[good_inds],
                                   wire_polydeg)
        trend_y = numpy.polyval(wire_trend, ccd_x)
        #pdb.set_trace()
        #pylab.figure()
        #pylab.plot(fit_x_arr, fit_y_arr, 'b')
        #pylab.plot(ccd_x, trend_y, 'r')
        #pylab.show()
        ctr_results[q,:] = trend_y
    f.close()
    results = pyfits.PrimaryHDU(data = ctr_results)
    g = pyfits.HDUList([results])
    g[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    g.writeto(out_file, overwrite=True)
    return

#-------------------------------------------------------------
# DATA CUBE!!!
def generate_wifes_cube(inimg, outimg,
                        wire_fn,
                        wsol_fn,
                        wmin_set = None, wmax_set = None,dw_set = None,
                        bin_x=None, bin_y=None,
                        ny_orig=76,
                        offset_orig=4,
                        multithread=False,
                        verbose=True,
                        adr=False):
    if multithread:
        generate_wifes_cube_multithread(
            inimg, outimg,
            wire_fn,
            wsol_fn,
            wmin_set = wmin_set, wmax_set = wmax_set,dw_set=dw_set,
            bin_x=bin_x, bin_y=bin_y,
            ny_orig=ny_orig,
            offset_orig=offset_orig,
            verbose=verbose,
            adr=adr)
    else:
        generate_wifes_cube_oneproc(
            inimg, outimg,
            wire_fn,
            wsol_fn,
            wmin_set = wmin_set, wmax_set = wmax_set, dw_set = dw_set,
            bin_x=bin_x, bin_y=bin_y,
            ny_orig=ny_orig,
            offset_orig=offset_orig,
            verbose=verbose,
            adr=adr)

# ------------------- Fred's update (2) -----------------------

from wifes_adr import ha_degrees
from wifes_adr import dec_dms2dd
from wifes_adr import adr_x_y

# ------------------- Fred's update (2) -----------------------
# Now takes into account ADR correction
# 3D interpolation is hard/slow/impossible. Because of the symmetry of the system

def generate_wifes_cube_oneproc(
    inimg, outimg,
    wire_fn,
    wsol_fn,
    wmin_set = None, wmax_set = None, dw_set = None,
    bin_x=None, bin_y=None,
    ny_orig=76,
    offset_orig=4,
    multithread=False,
    verbose=True,
    adr=False):
    #---------------------------
    # check if halfframe
    halfframe = is_halfframe(inimg)
    if halfframe:
        nslits = 12
    else:
        nslits = 25
    # setup base x/y array
    f3 = pyfits.open(inimg)
    ndy, ndx = numpy.shape(f3[1].data)
    xarr = numpy.arange(ndx)
    yarr = numpy.arange(ndy)
    full_x, full_y = numpy.meshgrid(xarr, yarr)
    obs_hdr = f3[0].header
    # figure out the binning!
    try:
        bin_temp = f3[1].header['CCDSUM'].split()
        default_bin_x = int(float(bin_temp[0]))
        default_bin_y = int(float(bin_temp[1]))
    except:
        default_bin_x = 1
        default_bin_y = 1
    if bin_x == None:
        bin_x = default_bin_x
    if bin_y == None:
        bin_y = default_bin_y
    # get the min/max wavelengths across all slits, and average dispersion
    frame_wmin = 0.0
    frame_wmax = 20000.0
    frame_wdisps = []
    for i in range(1,nslits+1):
        f4 = pyfits.open(wsol_fn)
        wave = f4[i].data
        f4.close()
        curr_wmin = numpy.max(numpy.min(wave,axis=1))
        curr_wmax = numpy.min(numpy.max(wave,axis=1))
        curr_wdisp = numpy.abs(numpy.mean(wave[:,1:]-wave[:,:-1]))
        if curr_wmin > frame_wmin:
            frame_wmin = curr_wmin
        if curr_wmax < frame_wmax:
            frame_wmax = curr_wmax
        frame_wdisps.append(curr_wdisp)
        #print numpy.shape(f3[i].data)
    if dw_set != None:
        disp_ave = dw_set 
    else : 
        disp_ave = numpy.mean(frame_wdisps)
    if verbose:
        print(' Data spectral resolution (min/max):',numpy.round(numpy.min(frame_wdisps),4),numpy.round(numpy.max(frame_wdisps),4))
        print(' Cube spectral resolution : ',disp_ave)
    # excise lowest pixel so interpolation doesn't fail
    frame_wmin += disp_ave
    # finally check against the user input value
    if wmin_set != None:
        final_frame_wmin = max(wmin_set, frame_wmin)
    else:
        final_frame_wmin = frame_wmin
    if wmax_set != None:
        final_frame_wmax = min(wmax_set, frame_wmax)
    else:
        final_frame_wmax = frame_wmax
    out_lambda = numpy.arange(final_frame_wmin, final_frame_wmax, disp_ave)
    # set up output data
    # load in spatial solutions
    try:
        f5 = pyfits.open(wire_fn)
        wire_trans = f5[0].data
        f5.close()
    except:
        wire_trans = numpy.zeros([ndy, ndx], dtype='d')+numpy.max(yarr)/2
    #ny=70//bin_y+1
    wire_offset = float(offset_orig)/float(bin_y)
    ny=ny_orig//bin_y
    nx=25
    nlam=len(out_lambda)
    # for each slitlet...
    init_out_y = numpy.arange(ny,dtype='d')
    out_y = init_out_y - numpy.median(init_out_y)
    out_y_full, out_lambda_full = numpy.meshgrid(out_y, out_lambda)
    #---------------------------
    # Prepare ADR corrections ...
    if adr :
        # observatory stuff
        lat = obs_hdr['LAT-OBS'] # degrees
        alt = obs_hdr['ALT-OBS'] # meters
        dec = dec_dms2dd(obs_hdr['DEC'])
        # want to calculate average HA...
        ha_start = obs_hdr['HA']
        ha_end   = obs_hdr['HAEND']
        ha = 0.5*(ha_degrees(ha_start)+
                  ha_degrees(ha_end))
        # and average ZD...
        zd_start = obs_hdr['ZD']
        zd_end   = obs_hdr['ZDEND']
        zd = numpy.radians(0.5*(zd_start+zd_end))
        secz = 1.0/numpy.cos(zd)
        tanz = numpy.tan(zd)
        # telescope PA!
        telpa = numpy.radians(obs_hdr['TELPAN'])
        # THIS SHOULD BE A FIXED VALUE!!
        #adr_ref = adr_x_y(numpy.array([out_lambda[-1]]),secz, ha, dec, lat, 
        #                  teltemp = 0.0, telpres=700.0, telpa=telpa)
        adr_ref = adr_x_y(numpy.array([5600.0]),
                          secz, ha, dec, lat, 
                          teltemp = 0.0, telpres=700.0, telpa=telpa)
    #---------------------------
    outfits = pyfits.HDUList(f3)
    # Create a temporary storage array for first iteration
    flux_data_cube_tmp = numpy.zeros([nx,ny,nlam])
    var_data_cube_tmp = numpy.ones([nx,ny,nlam])
    dq_data_cube_tmp = numpy.ones([nx,ny,nlam])
    orig_flux        = numpy.zeros([25,43,4096])
    # First interpolation : Wavelength + y (=wire & ADR)
    if verbose:
        print(' -> Step 1: interpolating along lambda and y (2D interp.)\r')
        sys.stdout.write('\r 0%')
        sys.stdout.flush()
    for i in range(1,nslits+1):
        f4 = pyfits.open(wsol_fn)
        wave = f4[i].data
        f4.close()
        #print ' -> generating cube data for slitlet %d' % i
        curr_flux = f3[i].data
        curr_var  = f3[25+i].data
        curr_dq   = f3[50+i].data
        # convert the *x* pixels to lambda
        dw = numpy.abs(wave[:,1:]-wave[:,:-1])
        full_dw = numpy.zeros(numpy.shape(wave))
        full_dw[:,1:] = dw
        full_dw[:,0] = dw[:,0]
        # and *y* to real y
        curr_wire = wire_trans[i-1,:]
        all_ypos = full_y - curr_wire - wire_offset # what is the wire offset doing ?
        # from the y-lambda-flux data, interpolate the flux
        # for the desired output y-lambda grid
        wave_flat = wave.flatten()
        all_ypos_flat = all_ypos.flatten()
        curr_flux_flat = (curr_flux/full_dw).flatten()
        curr_var_flat  = (curr_var/full_dw**2).flatten()
        curr_dq_flat   = curr_dq.flatten()
        # Calculate the ADR corrections (this is slow)
        if adr :
            adr = adr_x_y(wave_flat,secz, ha, dec, lat, teltemp = 0.0, telpres=700.0, 
                          telpa=telpa)
            adr_y = adr[1]-adr_ref[1]
            all_ypos_flat -= adr_y
        # Does the interpolation (this is equally slow ...)
        flux_data_cube_tmp[i-1,:,:] = disp_ave * scipy.interpolate.griddata(
            (wave_flat, all_ypos_flat),
            curr_flux_flat,
            (out_lambda_full, out_y_full),
            method='linear',
            fill_value=0.0).T
        ##two images that can be checked in testing
        #orig_flux[i-1,:,:] = curr_flux
        #ospat = numpy.sum(orig_flux,axis=2)
        #spat = numpy.sum(flux_data_cube_tmp,axis=2)
        var_data_cube_tmp[i-1,:,:] = (disp_ave**2) * scipy.interpolate.griddata(
            (wave_flat, all_ypos_flat),
            curr_var_flat,
            (out_lambda_full, out_y_full),
            method='linear',
            fill_value=0.0).T
        dq_data_cube_tmp[i-1,:,:] = scipy.interpolate.griddata(
            (wave_flat, all_ypos_flat),
            curr_dq_flat,
            (out_lambda_full, out_y_full),
            method='nearest',
            fill_value=1).T
        if verbose:
            sys.stdout.flush()
            sys.stdout.write('\r\r %d' % (i/(float(nslits))*100.) + '%')
            sys.stdout.flush()
            if i == nslits : sys.stdout.write('\n')
    # Second interpolation : x (=ADR)
    if adr :
        # To avoid interpolation issues at the edges,
        # add two extra values on either side (0 in this version).
        in_x = numpy.arange(-1,nslits+1,1, dtype='d') 
        out_x = numpy.arange(nslits, dtype='d')
        if verbose:
            print(' -> Step 2: interpolating along x (1D interp.)')
        for i in range(0,nlam) :
            adr = adr_x_y(numpy.array([out_lambda[i]]),
                          secz,ha,dec,lat, teltemp = 0.0, 
                          telpres=700.0, 
                          telpa=telpa)
            adr_x = adr[0] - adr_ref[0]
            for j in range(0,ny) :
                # here is the actual appending of 0 flux
                this_flux = numpy.append(
                    numpy.append(
                    flux_data_cube_tmp[0,j,i],
                    flux_data_cube_tmp[:nslits,j,i]),
                    flux_data_cube_tmp[nslits-1,j,i])
                this_var = numpy.append(
                    numpy.append(
                    var_data_cube_tmp[0,j,i],
                    var_data_cube_tmp[:nslits,j,i]),
                    var_data_cube_tmp[nslits-1,j,i])
                this_dq = numpy.append(
                    numpy.append(
                    dq_data_cube_tmp[0,j,i],
                    dq_data_cube_tmp[:nslits,j,i]),
                    dq_data_cube_tmp[nslits-1,j,i])
                # do interpolation
                f = scipy.interpolate.interp1d(
                    in_x-adr_x, this_flux, kind='linear', 
                    fill_value=0.0, bounds_error=False)
                g = scipy.interpolate.interp1d(
                    in_x-adr_x, this_var, kind='linear', 
                    fill_value=0.0, bounds_error=False)
                h = scipy.interpolate.interp1d(
                    in_x-adr_x, this_dq, kind='nearest', 
                    fill_value=3, bounds_error=False)
                flux_data_cube_tmp[:nslits,j,i] = f(out_x)
                var_data_cube_tmp[:nslits,j,i] = g(out_x)
                dq_data_cube_tmp[:nslits,j,i] = h(out_x)
            if verbose:
                if i > 0 :
                    sys.stdout.flush()
                sys.stdout.write('\r\r %d' % (i/(nlam-1.)*100.) + '%')
                if i == nlam-1 :
                    sys.stdout.write('\n')
    # All done, at last ! Now, let's save it all ...
    # ALWAYS ITERATE TO 25 EVEN IF HALFFRAME HERE!!
    for i in range(1,26):
        # save to data cube
        outfits[i].data = flux_data_cube_tmp[i-1,:,:]
        outfits[i].header.set('CRVAL1', final_frame_wmin)
        outfits[i].header.set('CDELT1', disp_ave)
        outfits[i+25].data = var_data_cube_tmp[i-1,:,:]
        outfits[i+50].data = dq_data_cube_tmp[i-1,:,:]
        #outfits[i].header.set('NAXIS1', len(out_lambda))
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f3.close()
    return 

###--------------------------------------------------------------
# multithread cube processes
def save_wifes_cube_single_thread(
    wave_in, ypos_in, flux_in,
    wave_out, ypos_out, disp,
    out_fn, method, hdr, fill_value=0):
    interp_flux = disp*scipy.interpolate.griddata(
        (wave_in, ypos_in),
        flux_in,
        (wave_out, ypos_out),
        method=method,
        fill_value=fill_value)
    # save to a fits file
    outfits = pyfits.HDUList([pyfits.PrimaryHDU(data = interp_flux,
                                                header = hdr)])
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(out_fn, overwrite=True)
    return

def generate_wifes_cube_multithread(
    inimg, outimg,
    wire_fn,
    wsol_fn,
    wmin_set = None, wmax_set = None,dw_set = None,
    bin_x=None, bin_y=None,
    ny_orig=76,
    offset_orig=4,
    verbose=True,
    adr=False):
    #---------------------------
    # check if halfframe
    halfframe = is_halfframe(inimg)
    if halfframe:
        nslits = 12
    else:
        nslits = 25
    # setup base x/y array
    f3 = pyfits.open(inimg)
    ndy, ndx = numpy.shape(f3[1].data)
    xarr = numpy.arange(ndx)
    yarr = numpy.arange(ndy)
    full_x, full_y = numpy.meshgrid(xarr, yarr)
    obs_hdr = f3[1].header
    # figure out the binning!
    try:
        bin_temp = f3[1].header['CCDSUM'].split()
        default_bin_x = int(float(bin_temp[0]))
        default_bin_y = int(float(bin_temp[1]))
    except:
        default_bin_x = 1
        default_bin_y = 1
    if bin_x == None:
        bin_x = default_bin_x
    if bin_y == None:
        bin_y = default_bin_y
    # get the min/max wavelengths across all slits, and average dispersion
    frame_wmin = 0.0
    frame_wmax = 20000.0
    frame_wdisps = []
    for i in range(1,nslits+1):
        f4 = pyfits.open(wsol_fn)
        wave = f4[i].data
        f4.close()
        curr_wmin = numpy.max(numpy.min(wave,axis=1))
        curr_wmax = numpy.min(numpy.max(wave,axis=1))
        curr_wdisp = numpy.abs(numpy.mean(wave[:,1:]-wave[:,:-1]))
        if curr_wmin > frame_wmin:
            frame_wmin = curr_wmin
        if curr_wmax < frame_wmax:
            frame_wmax = curr_wmax
        frame_wdisps.append(curr_wdisp)
        #print numpy.shape(f3[i].data)
    if dw_set != None:
        disp_ave = dw_set
    else:
        disp_ave = numpy.mean(frame_wdisps)
    frame_wmin += disp_ave
    # finally check against the user input value
    if wmin_set != None:
        final_frame_wmin = max(wmin_set, frame_wmin)
    else:
        final_frame_wmin = frame_wmin
    if wmax_set != None:
        final_frame_wmax = min(wmax_set, frame_wmax)
    else:
        final_frame_wmax = frame_wmax
    
    if verbose:
        print(' Data spectral resolution (min/max):',numpy.round(numpy.min(frame_wdisps),4),numpy.round(numpy.max(frame_wdisps),4))
        print(' Cube spectral resolution : ',disp_ave)
    
    out_lambda = numpy.arange(final_frame_wmin, final_frame_wmax, disp_ave)
    # set up output data
    # load in spatial solutions
    # Rajika's edit. Allows for reduction without wires
    try:
        f5 = pyfits.open(wire_fn)
        wire_trans = f5[0].data
        f5.close()
    except:
        wire_trans = numpy.zeros([ndy, ndx], dtype='d')+numpy.max(yarr)/2
    #ny=70//bin_y+1
    wire_offset = float(offset_orig)/float(bin_y)
    ny=ny_orig//bin_y
    nx=25
    nlam=len(out_lambda)
    # for each slitlet...
    init_out_y = numpy.arange(ny,dtype='d')
    out_y = init_out_y - numpy.median(init_out_y)
    out_y_full, out_lambda_full = numpy.meshgrid(out_y, out_lambda)
    #print numpy.shape(out_y_full)
    #print numpy.shape(full_y)
    #print squirrel
    #---------------------------
    # Prepare ADR corrections ...
    if adr :
        # observatory stuff
        lat = obs_hdr['LAT-OBS'] # degrees
        alt = obs_hdr['ALT-OBS'] # meters
        dec = dec_dms2dd(obs_hdr['DEC'])
        # want to calculate average HA...
        ha_start = obs_hdr['HA']
        ha_end   = obs_hdr['HAEND']
        ha = 0.5*(ha_degrees(ha_start)+
                  ha_degrees(ha_end))
        # and average ZD...
        zd_start = obs_hdr['ZD']
        zd_end   = obs_hdr['ZDEND']
        zd = numpy.radians(0.5*(zd_start+zd_end))
        secz = 1.0/numpy.cos(zd)
        tanz = numpy.tan(zd)
        # telescope PA!
        telpa = numpy.radians(obs_hdr['TELPAN'])
        # THIS MUST BE A FIXED VALUE
        #adr_ref = adr_x_y(numpy.array([out_lambda[-1]]),secz, ha, dec, lat, 
        #                  teltemp = 0.0, telpres=700.0, telpa=telpa)
        adr_ref = adr_x_y(numpy.array([5600.0]),
                          secz, ha, dec, lat, 
                          teltemp = 0.0, telpres=700.0, telpa=telpa)
    #---------------------------
    if verbose:
        print(' -> Step 1: interpolating along lambda and y (2D interp.) MULTITHREAD\r')
        sys.stdout.write('\r 0%')
        sys.stdout.flush()
    threads = []

    for i in range(1,nslits+1):
        f4 = pyfits.open(wsol_fn)
        wave = f4[i].data
        wave_hdr = f4[i].header
        f4.close()
        #print 'Generating cube data for slitlet %d' % i
        curr_flux = f3[i].data
        curr_var  = f3[25+i].data
        curr_dq   = f3[50+i].data
        # convert the *x* pixels to lambda
        dw = numpy.abs(wave[:,1:]-wave[:,:-1])
        full_dw = numpy.zeros(numpy.shape(wave))
        full_dw[:,1:] = dw
        full_dw[:,0] = dw[:,0]
        # and *y* to real y
        curr_wire = wire_trans[i-1,:]
        all_ypos = full_y - curr_wire - wire_offset
        #print curr_wire
        #if i == 2:
        #    print squirrel
        # from the y-lambda-flux data, interpolate the flux
        # for the desired output y-lambda grid
        wave_flat = wave.flatten()
        all_ypos_flat = all_ypos.flatten()
    
        # Calculate the ADR corrections (this is slow)
        if adr :
            adr = adr_x_y(wave_flat,secz, ha, dec, lat, teltemp = 0.0,
                          telpres=700.0, telpa=telpa)
            adr_y = adr[1]-adr_ref[1]
            all_ypos_flat -= adr_y
        curr_flux_flat = (curr_flux/full_dw).flatten()
        curr_var_flat  = (curr_var/full_dw**2).flatten()
        curr_dq_flat   = curr_dq.flatten()
        # HERE DO THE THREADING!!
        temp_data_fn = 'tmp_cubegen_wifes_s%02d.fits' % (i)
        new_data_thread = multiprocessing.Process(
            target=save_wifes_cube_single_thread,
            kwargs={'wave_in' : wave_flat,
                    'ypos_in' : all_ypos_flat,
                    'flux_in' : curr_flux_flat,
                    'wave_out' : out_lambda_full,
                    'ypos_out' : out_y_full,
                    'disp' : disp_ave,
                    'out_fn' : temp_data_fn,
                    'method' : 'linear',
                    'hdr' : wave_hdr})
        new_data_thread.start()
        threads.append(new_data_thread)
        # and for var...
        temp_var_fn  = 'tmp_cubegen_wifes_s%02d.fits' % (i+25)
        new_var_thread = multiprocessing.Process(
            target=save_wifes_cube_single_thread,
            kwargs={'wave_in' : wave_flat,
                    'ypos_in' : all_ypos_flat,
                    'flux_in' : curr_var_flat,
                    'wave_out' : out_lambda_full,
                    'ypos_out' : out_y_full,
                    'disp' : disp_ave**2,
                    'out_fn' : temp_var_fn,
                    'method' : 'linear',
                    'hdr' : wave_hdr})
        new_var_thread.start()
        threads.append(new_var_thread)
        # and for DQ
        temp_dq_fn   = 'tmp_cubegen_wifes_s%02d.fits' % (i+50)
        new_dq_thread = multiprocessing.Process(
            target=save_wifes_cube_single_thread,
            kwargs={'wave_in' : wave_flat,
                    'ypos_in' : all_ypos_flat,
                    'flux_in' : curr_dq_flat,
                    'wave_out' : out_lambda_full,
                    'ypos_out' : out_y_full,
                    'disp' : 1.0,
                    'out_fn' : temp_dq_fn,
                    'method' : 'nearest',
                    'fill_value' : 1,
                    'hdr' : wave_hdr})
        new_dq_thread.start()
        threads.append(new_dq_thread)

    for i in range(len(threads)):
        t = threads[i]
        t.join()
        if verbose:
            sys.stdout.flush()
            sys.stdout.write('\r\r %d' % ((i+1)/(float(len(threads)))*100.) + '%')
            sys.stdout.flush()
            if i == (len(threads)-1) : sys.stdout.write('\n')

    #---------------------------
    # GATHER ALL TRANSFORMED DATA
    # Create a temporary storage array for first iteration
    flux_data_cube_tmp = numpy.zeros([nx,ny,nlam])
    var_data_cube_tmp = numpy.ones([nx,ny,nlam])
    dq_data_cube_tmp = numpy.ones([nx,ny,nlam])
    for i in range(1,nslits+1):
        temp_data_fn = 'tmp_cubegen_wifes_s%02d.fits' % (i)
        fm = pyfits.open(temp_data_fn)
        curr_interp_flux = fm[0].data
        flux_data_cube_tmp[i-1,:,:] = curr_interp_flux.T
        fm.close()
        temp_var_fn  = 'tmp_cubegen_wifes_s%02d.fits' % (i+25)
        fm = pyfits.open(temp_var_fn)
        curr_interp_var = fm[0].data
        var_data_cube_tmp[i-1,:,:] = curr_interp_var.T
        fm.close()
        temp_dq_fn   = 'tmp_cubegen_wifes_s%02d.fits' % (i+50)
        fm = pyfits.open(temp_dq_fn)
        curr_interp_dq = fm[0].data
        dq_data_cube_tmp[i-1,:,:] = curr_interp_dq.T
        fm.close()
        subprocess.call(['rm', '-f', temp_data_fn])
        subprocess.call(['rm', '-f', temp_var_fn])
        subprocess.call(['rm', '-f', temp_dq_fn])
    #---------------------------
    # Second interpolation : x (=ADR)
    if adr :
        # To avoid interpolation issues at the edges,
        # add two extra values on either side (0 in this version).
        in_x = numpy.arange(-1,nslits+1,1, dtype='d') 
        out_x = numpy.arange(nslits, dtype='d')
        if verbose:
            print(' -> Step 2: interpolating along x (1D interp.)')
        for i in range(0,nlam) :
            adr = adr_x_y(numpy.array([out_lambda[i]]),
                          secz,ha,dec,lat, teltemp = 0.0, 
                          telpres=700.0, 
                          telpa=telpa)
            adr_x = adr[0] - adr_ref[0]
            for j in range(0,ny) :
                # here is the actual appending of 0 flux
                this_flux = numpy.append(
                    numpy.append(
                    flux_data_cube_tmp[0,j,i],
                    flux_data_cube_tmp[:nslits,j,i]),
                    flux_data_cube_tmp[nslits-1,j,i])
                this_var = numpy.append(
                    numpy.append(
                    var_data_cube_tmp[0,j,i],
                    var_data_cube_tmp[:nslits,j,i]),
                    var_data_cube_tmp[nslits-1,j,i])
                this_dq = numpy.append(
                    numpy.append(
                    dq_data_cube_tmp[0,j,i],
                    dq_data_cube_tmp[:nslits,j,i]),
                    dq_data_cube_tmp[nslits-1,j,i])
                # do interpolation
                f = scipy.interpolate.interp1d(
                    in_x-adr_x, this_flux, kind='linear', 
                    fill_value=0.0, bounds_error=False)
                g = scipy.interpolate.interp1d(
                    in_x-adr_x, this_var, kind='linear', 
                    fill_value=0.0, bounds_error=False)
                h = scipy.interpolate.interp1d(
                    in_x-adr_x, this_dq, kind='nearest', 
                    fill_value=3, bounds_error=False)
                flux_data_cube_tmp[:nslits,j,i] = f(out_x)
                var_data_cube_tmp[:nslits,j,i] = g(out_x)
                dq_data_cube_tmp[:nslits,j,i] = h(out_x)
            if verbose:
                if i > 0 :
                    sys.stdout.flush()
                sys.stdout.write('\r\r %d' % (i/(nlam-1.)*100.) + '%')
                if i == nlam-1 :
                    sys.stdout.write('\n')
    # All done, at last ! Now, let's save it all ...
    # ALWAYS ITERATE TO 25 EVEN IF HALFFRAME HERE!!
    outfits = pyfits.HDUList(f3)
    for i in range(1,26):
        # save to data cube
        outfits[i].data = flux_data_cube_tmp[i-1,:,:]
        outfits[i].header.set('CRVAL1', final_frame_wmin)
        outfits[i].header.set('CDELT1', disp_ave)
        outfits[i+25].data = var_data_cube_tmp[i-1,:,:]
        outfits[i+50].data = dq_data_cube_tmp[i-1,:,:]
        #outfits[i].header.set('NAXIS1', len(out_lambda))
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f3.close()
    return    

#------------------------------------------------------------------------
# scripts to derive a wavelength calibration from standard stars
from wifes_calib import derive_wifes_calibration
from wifes_calib import extract_wifes_stdstar
from wifes_calib import calibrate_wifes_cube
from wifes_calib import derive_wifes_telluric
from wifes_calib import apply_wifes_telluric
from wifes_calib import load_wifes_cube

#------------------------------------------------------------------------
def generate_wifes_3dcube(inimg, outimg):
    # load in data
    # assumes it is in pywifes format
    # otherwise why are you using this function
    f = pyfits.open(inimg)
    if len(f) == 76:
            ny,nlam = numpy.shape(f[1].data)
            nx=25
            # get wavelength array
            lam0 = f[1].header['CRVAL1']
            dlam = f[1].header['CDELT1']
            lam_array = lam0+dlam*numpy.arange(nlam,dtype='d')
            # get data, variance and data quality
            obj_cube_data = numpy.zeros([nlam,ny,nx],dtype='d')
            obj_cube_var  = numpy.zeros([nlam,ny,nx],dtype='d')
            obj_cube_dq   = numpy.zeros([nlam,ny,nx],dtype='d')
            for i in range(nx):
                curr_data = f[i+1].data
                curr_var = f[i+26].data
                curr_dq =f[i+51].data
                obj_cube_data[:,:,nx-i-1] = curr_data.T
                obj_cube_var[:,:,nx-i-1] = curr_var.T
                obj_cube_dq[:,:,nx-i-1] = curr_dq.T
    else:
            nlam, ny, nx = numpy.shape(f[1].data)
            # get wavelength array
            lam0 = f[1].header['CRVAL3']
            dlam = f[1].header['CDELT3']
            lam_array = lam0+dlam*numpy.arange(nlam,dtype='d')
            # get data and variance
            obj_cube_data = numpy.zeros([nlam,ny,nx],dtype='d')
            obj_cube_var  = numpy.zeros([nlam,ny,nx],dtype='d')
            obj_cube_dq   = numpy.zeros([nlam,ny,nx],dtype='d')
            for i in range(nx):
                curr_data = f[1].data[:,:,i]
                curr_var = f[2].data[:,:,i]
                curr_dq =f[3].data[:,:,i]
                obj_cube_data[:,:,nx-i-1] = curr_data
                obj_cube_var[:,:,nx-i-1] = curr_var
                obj_cube_dq[:,:,nx-i-1] = curr_dq
    # save to data cube
    # DATA
    cube_hdu = pyfits.PrimaryHDU(obj_cube_data, header = f[0].header)
    cube_hdu.header.set('CRVAL3',lam0, 'Reference wavelength')
    cube_hdu.header.set('CDELT3',dlam, 'Wavelength step')
    #cube_hdu.header.set('PYWIFES',__version__, 'Pywifes version'))
    outfits = pyfits.HDUList([cube_hdu])
    # VARIANCE
    var_hdu = pyfits.PrimaryHDU(obj_cube_var, header = f[0].header)
    var_hdu.header.set('CRVAL3',lam0, 'Reference wavelength')
    var_hdu.header.set('CDELT3',dlam, 'Wavelength step')
    #var_hdu.header.set('PYWIFES',__version__, 'Pywifes version')
    outfits.append(var_hdu)
    # DQ
    dq_hdu = pyfits.PrimaryHDU(obj_cube_dq, header = f[0].header)
    dq_hdu.header.set('CRVAL3',lam0, 'Reference wavelength')
    dq_hdu.header.set('CDELT3',dlam, 'Wavelength step')
    #dq_hdu.header.set('PYWIFES',__version__, 'Pywifes version')
    outfits.append(dq_hdu)
    # SAVE IT
    outfits[0].header.set('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f.close()
    return

