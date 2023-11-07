from __future__ import print_function
import re
from astropy.io import fits as pyfits
import numpy
import pylab
import pickle
import scipy.interpolate
import math
import multiprocessing
from itertools import cycle
import scipy.optimize as op
# Fred's upadate (wsol)
import os
import datetime

from . import optical_model as om
from .wifes_metadata import metadata_dir
from .wifes_metadata import __version__
from .mpfit import mpfit

#------------------------------------------------------------------------
f0 = open(os.path.join(metadata_dir,'basic_wifes_metadata.pkl'), 'rb')
try:
    wifes_metadata = pickle.load(f0, fix_imports=True, encoding='latin')
except:
    wifes_metadata = pickle.load(f0) # Python 2.7 can't handle fix_imports or encoding=latin
f0.close()
base_wsols = wifes_metadata['baseline_wsols']
all_ref_lines = wifes_metadata['ref_linelists']

#------------------------------------------------------------------------
# high-level functions to check if an observation is half-frame or N+S
def is_halfframe(inimg, data_hdu=0):
    f = pyfits.open(inimg)
    detsec = f[data_hdu].header['DETSEC']
    f.close()
    ystart = int(float(detsec.split(',')[1].split(':')[0]))
    if ystart == 2049:
        return True
    else:
        return False

#------------------------------------------------------------------------
def robust_median(dist,
                  return_sigma = False,
                  return_indices = False):
    c = 0.6745
    data = numpy.copy(dist)
    ind0 = numpy.arange(len(data))
    ind = numpy.arange(len(data))
    prev_npts = -1
    niter_max=10
    niter = 0
    while(len(ind) != prev_npts and niter<niter_max):
        prev_npts = len(ind)
        m = numpy.median(data)
        MAD = numpy.median(numpy.abs(data - m)/c)
        ind = numpy.where(abs(data-m)<4*MAD)[0]
        data = data[ind]
        ind0 = ind0[ind]
        niter += 1
    final_median = numpy.median(data)
    final_sigma = numpy.median(numpy.abs(data - final_median)/c)
    if return_sigma == True and return_indices == True:
        return final_median, final_sigma, ind0
    elif return_sigma == True:
        return final_median, final_sigma
    elif return_indices == True:
        return final_median, ind0
    else:
        return final_median

#------------------------------------------------------------------------
# function to predict baseline wavelength value for x,y
# in a given slitlet for a given grating
def wavelength_guess(x_array, y_array,
                     chosen_slitlet, grating):
    # **NEW** : use default optical model parameters!
    s_array = int(chosen_slitlet)*numpy.ones(len(x_array), dtype='i')
    params = wifes_metadata['optical_params'][grating]
    return om.evaluate_optical_model(
        x_array, y_array, s_array, grating.lower(),
        bin_x=1, bin_y=1, params=params)

def wavelength_guess_poly(x_array, y_array,
                          chosen_slitlet, grating):
    slit_str = '%d' % chosen_slitlet
    results = numpy.zeros(len(x_array), dtype='d')
    ylist = numpy.unique(y_array)
    ylist.sort()
    for y in ylist:
        y_str = '%d' % y
        yinds = numpy.nonzero(y_array==y)[0]
        x = x_array[yinds]
        # get the polynomial coeffs for this slitet and y
        try:
            curr_poly = base_wsols[grating][slit_str][y_str]
            curr_results = numpy.polyval(curr_poly, x)
            results[yinds]=curr_results
        except:
            results[yinds] = -1.0
    return results

# analogous functions to associate found lines with a reference list
def associate_linelists(all_found_lines, # len NF
                        ref_lines,       # len NR
                        y_array = None,  # len NF
                        dlam_cut = 100.0):
    # should only be operated on a single row of data at a time!
    if y_array == None:
        y_array = numpy.ones(len(all_found_lines))
    ylist = numpy.unique(y_array)
    nr = len(ref_lines)
    found_results = -1.0*numpy.ones(len(all_found_lines), dtype='d')
    for y in ylist:
        yinds = numpy.nonzero(y_array == y)[0]
        found_lines = all_found_lines[yinds]
        nf = len(found_lines)
        # get distances between pairs
        td_base = numpy.ones([nr,nf],dtype='d')
        d2 = (((found_lines*td_base).T-ref_lines)**2).T
        # find closest reference line for each found line
        best_ref_inds = d2.argmin(axis=1)
        ref_mask   = numpy.zeros([nr,nf],dtype='d')
        ref_match_inds = (numpy.arange(nr), best_ref_inds)
        ref_mask[ref_match_inds] = 1
        # find closest foundline for each reference line
        best_found_inds = d2.argmin(axis=0)
        found_mask = numpy.zeros([nr,nf],dtype='d')
        found_match_inds = (best_found_inds, numpy.arange(nf))
        found_mask[found_match_inds] = 1
        # find where those two agree!!!
        good_matches = numpy.nonzero(found_mask*ref_mask*(d2<dlam_cut**2))
        good_ref_inds = good_matches[0]
        good_found_inds = good_matches[1]
        # set the output list to a real found wavelength!
        found_results[yinds[good_found_inds]] = ref_lines[good_ref_inds]
    return found_results

#------------------------------------------------------------------------
# arc line fitting functions

#Fred's update ... mpfit and multiprocessing
def gauss_line(p,x):
    return p[0]*numpy.exp(-(x-p[1])**2/(2*p[2]**2))

def err_gauss_line(p,x,y,fjac=None):
    status = 0
    return [status,y - gauss_line(p,x)]

def gauss_line_resid(p,x,y, gain=None, rnoise=10.0):
    if gain==None:
        return (gauss_line(p,x) - y)/rnoise**2
    else:
        return np.sqrt(numpy.maximum(y,0) + rnoise**2)

def lsq_gauss_line( args ):
    """
    Fit a Gaussian to data y(x)
    
    Parameters
    ----------
    args: tuple
        nline, guess_center, width_guess, xfit, yfit
    
    Notes
    -----
    nline: int
        index of this line
    guess_center: float
        initial guess position
    """
    fit = op.least_squares(gauss_line_resid, [numpy.max(args[4]), args[1], args[2]], method='lm', \
            xtol=1e-04, ftol=1e-4, f_scale=[3.,1.,1.], args=(args[3], args[4]))
    #Look for unphysical solutions
    if (fit.x[2]<0.5) or (fit.x[2]>10) or (fit.x[1]<args[3][0]) or (fit.x[1]>args[3][-1]): 
        return args[0], float('nan'), 0., fit.x[2], 0.
    else:
        cov = numpy.linalg.inv(fit.jac.T.dot(fit.jac))
        return args[0], fit.x[1], 1/cov[1,1], fit.x[2], 1/cov[2,2]

def scipy_gauss_line(nline, guess_center, width_guess, xfit, yfit):
    """
    Parameters
    ----------
    nline: int
        index of this line
    guess_center: float
        initial guess position
    """
    return None

def mpfit_gauss_line( packaged_args ):
    (nline, guess_center, width_guess, xfit, yfit) = packaged_args
    # choose x,y subregions for this line
    fa = {'x':xfit, 'y':yfit}
    parinfo = [{'value':yfit[xfit==guess_center], 'fixed':0, 
                'limited':[1,0], 'limits':[0.,0.]},
               {'value':guess_center, 'fixed':0, 
                'limited':[1,1], 
                'limits':[guess_center-width_guess,
                          guess_center+width_guess]},
               {'value':width_guess/2., 'fixed':0, 
                'limited':[1,1], 'limits':[width_guess/20.,width_guess]},
               ]

    my_fit = mpfit(err_gauss_line,functkw=fa, parinfo=parinfo, 
                             quiet=True)
    p1 = my_fit.params
    #print p1, my_fit.status

    # plot for check purposes
    '''
    print my_fit.params, my_fit.status
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(xfit,yfit,'ko-')
    plt.plot(xfit, gauss_line(my_fit.params,xfit),'ro-')
    plt.axvline(my_fit.params[1],c='g')
    plt.axvline(guess_center,c='b',linestyle='--')
    plt.show()
    '''
    if p1[2] == parinfo[2]['limits'][1]:
        # Hum ... line too wide = problem
        return[nline,float('nan')]
    else:
        return[nline,p1[1]]
    # return desired values, for now only care about centers

def weighted_loggauss_arc_fit(subbed_arc_data,
                              peak_centers,
                              width_guess,
                              find_method = 'mpfit',
                              multithread = True): 
    N = len(subbed_arc_data)
    x = numpy.arange(N,dtype='d')
    y = subbed_arc_data
    narc = len(peak_centers)
    fitted_centers = []
    # So, Fred's update (via mpfit) is very much slower ... :-(
    # To make things better, let's do some multiprocessing ...
    if find_method == 'loggauss': #Mike C.'s original code using log(gaussian)
        for i in range(narc):
            curr_ctr_guess = peak_centers[i]
            # choose x,y subregions for this line
            ifit_lo = int(curr_ctr_guess-5*width_guess)
            ifit_hi = int(curr_ctr_guess+5*width_guess)
            xfit = x[ifit_lo:ifit_hi]
            yfit = y[ifit_lo:ifit_hi]
            try:
                good_pix = numpy.nonzero(yfit > 0.2*yfit.max())[0]
            except:
                fitted_centers.append(float('nan'))
                continue
            np = len(good_pix)
            P = numpy.ones([np,3],dtype='d')
            P[:,0] = xfit[good_pix]**2
            P[:,1] = xfit[good_pix]
            weights = yfit[good_pix]
            A_mat = numpy.transpose(numpy.transpose(P)*(weights**2))
            B_mat = numpy.log(yfit[good_pix])*(weights**2)
            results = numpy.linalg.lstsq(A_mat, B_mat)
            # fit parabola to log of flux
            A, B, C = results[0]
            #print 2.0*numpy.sqrt(2.0*numpy.log(2.0)*(-2.0*A)**-1)
            new_xctr = (-1.0*B)/(2.0*A)
            if numpy.abs(new_xctr-curr_ctr_guess) > 2.0*width_guess:
                fitted_centers.append(float('nan'))
            else:
                fitted_centers.append(new_xctr)
        # return desired values, for now only care about centers
    elif find_method =='mpfit' or find_method =='least_squares':
        # Fred's update : fitting a log(gaussian) fails miserably if other 
        # lines are close-by (not even blended !). 
        # Do an mpfit with gaussian function instead (can limit width!)
        # Use multicore to speed things up
        # Mike I : change in a minimal way to use scipy least squares.
        jobs = []
        fitted_centers = numpy.zeros_like(peak_centers, dtype = float)
        for i in range(narc):
            curr_ctr_guess = peak_centers[i]
            ifit_lo = int(curr_ctr_guess-5*width_guess)
            ifit_hi = int(curr_ctr_guess+5*width_guess)
            xfit = x[ifit_lo:ifit_hi]
            yfit = y[ifit_lo:ifit_hi]
            try:
                good_pix = numpy.nonzero(yfit > 0.2*yfit.max())[0]
            except:
                fitted_centers[i] = float('nan')
                continue
            jobs.append( (i,curr_ctr_guess, width_guess, xfit, yfit) )
        if len(jobs)>0:
            if multithread:
                num_available_cpus = len(os.sched_getaffinity(0))
                with multiprocessing.Pool(num_available_cpus) as mypool:
                    chunksize = math.ceil(len(jobs) / num_available_cpus)
                    if find_method =='mpfit':
                        results = mypool.imap_unordered(mpfit_gauss_line, jobs, chunksize=chunksize)
                    else:
                        results = mypool.imap_unordered(lsq_gauss_line, jobs, chunksize=chunksize)

                    # Process the results
                    for r in results:
                        try:
                            this_line = r[0]
                            this_center = r[1]
                        except:
                            continue
                        fitted_centers[this_line] = float(this_center)
            else:
                for i in range(narc):
                    if find_method =='mpfit':
                        # If sentence added by MZ
                        #~ print('FITTED CENTERS', len(fitted_centers))
                        if i<len(jobs):
                            fitted_centers[i] = mpfit_gauss_line(jobs[i])[1]
                        else:
                            pass
                    else:
                        fitted_centers[i] = lsq_gauss_line(jobs[i])[1]

    if find_method == 'loggauss':
        return numpy.array(fitted_centers)
    elif find_method == 'mpfit' or find_method =='least_squares':
        return fitted_centers

def quick_arcline_fit(arc_data,
                      find_method='loggauss',
                      flux_threshold=100.0,
                      deriv_threshold=40.0,
                      width_guess=2.0,
                      flux_saturation = 50000.0,
                      prev_centers = None, # Fred's update (wsol)
                      multithread=False
                      ):
    N = len(arc_data)
    arc_deriv = arc_data[1:]-arc_data[:-1]
    p1_data  = numpy.zeros(N,dtype='d')
    p1_data[:-1] = arc_data[1:]
    m1_data = numpy.zeros(N,dtype='d')
    m1_data[1:] = arc_data[:-1]
    p1_deriv  = numpy.zeros(N,dtype='d')
    p1_deriv[:-1] = arc_deriv
    m1_deriv = numpy.zeros(N,dtype='d')
    m1_deriv[1:] = arc_deriv
    p2_deriv  = numpy.zeros(N,dtype='d')
    p2_deriv[:-2] = arc_deriv[1:]
    m2_deriv = numpy.zeros(N,dtype='d')
    m2_deriv[2:] = arc_deriv[:-1]
    init_inds = numpy.nonzero(
        (arc_data <= flux_saturation)*
        (arc_data >= flux_threshold)*
        (m1_data >= flux_threshold)*
        (p1_data >= flux_threshold)*
        ((m1_deriv >= deriv_threshold)+
         (m2_deriv >= deriv_threshold))*
        ((p1_deriv <= -1.0*deriv_threshold)+
         (p2_deriv <= -1.0*deriv_threshold)))[0]
    # excise adjacent indices
    ind_diffs = init_inds[1:]-init_inds[:-1]
    full_diffs = numpy.zeros(len(init_inds))
    full_diffs[:-1] = ind_diffs
    potential_line_inds = init_inds[
        numpy.nonzero(ind_diffs != 1)[0]]
    
    # Fred's update (wsol)
    # Ok, so roughly the same lines will be selected in each row.
    # 30-50% are rubbsih ... so, if we do it once, we could avoid to do it
    # again and again ...
    # Use Xcorrelation to do that ...
    if (find_method == 'mpfit' or find_method =='least_squares')  and not prev_centers is None:
        # Find the shift between this set of lines and the previous ones
        prev_lines = numpy.zeros(N)
        curr_lines = numpy.zeros(N)
        for j in prev_centers:
            prev_lines[int(j)] = 1.
        for j in potential_line_inds:
            curr_lines[int(j)] = 1.
            
        corr_out = numpy.correlate(prev_lines,curr_lines, 
                                   mode='full')
        shift = (numpy.where(corr_out == numpy.max(corr_out))[0] - N)[0]+1
        # Given the best shift, just select the lines that we could fit with mpfit previously !
        checked_ones = []
        for j in potential_line_inds:
            if (1 <= j+shift < len(prev_lines)-1) and \
               ((prev_lines[j+shift]+curr_lines[j] == 2.0) or \
                    (prev_lines[j+shift+1]+curr_lines[j] == 2.0) or \
                    (prev_lines[j+shift-1]+curr_lines[j] == 2.0)):
                 # these two lines should be in to get 'all' the good lines
                 # in practice, it doesn't matter .. we have enough lines !
                 # This saves ~200 s or so ...
                checked_ones = numpy.append(checked_ones,[j])

        potential_line_inds = numpy.array(checked_ones)
    
    # fit the centers, make sure it didn't fit noise
    next_ctrs = weighted_loggauss_arc_fit(arc_data,
                                          potential_line_inds,
                                          width_guess,
                                          find_method=find_method,
                                          multithread=multithread)
    next_peak_list = []
    for x in next_ctrs:
        if x == x:
            next_peak_list.append(arc_data[int(x)])
        else:
            next_peak_list.append(0.0)
    next_peaks = numpy.array(next_peak_list)
    good_inds = numpy.nonzero(next_peaks >= flux_threshold)[0]
    return next_ctrs[good_inds]

#------------------------------------------------------------------------
# functions for getting wavelength solution for a given set of data
def fit_wsol_poly(x_array, y_array, ref_array, x_polydeg, y_polydeg):
    ncoeffs = x_polydeg+y_polydeg+1
    np = len(x_array)
    base_array = numpy.ones([np, ncoeffs], dtype='d')
    for i in range(x_polydeg):
        base_array[:,x_polydeg-i-1] = x_array**(i+1)
    for j in range(y_polydeg):
        base_array[:,x_polydeg+y_polydeg-j] = y_array**(j+1)
    best_coeffs = numpy.linalg.lstsq(base_array,
                                     numpy.transpose(ref_array))[0]
    x_poly = best_coeffs[:x_polydeg+1]
    y_poly = numpy.zeros(y_polydeg+1,dtype='d')
    y_poly[:y_polydeg] = best_coeffs[x_polydeg+1:]
    return x_poly, y_poly

def evaluate_wsol_poly(x_array, y_array, xpoly, ypoly):
    wave_array = (numpy.polyval(xpoly, x_array)+
                  numpy.polyval(ypoly, y_array))
    return wave_array

def _set_refarray_from_fits(ref_array, init_y_array, fits):
    for this_fit in fits:
        try:
            row_ind = this_fit[0]
            ref_row_inds = (init_y_array == row_ind + 1)
            ref_array[ref_row_inds] = this_fit[2]
        except:
            continue

#------------------------------------------------------------------------
# FUNCTION TO FIND LINES AND GUESS THEIR WAVELENGTHS
# NOTE: OPERATES ON DATA FROM A SINGLE SLITLET!!!
def find_lines_and_guess_refs(slitlet_data,
                              chosen_slitlet,
                              grating,
                              arc_name,
                              find_method='mpfit',
                              shift_method='xcorr_all',
                              ref_arclines=None,
                              dlam_cut_start=5.0,
                              bin_x=1, bin_y=1,
                              verbose=False,
                              yzp=0,
                              flux_threshold_nsig=3.0,
                              deriv_threshold_nsig=1.0,
                              multithread = True,
                              plot=False):
    #-----------------------------------
    # get arclines
    if ref_arclines == None:
        ref_arclines = all_ref_lines[arc_name]
    #-----------------------------------
    # 1 - get some dimensions
    slit_dims = numpy.shape(slitlet_data)
    nrows = slit_dims[0]
    ncols = slit_dims[1]
    linx = numpy.arange(ncols)
    liny = numpy.arange(nrows)
    gridx, gridy = numpy.meshgrid(linx, liny)
    #-----------------------------------
    # 1b - based on pixels statistics get the thresholds
    #      for flux and derivative
    bg_mean, bg_sig = robust_median(slitlet_data.flatten(),
                                    return_sigma=True)
    flux_threshold = bg_mean + flux_threshold_nsig*bg_sig
    # same for x-derivative thresholding
    deriv_data = slitlet_data[:,1:]-slitlet_data[:,:-1]
    df_mean, df_sig = robust_median(deriv_data.flatten(),
                                    return_sigma=True)
    deriv_threshold = df_mean + deriv_threshold_nsig*df_sig
    #-----------------------------------
    # 2 - fit arclines in all rows
    #     calculate expected wavelengths from known wave solution
    #     associate those wavelengths to known lines
    full_fitted_x = []
    full_fitted_y = []
    full_fitted_lam = []
    full_ref_lam = []
    if verbose:
        print(' Slitlet', chosen_slitlet)
        print('  ... detecting arc lines with',find_method,'...')
        start = datetime.datetime.now()

    # Fred's update (wsol)
    # Avoid running mpfit where it is unuseful.
    # Run it once for all detected line in the middle of the slice ...
    # then for all other slices, just re-fit the lines that are real.
    # Do this, and you reduce the total time by 50% for this step !
    if find_method == 'mpfit':
        mid_slit = int(nrows/2./bin_y)
        #~ print('SLITLET_DATA', mid_slit, slitlet_data[mid_slit,:]) # MZ
        mid_fit_centers = quick_arcline_fit(slitlet_data[mid_slit,:],
                                            find_method = find_method,
                                            flux_threshold = flux_threshold,
                                            deriv_threshold = deriv_threshold,
                                            width_guess=2.0,
                                            multithread=multithread)
    else :
        mid_fit_centers = None # don't do it for the loggauss method ...

    #!!! MJI If we cared, this is the part here that should be parallelised.
    for i in range(8//bin_y, nrows-8//bin_y):
        test_z = slitlet_data[i,:]
        fitted_ctrs = quick_arcline_fit(
            test_z,
            find_method=find_method,
            flux_threshold=flux_threshold,
            deriv_threshold=deriv_threshold,
            width_guess=2.0,prev_centers = mid_fit_centers)
        for ctr in fitted_ctrs:
            full_fitted_x.append(ctr)
            full_fitted_y.append(i+1)
    init_y_array = numpy.array(full_fitted_y)
    init_x_array = numpy.array(full_fitted_x)
    if verbose:
        print('  done in',datetime.datetime.now()-start)
    #------------------------------
    # 2b - use xcorr to get shift of the wavelength guess!
    ref_interp = scipy.interpolate.interp1d(
        wifes_metadata['ref_arc_spectra'][grating][arc_name]['x'],
        wifes_metadata['ref_arc_spectra'][grating][arc_name]['y'],
        bounds_error=False,
        fill_value=0.0)
    #---------
    # no xcorr shift
    if shift_method == None:
        temp_wave_array = wavelength_guess(
        init_x_array*bin_x,
        init_y_array*bin_y + yzp*bin_y,
        chosen_slitlet, grating)
    #---------
    # single xcorr version
    elif shift_method == 'xcorr_single':
        mid_flux = slitlet_data[nrows/2,:]
        wave_guess = wavelength_guess(
            numpy.arange(ncols)*bin_x,
            (nrows/2)*bin_y*numpy.ones(ncols) + yzp*bin_y,
            chosen_slitlet, grating)
        wshift = xcorr_shift_single(mid_flux, wave_guess,
                                  ref_interp)
        if verbose:
            print('Slitlet %02d xcorr shift = %.03f' % (chosen_slitlet, wshift))
        temp_wave_array = wavelength_guess(
            init_x_array*bin_x,
            init_y_array*bin_y + yzp*bin_y,
            chosen_slitlet, grating) + wshift
    #---------
    # grid xcorr version
    elif shift_method == 'xcorr_grid':
        wave_guess = numpy.reshape(wavelength_guess(
            gridx.flatten()*bin_x, gridy.flatten()*bin_y,
            chosen_slitlet, grating),
                                   numpy.shape(gridx))
        shift_coeffs, best_stretch = xcorr_shift_grid(
            slitlet_data,
            wave_guess,
            ref_interp)
        w0 = numpy.median(wave_guess)
        init_wguess = wavelength_guess(
            init_x_array*bin_x, init_y_array*bin_y,
            chosen_slitlet, grating)
        temp_wave_array = (
            best_stretch*(init_wguess-w0) + w0 +
            shift_coeffs[1] +
            shift_coeffs[0]*(init_y_array-1))
        if verbose:
            mean_wshift = numpy.mean(temp_wave_array - init_wguess)
            print('Slitlet %02d mean xcorr shift = %.03f' % (
                chosen_slitlet,
                mean_wshift))
    #---------
    # full xcorr version
    elif shift_method == 'xcorr_all':  
        # Fred's update (wsol)
        # First, load the initial guess - currently, only for:
        # B3000, R3000, R7000 with NeAr lamp
        # B3000, R3000, B7000, I7000 with CuAr lamp
        if arc_name != 'NeAr' and arc_name != 'CuAr' :
            print(' Arc lamp not supported for Xcorr identification method !')
            print(" I will crash now... bye !")
            print(' ')
        ref_fn = os.path.join(metadata_dir,
                              'arclines.'+grating+'.'+arc_name+'.txt')
        if os.path.exists(ref_fn) :
            ref_arc = numpy.loadtxt(ref_fn,skiprows=1)
        else : 
            print(' Ref. file for the current arc lamp + grating unavailable.')

        # Create the storage array for the reference wavelength
        ref_array = numpy.zeros_like(init_x_array)
        # Careful here ... some reference files have line position inverted
        # Check which one
        f = open(ref_fn,'r')
        if 'Channel' in f.readline() and 'B' in grating:
            ref_arc[:,0]= ncols - ref_arc[:,0]
        f.close() 

        # Loop through each row (each spectrum of the slice)
        # It is slow, so, I use the multiprocessing package to speed things up
        # Package is installed by default for Python v.2.6 and above.
        # I use the 'Pool' concept - check automatically the number of CPU
        # available. No need to set mutlithread at all ! 
        # Some quick tests showed that the 'Pool' approach is much faster than
        # simply pilling up many 'Process' ... 
        #
        # Fred, 04.2013

        # Stretch value is not varying much over 1 slice.
        # So, get it in the middle, and use it throughout. 
        # Gain some ~6.3 sec per slice by doing this !
        mid_row = int(nrows/2)
        mid_row_ind = (init_y_array == mid_row)
        best_stretch = xcorr_shift_all( (chosen_slitlet, 
                                       mid_row,ncols,mid_row_ind,
                                       init_x_array[mid_row_ind],
                                       ref_arc,None, True,
                                       verbose) )
        jobs = []
        for i in range(nrows):
            row_inds = (init_y_array == i+1)
            # only use rows where enough lines are found
            if len(init_x_array[row_inds]) >= 0.5*len(ref_arc[:,0]):    
                jobs.append( (chosen_slitlet,i,ncols,row_inds,
                              init_x_array[row_inds],ref_arc, 
                              [best_stretch], False, verbose) )
        if multithread:
            num_available_cpus = len(os.sched_getaffinity(0))
            chunksize = math.ceil(len(jobs) / num_available_cpus)
            if verbose:
                print(f'  ... assign lambdas using up to {num_available_cpus} cpu(s) ...')
            with multiprocessing.Pool(num_available_cpus) as pool:
                results = pool.imap_unordered(xcorr_shift_all, jobs, chunksize=chunksize)
                _set_refarray_from_fits(ref_array, init_y_array, results)
        else :
            #Completely avoid using Pool in case of an ipython debugging environment.
            results = [xcorr_shift_all(job) for job in jobs] 
            _set_refarray_from_fits(ref_array, init_y_array, results)

        # Create array with only detected lines
        good_inds = ref_array > 0
        iter_x_array = init_x_array[good_inds]
        iter_y_array = init_y_array[good_inds]
        iter_ref_array = ref_array[good_inds]
        
        if plot :
            plot_detected_lines(chosen_slitlet, iter_x_array, 
                                iter_y_array, iter_ref_array, ncols)
        
        return  iter_x_array,iter_y_array, iter_ref_array
            
    #---------
    else:
        raise ValueError('Xcorr shift method not recognized.')
    #---------
    # 3 - get starting wavelengths!
    init_winds = numpy.nonzero(temp_wave_array > 0.0)[0]
    start_wave_array = temp_wave_array[init_winds]
    start_full_x_array = init_x_array[init_winds]
    start_full_y_array = init_y_array[init_winds]
    niden = len(start_full_x_array)
    #------------------------------
    # 4 - figure out which lines these are, excise outliers!
    start_full_ref_array = associate_linelists(
        start_wave_array,
        ref_arclines,
        y_array = start_full_y_array,
        dlam_cut = dlam_cut_start)
    iter0_good_inds = numpy.nonzero(start_full_ref_array > 0)[0]
    iter_ref_array = start_full_ref_array[iter0_good_inds]
    iter_x_array   = start_full_x_array[iter0_good_inds]
    iter_y_array   = start_full_y_array[iter0_good_inds]

    if plot :
        plot_detected_lines(chosen_slitlet, iter_x_array, iter_y_array, 
                            iter_ref_array,ncols)

    return iter_x_array, iter_y_array, iter_ref_array

def plot_detected_lines(slitlet, x, y, ref,ncols):
      pylab.figure()
      pylab.plot(x, ref, 'k.')
      pylab.xlim([0,ncols])
      pylab.xlabel('Detected arc lines position (spectral dir.) [pixel]')
      pylab.ylabel('Associated wavelength [Angstroem]')
      pylab.title('Slice '+str(slitlet))

      pylab.figure()
      linestyles = ["k.-","r.-","g.-","b.-","c.-","m.-"]
      linestylecycler = cycle(linestyles)
      for a in numpy.unique(ref):
          args = ref == a
          pylab.plot(x[args], y[args], next(linestylecycler))
      pylab.xlim([0,ncols])
      pylab.xlabel('Detected arc lines position (spectral dir.) [pixel]')
      pylab.ylabel('Detected arc lines position (spatial dir.) [pixel]')
      pylab.title('Slice '+str(slitlet))
      pylab.show()

#------------------------------------------------------------------------
# XCORR SHIFT FUNCTIONS
def xcorr_shift_single(flux_data, wave_guess,
                       ref_interp):
    #------------------------------------
    # evaluate the reference spectra at that wavelength array
    ref_spec = ref_interp(wave_guess)
    #------------------------------------
    # do xcorr (later: and stretch?)
    corr = numpy.correlate(flux_data, ref_spec, mode='same')
    corr_lim = 20
    corr[2048+corr_lim+1:] = 0
    corr[:2048-corr_lim] = 0
    init_shift = float(corr.argmax())-2048
    #------------------
    ##print init_shift
    ##pylab.figure()
    ##pylab.plot(flux_data, color='b')
    ##pylab.plot(ref_spec, color='r')
    ##new_mask = numpy.concatenate([flux_data[init_shift:],
    ##                              flux_data[:init_shift]])
    ##pylab.plot(new_mask, color='g')
    ##pylab.show()
    #------------------
    if init_shift > corr_lim:
        shift = init_shift - 4096
    else:
        shift = init_shift
    dw = numpy.mean(wave_guess[1:]-wave_guess[:-1])
    return -1.0*dw*shift

def xcorr_shift_grid(slitlet_data,
                     wave_guess,
                     ref_interp):
    #------------------------------------
    # decide what y-values to fit by xcorr
    ny, nx = numpy.shape(slitlet_data)
    bin_y = int(86.0/float(ny))
    y_samp = numpy.array([14, 24, 34, 44, 54, 64, 74])//bin_y
    n_samp = len(y_samp)
    w0 = numpy.median(wave_guess)
    #------------------------------------
    # for each y-value, find stretch and shift
    #stretches = [0.9980, 0.9985, 0.9990, 0.9995, 1.0000,
    #             1.0005, 1.0010, 1.0015, 1.0020]
    stretches = [1.0000]
    ns = len(stretches)
    all_shifts = numpy.zeros(n_samp, dtype='f')
    all_stretches = numpy.zeros(n_samp, dtype='f')
    for q in range(n_samp):
        y = y_samp[q]
        flux = slitlet_data[y,:]
        wave = wave_guess[y,:]
        orig_dw = numpy.mean(wave[1:]-wave[:-1])
        shifts  = numpy.zeros(ns, dtype='f')
        metrics = numpy.zeros(ns, dtype='f')
        # cover all stretches
        for p in range(ns):
            # reference flux
            ref_flux = ref_interp(stretches[p]*(wave-w0)+w0)
            # correlate it
            corr = numpy.correlate(flux, ref_flux, mode='same')
            corr_lim = 20
            corr[2048+corr_lim+1:] = 0
            corr[:2048-corr_lim] = 0
            init_shift = float(corr.argmax())-2048
            if init_shift > corr_lim:
                shift = init_shift - 4096
            else:
                shift = init_shift
            shifts[p] = stretches[p]*orig_dw
            # metric - just correlation peak?
            metrics[p] = corr.max()
            # metric - resids rms
            #shifted_flux = numpy.concatenate([flux[init_shift:],
            #                                  flux[:init_shift]])
            #resids = ref_flux - shifted_flux
            #metrics[p] = 1.0/numpy.std(resids)
        # find best stretch, save that and shift
        best_s_ind = metrics.argmax()
        all_shifts[q] = shifts[best_s_ind]
        all_stretches[q] = stretches[best_s_ind]
    #------------------------------------
    # based on fitted stretches and shifts,
    # fit dw = Ax+By+C
    #fval_y = numpy.arange(10/bin_y, 80/bin_y)
    #shift_poly = numpy.polyfit(y_samp, all_shifts, 1)
    #shift_fvals = numpy.polyval(shift_poly, fval_y)
    #pylab.figure()
    #pylab.plot(y_samp, all_shifts, 'bo')
    #pylab.plot(fval_y, shift_fvals, 'g')
    #stretch_poly = numpy.polyfit(y_samp, all_stretches, 1)
    #stretch_fvals = numpy.polyval(stretch_poly, fval_y)
    #pylab.figure()
    #pylab.plot(y_samp, all_stretches, 'bo')
    #pylab.plot(fval_y, stretch_fvals, 'g')
    #pylab.show()
    #print shift_poly, stretch_poly
    fval_y = numpy.arange(10//bin_y, 80//bin_y)
    best_stretch = numpy.median(all_stretches)
    good_inds = numpy.nonzero(all_stretches == best_stretch)[0]
    shift_poly = numpy.polyfit(y_samp[good_inds], all_shifts[good_inds], 1)
    #shift_fvals = numpy.polyval(shift_poly, fval_y)
    #pylab.figure()
    #pylab.plot(y_samp, all_shifts, 'ro')
    #pylab.plot(y_samp[good_inds], all_shifts[good_inds], 'bo')
    #pylab.plot(fval_y, shift_fvals, 'g')
    #pylab.show()
    #print shift_poly, best_stretch
    #------------------------------------
    return shift_poly, best_stretch

# Fred's update (wsol)

def xcorr_shift_all( packaged_args ):
    (slitlet_data, this_row, ncols, row_inds, this_init_x_array, this_ref_arc, \
    stretches, # if provided, will only try these ones !
    get_stretch, # if yes, get the best stretch+exit
    verbose ) = packaged_args

    if stretches == None:
        stretches = numpy.arange(0.98,1.03,0.001)

    corrs = numpy.zeros_like(stretches)
    shifts = numpy.zeros_like(stretches)
    # Try different 'stretch and shift values to match the 
    # reference line position. This can be approximative, but not too much !
    for (k,stretch) in enumerate(stretches):
        # Create stretched, temporary "spectrum"
        pseudo_x_obs = numpy.zeros(ncols)
        pseudo_x_bes = numpy.zeros(ncols)
        for j in this_init_x_array:
            pseudo_x_obs[int(j)]=1.
        for j in this_ref_arc[:,0]:
            if int(j*stretch) < len(pseudo_x_bes):
                pseudo_x_bes[int(j*stretch)]=1.

        # Cross-correlate pseudo-ref spectrum with real arc lines
        corr_out = numpy.correlate(pseudo_x_obs,pseudo_x_bes, mode='full')
        corrs[k] = numpy.max(corr_out)
        shifts[k] = numpy.argmax(corr_out)       
        #pylab.plot(pseudo_x_obs,'ko-', markerfacecolor='w')
        #pylab.plot(pseudo_x_bes,'r-')
        #pylab.show()

        #pylab.plot(stretches,corrs, 'k.-')
        #pylab.show()
      
    # Find best shift and stretch based on X correlation results
    best_shift = int(shifts[numpy.argmax(corrs)])
    best_stretch = stretches[numpy.argmax(corrs)]
    
    # Just get the best stretch and return, or keep going.
    if get_stretch:
        return best_stretch
    
    #if verbose :
    #    print '   row:',this_row+1,\
    #        'best_shift/stretch:',best_shift,best_stretch

    # Some more temporary storage structures
    x_obs = numpy.zeros(ncols)
    pseudo_x_obs = numpy.zeros(ncols)
    pseudo_x_bes = numpy.zeros(ncols)
    final_x_bes = numpy.zeros(ncols)
    pseudo_lam_bes = numpy.zeros(ncols)
    final_lam_bes = numpy.zeros(ncols)
         
    # re-create best matching pseudo reference spectrum
    for j in this_init_x_array:
        pseudo_x_obs[int(j)]=1.
        x_obs[int(j)] = j
    for (line,j) in enumerate(this_ref_arc[:,0]):
        # Stretch it ...
        if int(j*best_stretch) < len(pseudo_x_bes):
            pseudo_x_bes[int(j*best_stretch)]=1.
            pseudo_lam_bes[int(j*best_stretch)] = this_ref_arc[:,1][line]
        # Shift it ...
        if best_shift <= ncols :
            final_x_bes[0:best_shift] = pseudo_x_bes[-best_shift:]
            final_lam_bes[0:best_shift] = pseudo_lam_bes[-best_shift:]
        else :
            final_x_bes[best_shift-ncols:] = \
                pseudo_x_bes[0:-(best_shift-ncols)]
            final_lam_bes[best_shift-ncols:] = \
                pseudo_lam_bes[0:-(best_shift-ncols)]
     
    # Make sure there are no other line within 10 Angstroem 
    # on either side ... should maybe be a parameter ...    
    this_ref_array = numpy.zeros_like(this_init_x_array)
    for (j,item) in enumerate(this_init_x_array):
        if j > 0 :
            cond1 = int(item) > ( int(this_init_x_array[j-1]) \
                                      + 10)
        else :
            cond1 = True
        if j < len(this_init_x_array)-1:
            cond2 = int(item) < ( int(this_init_x_array[j+1]) \
                                      - 10)
        else :
            cond2 = True
        # Finally, associate each identified line 
        # with appropriate wavelength.
        if cond1 and cond2 :
            loc = int(numpy.where(x_obs == item)[0])
            if len(final_lam_bes[loc-2:loc+3]
                   [final_lam_bes[loc-2:loc+3]>0])==1:
                this_ref_array[numpy.where(this_init_x_array==item)] = \
                    numpy.max(final_lam_bes[loc-2:loc+3])

    return [this_row,this_init_x_array,this_ref_array]

#------------------------------------------------------------------------
def slitlet_wsol(slitlet_data,
                 chosen_slitlet,
                 grating,
                 arc_name,
                 ref_arclines=None,
                 n_iter=2,
                 dlam_cut_start=7.0,
                 dlam_cut=3.0,
                 bin_x=1, bin_y=1,
                 return_poly=False,
                 x_polydeg=4, y_polydeg=2,
                 flux_threshold_nsig=3.0,
                 deriv_threshold_nsig=1.0,
                 verbose=False):
    #-----------------------------------
    # get arclines
    if ref_arclines == None:
        ref_arclines = all_ref_lines[arc_name]
    #-----------------------------------
    # 1 - get some dimensions
    slit_dims = numpy.shape(slitlet_data)
    nrows = slit_dims[0]
    ncols = slit_dims[1]
    linx = numpy.arange(ncols)
    liny = numpy.arange(nrows)
    gridx, gridy = numpy.meshgrid(linx, liny)
    #-----------------------------------
    # 1b - based on pixels statistics get the thresholds
    #      for flux and derivative
    bg_mean, bg_sig = robust_median(slitlet_data.flatten(),
                                    return_sigma=True)
    flux_threshold = bg_mean + flux_threshold_nsig*bg_sig
    # same for x-derivative thresholding
    deriv_data = slitlet_data[:,1:]-slitlet_data[:,:-1]
    df_mean, df_sig = robust_median(deriv_data.flatten(),
                                    return_sigma=True)
    deriv_threshold = df_mean + deriv_threshold_nsig*df_sig
    #-----------------------------------
    # 2 - fit arclines in all rows
    #     calculate expected wavelengths from known wave solution
    #     associate those wavelengths to known lines
    full_fitted_x = []
    full_fitted_y = []
    full_fitted_lam = []
    full_ref_lam = []
    for i in range(nrows):
        test_z = slitlet_data[i,:]
        fitted_ctrs = quick_arcline_fit(
            test_z,
            flux_threshold=flux_threshold,
            deriv_threshold=deriv_threshold,
            width_guess=2.0)
        for ctr in fitted_ctrs:
            full_fitted_x.append(ctr)
            full_fitted_y.append(i+1)
    init_y_array = numpy.array(full_fitted_y)
    init_x_array = numpy.array(full_fitted_x)
    #------------------------------
    # 2b - use xcorr to get shift of the wavelength guess!
    ref_interp = scipy.interpolate.interp1d(
        wifes_metadata['ref_arc_spectra'][grating][arc_name]['x'],
        wifes_metadata['ref_arc_spectra'][grating][arc_name]['y'],
        bounds_error=False,
        fill_value=0.0)
    #---------
    # single xcorr version
    #mid_flux = slitlet_data[nrows/2,:]
    #wave_guess = wavelength_guess(
    #    numpy.arange(ncols)*bin_x, (nrows/2)*bin_y*numpy.ones(ncols),
    #    chosen_slitlet, grating)
    #wshift = xcorr_shift_single(mid_flux, wave_guess,
    #                          ref_interp)
    #temp_wave_array = wavelength_guess(
    #    init_x_array*bin_x, init_y_array*bin_y,
    #    chosen_slitlet, grating) + wshift
    #---------
    # global xcorr version
    wave_guess = numpy.reshape(wavelength_guess(
        gridx.flatten()*bin_x,
        gridy.flatten()*bin_y,
        chosen_slitlet, grating),
                               numpy.shape(gridx))
    shift_coeffs, best_stretch = xcorr_shift_grid(
        slitlet_data,
        wave_guess,
        ref_interp)
    init_wguess = wavelength_guess(
        init_x_array*bin_x, init_y_array*bin_y,
        chosen_slitlet, grating)
    temp_wave_array = (
        best_stretch*init_wguess + 
        shift_coeffs[1] +
        shift_coeffs[0]*(init_y_array-1))
    #print numpy.mean(temp_wave_array - init_wguess)
    #print squirrel
    #------------------------------
    # 3 - get starting wavelengths!
    init_winds = numpy.nonzero(temp_wave_array > 0.0)[0]
    start_wave_array = temp_wave_array[init_winds]
    start_full_x_array = init_x_array[init_winds]
    start_full_y_array = init_y_array[init_winds]
    niden = len(start_full_x_array)
    #------------------------------
    # 4 - figure out which lines these are, excise outliers!
    start_full_ref_array = associate_linelists(
        start_wave_array,
        ref_arclines,
        y_array = start_full_y_array,
        dlam_cut = dlam_cut_start)
    iter0_good_inds = numpy.nonzero(start_full_ref_array > 0)[0]
    # START OLD CODE
    #nref = len(ref_arclines)
    #td_base = numpy.ones([nref,niden],dtype='d')
    #d2 = ((start_wave_array*td_base).T-ref_arclines)**2
    #iden_inds = d2.argmin(axis=1)
    #start_full_ref_array = ref_arclines[iden_inds]
    #dlam = (start_wave_array - start_full_ref_array)
    #iter0_good_inds = numpy.nonzero(numpy.abs(dlam) < dlam_cut_start)[0]
    # END OLD CODE
    iter_ref_array = start_full_ref_array[iter0_good_inds]
    iter_x_array   = start_full_x_array[iter0_good_inds]
    iter_y_array   = start_full_y_array[iter0_good_inds]
    #-----------------------------------------------
    # 5 - iterate!
    for i in range(n_iter):
        iter_xpoly, iter_ypoly = fit_wsol_poly(
            iter_x_array,
            iter_y_array,
            iter_ref_array,
            x_polydeg=x_polydeg, y_polydeg=y_polydeg)
        iter_fit_lam = (numpy.polyval(iter_xpoly, iter_x_array)+
                         numpy.polyval(iter_ypoly, iter_y_array))
        iter_resids = iter_fit_lam - iter_ref_array
        # try to bring back all line with the new wavelength solution
        next_wave_array = (numpy.polyval(iter_xpoly, start_full_x_array)+
                           numpy.polyval(iter_ypoly, start_full_y_array))
        next_dlam = (next_wave_array - start_full_ref_array)
        iter_good_inds = numpy.nonzero(numpy.abs(next_dlam) < dlam_cut)[0]
        iter_x_array   = start_full_x_array[iter_good_inds]
        iter_y_array   = start_full_y_array[iter_good_inds]
        iter_ref_array = start_full_ref_array[iter_good_inds]
    #------------------
    if verbose:
        resids = next_dlam[iter_good_inds]
        resid_rms = numpy.mean(resids**2)**0.5
        resid_MAD = numpy.median(numpy.abs(resids)/0.6745)
        full_x = numpy.arange(ncols)
        full_y = numpy.arange(nrows)
        mesh_x, mesh_y = numpy.meshgrid(full_x, full_y)
        iter_fvals = (numpy.polyval(iter_xpoly, mesh_x)+
                      numpy.polyval(iter_ypoly, mesh_y))
        wmin = numpy.min(iter_fvals)
        wmax = numpy.max(iter_fvals)
        #print chosen_slitlet, resid_rms, len(iter_good_inds), wmin, wmax
        print(chosen_slitlet, resid_MAD, len(iter_good_inds), wmin, wmax)
    #------------------
    # FINAL RESULTS!
    full_x_array = iter_x_array
    full_y_array = iter_y_array
    full_ref_array = iter_ref_array
    full_xpoly = iter_xpoly
    full_ypoly = iter_ypoly
    if return_poly:
        return full_y_array, full_x_array, full_ref_array, full_xpoly, full_ypoly
    else:
        return full_y_array, full_x_array, full_ref_array

#------------------------------------------------------------------------
def derive_wifes_polynomial_wave_solution(inimg,
                                          out_file,
                                          dlam_cut_start=7.0,
                                          dlam_cut=3.0,
                                          ref_arclines=None,
                                          ref_arcline_file=None,
                                          arc_name=None,
                                          grating=None,
                                          bin_x=None, bin_y=None,
                                          x_polydeg=4, y_polydeg=2,
                                          flux_threshold_nsig=3.0,
                                          deriv_threshold_nsig=1.0,
                                          verbose=False):
    # check if halfframe
    halfframe = is_halfframe(inimg)
    # MUST HAVE MEF FILE AS INPUT
    a = pyfits.open(inimg)
    arc_hdr = a[0].header
    outfits = pyfits.HDUList([pyfits.PrimaryHDU(header=arc_hdr)])
    ndy, ndx = numpy.shape(a[1].data)
    xarr = numpy.arange(ndx)
    yarr = numpy.arange(ndy)
    full_x, full_y = numpy.meshgrid(xarr, yarr)
    # check binning
    try:
        bin_temp = arc_hdr['CCDSUM'].split()
        default_bin_x = int(float(bin_temp[0]))
        default_bin_y = int(float(bin_temp[1]))
    except:
        default_bin_x = 1
        default_bin_y = 1
    if bin_x == None:
        bin_x = default_bin_x
    if bin_y == None:
        bin_y = default_bin_y
    # get grating, and arc name
    if grating == None:
        if arc_hdr['CAMERA'] == 'WiFeSRed':
            grating = arc_hdr['GRATINGR']
        else:
            grating = arc_hdr['GRATINGB']
    if arc_name == None:
        init_arc_name = arc_hdr['LAMP']
        next_arc_name = re.sub('-', '', init_arc_name)
        again_arc_name = re.sub(' ', '', next_arc_name)
        arc_name = re.sub('_', '', again_arc_name)
    # set the arc linelist!
    if ref_arcline_file != None:
        f1 = open(ref_arcline_file, 'r')
        ref_arclines = numpy.array([float(line.split()[0])
                                    for line in f1.readlines()])
        f1.close()
    # fit each slitlet
    for i in range(1,26):
        if halfframe and i>12:
            wave_data = numpy.zeros(numpy.shape(a[i].data), dtype='d')
        else:
            new_y, new_x, new_ref, xpoly, ypoly = slitlet_wsol(
                a[i].data,
                i,
                grating,
                arc_name,
                ref_arclines=ref_arclines,
                dlam_cut_start=dlam_cut_start,
                dlam_cut=dlam_cut,
                return_poly=True,
                bin_x=bin_x,
                bin_y=bin_y,
                x_polydeg=x_polydeg,
                y_polydeg=y_polydeg,
                flux_threshold_nsig=flux_threshold_nsig,
                deriv_threshold_nsig=deriv_threshold_nsig,
                verbose=verbose)
            wave_data = (numpy.polyval(xpoly, full_x)+
                         numpy.polyval(ypoly, full_y))
        # calculate lambda for all x/y!!
        hdu_name = 'WSOL%d' % i
        # maybe save poly coeffs as header keywords??
        new_hdu = pyfits.ImageHDU(wave_data, arc_hdr, name=hdu_name)
        outfits.append(new_hdu)
    outfits[0].header.update('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(out_file, overwrite=True)
    a.close()
    return

def save_found_lines(inimg,
                     out_file,
                     dlam_cut_start=7.0,
                     dlam_cut=3.0,
                     ref_arclines=None,
                     ref_arcline_file=None,
                     arc_name=None,
                     grating=None,
                     bin_x=None, bin_y=None,
                     x_polydeg=4, y_polydeg=2,
                     flux_threshold_nsig=3.0,
                     deriv_threshold_nsig=1.0,
                     verbose=False):
    # check if halfframe
    halfframe = is_halfframe(inimg)
    # MUST HAVE MEF FILE AS INPUT
    a = pyfits.open(inimg)
    arc_hdr = a[0].header
    ndy, ndx = numpy.shape(a[1].data)
    xarr = numpy.arange(ndx)
    yarr = numpy.arange(ndy)
    full_x, full_y = numpy.meshgrid(xarr, yarr)
    # check binning
    try:
        bin_temp = arc_hdr['CCDSUM'].split()
        default_bin_x = int(float(bin_temp[0]))
        default_bin_y = int(float(bin_temp[1]))
    except:
        default_bin_x = 1
        default_bin_y = 1
    if bin_x == None:
        bin_x = default_bin_x
    if bin_y == None:
        bin_y = default_bin_y
    # get grating, and arc name
    if grating == None:
        if arc_hdr['CAMERA'] == 'WiFeSRed':
            grating = arc_hdr['GRATINGR']
        else:
            grating = arc_hdr['GRATINGB']
    if arc_name == None:
        init_arc_name = arc_hdr['LAMP']
        next_arc_name = re.sub('-', '', init_arc_name)
        again_arc_name = re.sub(' ', '', next_arc_name)
        arc_name = re.sub('_', '', again_arc_name)
    # set the arc linelist!
    if ref_arcline_file != None:
        f1 = open(ref_arcline_file, 'r')
        ref_arclines = numpy.array([float(line.split()[0])
                                    for line in f1.readlines()])
        f1.close()
    # fit each slitlet
    fitted_lines = []
    for i in range(1,26):
        if halfframe and i>12:
            continue
        else:
            new_y, new_x, new_ref, xpoly, ypoly = slitlet_wsol(
                a[i].data,
                i,
                grating,
                arc_name,
                ref_arclines=ref_arclines,
                dlam_cut_start=dlam_cut_start,
                dlam_cut=dlam_cut,
                return_poly=True,
                bin_x=bin_x,
                bin_y=bin_y,
                x_polydeg=x_polydeg,
                y_polydeg=y_polydeg,
                flux_threshold_nsig=flux_threshold_nsig,
                deriv_threshold_nsig=deriv_threshold_nsig,
                verbose=verbose)
            new_dict = {
                'slitlet_number' : i*numpy.ones(len(new_x)),
                'found_x' : new_x,
                'found_y' : new_y,
                'found_ref' : new_ref,
                'xpoly' : xpoly,
                'ypoly' : ypoly}
            fitted_lines.append(new_dict)
    f = open(out_file, 'wb')
    pickle.dump(fitted_lines, f)
    f.close()
    return

#------------------------------------------------------------------------
# wrapper for fitting sky lines!
def derive_wifes_skyline_solution(inimg, out_file,
                                  dlam_cut=5.0,
                                  grating=None,
                                  bin_x=None, bin_y=None,
                                  x_polydeg=3,
                                  y_polydeg=1,
                                  verbose=False):
    # set the line list!
    sky_lines = wifes_metadata['ref_linelists']['sky']
    # run the normal script
    return derive_wifes_wave_solution(inimg, out_file,
                                      dlam_cut=dlam_cut,
                                      ref_arclines=sky_lines,
                                      grating=grating,
                                      bin_x=bin_y, bin_y=bin_y,
                                      x_polydeg=x_polydeg,
                                      y_polydeg=y_polydeg,
                                      verbose=verbose)

#------------------------------------------------------------------------
# OPTICAL MODEL FUNCTION(S)

# Tolerance for stopping the fit.  0.001 rmse should be good enough I think.
FTOL = 1e-3

def excludeLines(lines, exclude, index=3, epsilon=0.05):
  """ Work out if any of the lines we read in are to be excluded
      We do this by constructing a matrix of differences between the input lines
      and the exclude lines, then determining whether any are close enough to
      count as a match.  We then exclude those. """
  print('MZ: EXCLUDE LINES', lines, exclude)
  if ((not exclude is None) and len(exclude) > 0):
    print('Excluding', exclude, 'with tolerance of', epsilon)
    exclude = numpy.asarray(exclude)
    nlines = lines.shape[0]
    nexclude = exclude.shape[0]
    keeplines = (numpy.abs(numpy.repeat(lines[:,index], nexclude).reshape(nlines,nexclude) - numpy.tile(exclude, (nlines,1))) > epsilon).all(axis=1)
    return lines[keeplines]
  else:
    return lines

def _fit_optical_model(title, grating, bin_x, bin_y, lines, alphap, doalphapfit, doplot, automatic, sigma, verbose, decimate, save_prefix=None):

  # Don't do the alphap fit initially
  doalphapfit_thistime = False

  # Set up the initial set of parameters
  plorig = om.defaultParams(grating)

  # Extract individual arrays from the main array
  alls, ally, allx, allarcs = om.extractArrays(lines, grating, bin_x, bin_y)

  lambda0 = plorig[13]
  print('lambda0=',lambda0)

  # Get sensible values for xdc and ydc
  # by looking at the central slitlet
  args = (alls == 13)
  tmparcs = allarcs[args]
  if (len(tmparcs) > 0):
    largs = (tmparcs <= lambda0)
    rargs = (tmparcs >= lambda0)
    if (numpy.any(tmparcs[largs]) and numpy.any(tmparcs[rargs])):
      left = numpy.argmax(tmparcs[largs])
      right = numpy.argmin(tmparcs[rargs])
      plorig[11] = (allx[args][largs][left]+allx[args][rargs][right])/2
      print('xdc=',plorig[11])

    plorig[12] = (ally[args].max() + ally[args].min()) / 2.
    print('ydc=',plorig[12])

  # Initial residuals
  resid = om.errfunc(grating, plorig, alphap, alls, ally, allx, allarcs)

  # Do an initial plot of the residuals if required
  if doplot:
    om.plotFunc(title,allx,ally,allarcs,om.fitfunc(grating,plorig,alphap,alls,ally,allx))
    om.plotResid(title,allx,ally,allarcs,resid)

  # The initial RMSE before fitting
  var = numpy.sum(resid**2) / len(allx)
  bias = numpy.sum(resid) / len(allx)
  rmse = math.sqrt(var + bias**2)
  print('MZ: diagnostics', var, bias, rmse)
  print("Initial RMSE",rmse)

  # Set up parameter info ready for fitting
  parinfo = [{'value':0., 'fixed':1, 'limited':[0,0], 'limits':[0.,0.]} for i in range(om.nparams + len(alphap))]

  # Set the parameter values
  for i,v in enumerate(numpy.concatenate((plorig,alphap))):
    parinfo[i]['value']=v

  # Limit the input alpha angle
  parinfo[1]['limited']=[1,1]
  parinfo[1]['limits']=[0, math.pi/4]

  # Limit phi
  parinfo[2]['limited']=[1,1]
  parinfo[2]['limits']=[-math.pi/2,math.pi/2]

  # The sets of parameters to fit in each iteration
  if (grating[1:] == '3000'):
    paramlist = [(0,13),(1,),(2,3,4,5,6,7,8,9,10,11,12,14,15,16,17)]
  else:
    paramlist = [(1,),(2,3,4,5,6,7,8,9,10,11,12,16,17)]

  # Work with decimated data first, if asked to
  if (decimate):
    if (verbose):
      print('Working with decimated data')
    origLines = lines
    lines = lines[::10]

  finished = False
  while (not finished):

    # Check if this is to be the last run through
    if automatic <= 0:
      # Not automatic or interactive, so we're done
      finished = True
      if (doalphapfit and not doalphapfit_thistime):
        # Unless we want an alphapfit
        doalphapfit_thistime = True
        paramlist.append(tuple(range(om.nparams,om.nparams+len(alphap)),))

    # Report on how many points are in the data set
    if (verbose):
      print("Working with",len(lines),"data points")

    # Check that we actually have a sensible number of lines
    if (len(lines) < 100):
      return (None, None, None, None, None)

    # Extract the columns we need
    alls, ally, allx, allarcs = om.extractArrays(lines, grating, bin_x, bin_y)

    # FIXME what is the real error?
    err = numpy.ones_like(ally)
    fa = {'s':alls, 'y':ally, 'x':allx, 'grating':grating, 'arc':allarcs, 'err':err}

    for params in paramlist:
      for r in range(1):
        if (verbose):
          print('Fitting for parameters',params)
        # Fix all parameters
        for p in parinfo:
          p['fixed'] = 1
        # Unfix only those we want to fit this time around
        for i in params:
          parinfo[i]['fixed'] = 0
        # Do the fit
        fitdone = False
        fitcount = 0
        # Maximum number of times to repeat the fit
        MAXFITS = 1
        while (not fitdone):
          fitcount += 1
          # Actually do the fit
          m = mpfit(om.mpfitfunc, functkw=fa, parinfo=parinfo, iterfunct=None, ftol=FTOL)
          # Report on it
          if (verbose):
            print('status = ', m.status)
          if (m.status <= 0):
             print('error message = ', m.errmsg)

          if (verbose):
            # Work out the RMSE
            chisq = m.fnorm
            dof=len(allx)-len(m.params)
            rmse=numpy.sqrt(chisq/dof)
            print("RMSE",rmse)

          # Copy back fitted parameters into parinfo structure
          for i,v in enumerate(m.params):
            parinfo[i]['value'] = v
            # Report the ones we just changed
            if (parinfo[i]['fixed'] == 0) and verbose:
              print(i, parinfo[i]['value'])

          # Repeat the fit if we need more steps
          if ((m.status == 5) and (fitcount < MAXFITS)):
            print('mpfit needs more steps; repeating fit')
          else:
            fitdone = True

        if doplot:
          pl = numpy.asarray(m.params)[:om.nparams]
          resid = om.errfunc(grating, pl, alphap, alls, ally, allx, allarcs)
          om.plotResid(title,allx,ally,allarcs,resid)
          #om.plotFunc(title,allx,ally,allarcs,om.fitfunc(grating, pl,alphap,alls,ally,allx))

    # Extract the fitted parameters
    pl = numpy.asarray(m.params)
    if (verbose):
      print('Fit complete')
      om.printParams(grating, pl[:om.nparams], pl[om.nparams:])

    # Generate the final residuals
    resid = om.errfunc(grating, pl[:om.nparams],pl[om.nparams:],alls,ally,allx,allarcs)

    # Calculate how good the fit was
    var = numpy.sum(resid**2) / len(allx)
    bias = numpy.sum(resid) / len(allx)
    rmse = math.sqrt(var + bias**2)

    if (verbose):
      print("VAR",var)
      print("BIAS",bias)
      print("RMSE",rmse)

    # If we were working with decimated data, we now go back to using
    # the full data set for the next run.
    if (decimate):
      if (verbose):
        print('Working with full data set')
      lines = origLines
      decimate = False
      alls, ally, allx, allarcs = om.extractArrays(lines, grating, bin_x, bin_y)
      resid = om.errfunc(grating, pl[:om.nparams],pl[om.nparams:],alls,ally,allx,allarcs)
      if (doplot):
        om.plotResid(title,allx,ally,allarcs,resid)

    # In automatic mode we select some of the lines to be removed
    if automatic > 0:
      lines = om.excludeAuto(lines, grating, bin_x, bin_y, resid, sigma, doplot, verbose)
      alls, ally, allx, allarcs = om.extractArrays(lines, grating, bin_x, bin_y)
      err = numpy.ones_like(ally)
      fa = {'s':alls, 'y':ally, 'x':allx, 'arc':allarcs, 'err':err}
      automatic -= 1

    elif doplot:
      # Do a plot of the residuals if required
      if save_prefix != None:
          save_fn = save_prefix+'final_resids.png'
      else:
          save_fn = None
      om.plotResid(title,allx,ally,allarcs,resid, save_fn=save_fn)

  # Do a final plot of the data if required
  if doplot:
    alls, ally, allx, allarcs = om.extractArrays(lines, grating, bin_x, bin_y)
    if save_prefix != None:
        save_fn = save_prefix+'final_wsol.png'
    else:
        save_fn = None
    om.plotLines(title,allx,ally, save_fn=save_fn)

  print("Final RMSE",rmse)
  return (allx, ally, alls, allarcs, pl)

def derive_wifes_optical_wave_solution(inimg,
                                       outfn,
                                       # line finding parameters
                                       arc_name=None,
                                       ref_arclines=None,
                                       ref_arcline_file=None,
                                       dlam_cut_start=5.0,
                                       flux_threshold_nsig=3.0,
                                       find_method='loggauss',
                                       shift_method='xcorr_single',
                                       # optical model parameters
                                       exclude_from=None,
                                       exclude=None,
                                       epsilon=0.005,
                                       doalphapfit=False,
                                       automatic=False,
                                       verbose=False,
                                       decimate=False,
                                       sigma=1.0,
                                       alphapfile=None,
                                       #global parameters
                                       doplot=False,
                                       savefigs=False,
                                       save_prefix='wsol_',
                                       multithread=True):
  """ The main user-callable function that performs the fit"""
  #------------------------------------------------------
  # *** Mike's edits: operate on PyWiFeS MEF files ***
  #------------------------------------------------------
  # step 1 - gather metadata from header
  f = pyfits.open(inimg, ignore_missing_end=True) # MZ: added ignore_missing_end=True, but this is supposed to work only for python 2 but not 3.
  camera = f[1].header['CAMERA']
  if camera == 'WiFeSRed':
      grating = f[1].header['GRATINGR']
  else:
      grating = f[1].header['GRATINGB']
  ccdsum = f[1].header['CCDSUM']
  bin_x = int(float(ccdsum.split()[0]))
  bin_y = int(float(ccdsum.split()[1]))

  # Get some optional meta-data
  dateobs = f[1].header.get('DATE-OBS')
  tdk = f[1].header.get('TDK')
  pmb = f[1].header.get('PMB')
  rh = f[1].header.get('RH')
  rma = f[1].header.get('ROTSKYPA') # Dumb name for rotator mechanical angle
  if arc_name == None:
      init_arc_name = f[1].header['LAMP']
      next_arc_name = re.sub('-', '', init_arc_name)
      again_arc_name = re.sub(' ', '', next_arc_name)
      arc_name = re.sub('_', '', again_arc_name)
  # set the arc linelist!
  if ref_arcline_file != None:
      f1 = open(ref_arcline_file, 'r')
      ref_arclines = numpy.array([float(line.split()[0])
                                  for line in f1.readlines()])
      f1.close()
  # step 2 - find lines!
  found_x_lists = []
  found_y_lists = []
  found_s_lists = []
  found_r_lists = []
  yrange = []
  for i in range(1,26):
      # and the yrange...
      detsec = f[i].header['DETSEC']
      y0 = int(detsec.split(',')[1].split(':')[0])
      y1 = int(detsec.split(',')[1].split(':')[1].split(']')[0])
      # Make the values work nicely with the range() command
      if (y0 > y1):
          ystop = y0 + 1
          ystart = y1
      else:
          ystop = y1 + 1
          ystart = y0
      # Plot them or not ?
      if doplot == True or ((doplot != False) and ('step1' in doplot)):
          step1plot = True
      else:
          step1plot = False
      if doplot == True or (doplot != False) and (('step2' in doplot)):
          step2plot = True
      else:
          step2plot = False    

      # guess the reference wavelengths
      new_x, new_y, new_r = find_lines_and_guess_refs(
          f[i].data,
          i,
          grating,
          arc_name,
          find_method=find_method,
          shift_method=shift_method,
          ref_arclines=ref_arclines,
          dlam_cut_start=dlam_cut_start,
          bin_x=bin_x, bin_y=bin_y,
          yzp=ystart,
          verbose=verbose,
          flux_threshold_nsig=flux_threshold_nsig,
          deriv_threshold_nsig=1.0,
          multithread=multithread,
          plot=step1plot)
      nl = len(new_x)
      found_x_lists.append(new_x)
      found_r_lists.append(new_r)
      found_s_lists.append(i*numpy.ones(nl))
      yrange.append((ystart,ystop))
      found_y_lists.append(new_y+ystart)
  f.close()
  all_x = numpy.concatenate(found_x_lists)
  all_y = numpy.concatenate(found_y_lists)
  all_s = numpy.concatenate(found_s_lists)
  all_r = numpy.concatenate(found_r_lists)
  if verbose:
      print('Line finding complete')
  # NEED TO (FOR NOW) HAVE COMPLIANCE WITH NIELSEN 'LINES' TEMPLATE
  lines = numpy.column_stack((all_s, all_y, all_x, all_r))
  grating = grating.lower()
  alls, ally, allx, allarcs = om.extractArrays(lines, grating, bin_x, bin_y)
  #------------------------------------------------------
  # NOTE: THIS CODE BELOW IS DEPRECATED
  #------------------------------------------------------
  # Read in file containing yrange and all line data
  #try:
  #  with open(idfile, "r") as f:
  #    (grating, (bin_x, bin_y), yrange, lines) = pickle.load(f)
  #except IOError:
  #  print "Warning:", idfile, "could not be opened"
  #  return
  #
  #grating = grating.lower()
  #
  ## Extract the separate arrays from the line data
  #alls, ally, allx, allarcs = om.extractArrays(lines, grating, bin_x, bin_y)
  #
  #------------------------------------------------------
  # Read in lines to exclude
  if exclude_from is None:
    exclude = numpy.empty(0)
  else:
    exclude = numpy.loadtxt(exclude_from)

  # Add explicitly specified lines
  if (not exclude is None):
    exclude = numpy.append(exclude, exclude)

  # Exclude the lines from the data set
  lines = excludeLines(lines, exclude, index=3, epsilon=epsilon)
  print('MZ: LINES', lines)

  # Set parameters
  alphap = None

  # Read in alphap
  if (not alphapfile is None):
    try:
      alphap = numpy.loadtxt(alphapfile)
      print('Using alphap',alphap)
    except IOError:
      pass

  if alphap is None:
    if (verbose):
      print('Using default alphap (all zero)')
    alphap = numpy.zeros(25)

  if (verbose):
    print('Grating',grating)

  title = inimg.split('/')[-1][:-5] + ' - '+grating

  if savefigs:
      final_save_prefix = save_prefix
  else:
      final_save_prefix = None
  allx, ally, alls, allarcs, params = _fit_optical_model(title, grating, bin_x, bin_y, lines, alphap, doalphapfit, step2plot, automatic, sigma, verbose, decimate, final_save_prefix)

  if not (params is None):
    # Dump some output
    newlines = numpy.column_stack((alls,ally,allx,allarcs,om.fitfunc(grating, params[:om.nparams],params[om.nparams:],alls,ally,allx)))
    om.saveData(outfn+"_extra.pkl", grating, params, newlines, (dateobs, tdk, pmb, rh, rma))

    # And the resampling data
    om.saveResamplingData(outfn, yrange, grating, bin_x, bin_y, params)

  return

#------------------------------------------------------------------------
def derive_wifes_wave_solution(inimg,
                               out_file,
                               method='optical',
                               **args):
    if method == 'poly':
        derive_wifes_polynomial_wave_solution(
            inimg, out_file, **args)
    elif method == 'optical':
        derive_wifes_optical_wave_solution(
            inimg, out_file, **args)
    else:
        raise ValueError('Wavelength solution method not recognized')
