import numpy
from scipy.interpolate import splrep
import pickle
from astropy.io import fits as pyfits
import scipy.interpolate
import wifes_ephemeris
from wifes_metadata import metadata_dir
from wifes_metadata import __version__

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# reference star information!
stdstar_fn = metadata_dir+'stdstar_lookup_table.dat'
f1 = open(stdstar_fn, 'r')
stdstar_lines = f1.readlines()[1:]
f1.close()

ref_fname_lookup = {}
ref_coords_lookup = {}
for line in stdstar_lines:
    lns = line.split()
    fn = lns[0]
    name = lns[1]
    radec = '%s %s' % (lns[2], lns[3])
    ref_fname_lookup[name] = fn
    ref_coords_lookup[name] = radec

# extinction interpolation object
extinct_fn = metadata_dir+'sso_extinction.dat'
extinct_data = numpy.loadtxt(extinct_fn)
sso_extinct_interp = scipy.interpolate.interp1d(extinct_data[:,0],
                                                extinct_data[:,1],
                                                bounds_error=False,
                                                fill_value=numpy.nan)

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# high-level function to find nearest standard star for a given frame!
stdstar_list = ref_coords_lookup.keys()
stdstar_list.sort()
nstds = len(stdstar_list)
stdstar_ra_array  = numpy.zeros(nstds, dtype='d')
stdstar_dec_array = numpy.zeros(nstds, dtype='d')
for i in range(nstds):
    stdstar_radec = ref_coords_lookup[stdstar_list[i]]
    stdstar_ra, stdstar_dec = wifes_ephemeris.sex2dd(stdstar_radec)
    stdstar_ra_array[i] = stdstar_ra
    stdstar_dec_array[i] = stdstar_dec

def find_nearest_stdstar(inimg, data_hdu=0):
    f = pyfits.open(inimg)
    radec = '%s %s' % (f[data_hdu].header['RA'],
                       f[data_hdu].header['DEC'])
    f.close()
    ra, dec = wifes_ephemeris.sex2dd(radec)
    angsep_array = 3600.0*(
        (dec-stdstar_dec_array)**2+
        (numpy.cos(numpy.radians(dec))*(ra-stdstar_ra_array))**2)**0.5
    best_ind = numpy.argmin(angsep_array)
    return stdstar_list[best_ind], angsep_array[best_ind]

#------------------------------------------------------------------------
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

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# scripts for masking out certain wavelength regions
def wavelength_mask(wave_array, band_list):
    mask = numpy.ones(len(wave_array))
    for band in band_list:
        mask *= ((wave_array <= band[0])+
                 (wave_array >= band[1]))
    return mask

# O2  bands are saturated - don't depend on airmass
# H2O bands DO depend on airmass!!
O2_telluric_bands = [
    [6856.0, 6956.0],
    [7584.0, 7693.0]]
    #[7547.0, 7693.0]]

H2O_telluric_bands = [
    [6270.0, 6290.0],
    [7154.0, 7332.0],
    [8114.0, 8344.0],
    [8937.0, 9194.0],
    [9270.0, 9776.0]]

strong_H2O_telluric_bands = [
    [6270.0, 6290.0],
    [7154.0, 7332.0],
    [8114.0, 8344.0],
    [8937.0, 9194.0],
    [9270.0, 9400.0]]

master_H2O_telluric_bands = [
    [5870.0, 6000.0],
    [6270.0, 6290.0],
    [6459.0, 6598.0],
    [7154.0, 7332.0],
    [8114.0, 8344.0],
    [8937.0, 9194.0],
    [9270.0, 9776.0]]

# functions for masking telluric and halpha features
def strong_telluric_mask(wave_array):
    return wavelength_mask(wave_array,
                           O2_telluric_bands + strong_H2O_telluric_bands)

def telluric_mask(wave_array):
    return wavelength_mask(wave_array,
                           O2_telluric_bands + H2O_telluric_bands)

def halpha_mask(wave_array):
    ha_bands = [[6550.0, 6575.0]]
    return wavelength_mask(wave_array, ha_bands)

#------------------------------------------------------------------------
def load_wifes_cube(cube_fn, ytrim=[0,0]):
    f = pyfits.open(cube_fn)
    ny,nlam = numpy.shape(f[1].data)
    nx=25
    # get wavelength array
    lam0 = f[1].header['CRVAL1']
    dlam = f[1].header['CDELT1']
    lam_array = lam0+dlam*numpy.arange(nlam,dtype='d')
    # get data and variance
    obj_cube_data = numpy.zeros([nlam,ny-sum(ytrim),nx],dtype='d')
    obj_cube_var  = numpy.zeros([nlam,ny-sum(ytrim),nx],dtype='d')
    for i in range(nx):
        curr_data = f[i+1].data[ytrim[0]:ny-ytrim[1],:]
        curr_var = f[i+26].data[ytrim[0]:ny-ytrim[1],:]
        obj_cube_data[:,:,i] = curr_data.T
        obj_cube_var[:,:,i] = curr_var.T
    f.close()
    # return flux, variance, wavelength
    return obj_cube_data, obj_cube_var, lam_array

#------------------------------------------------------------------------
def extract_wifes_stdstar(cube_fn,
                          x_ctr = None, y_ctr = None,
                          extract_radius = 5.0,
                          sky_radius = 8.0,
                          ytrim=0,
                          save_mode=None,
                          save_fn=None):
    # extraction parameters...
    slice_size_arcsec = 1.0
    # check if halfframe
    halfframe = is_halfframe(cube_fn)
    # check spatial binning!
    f = pyfits.open(cube_fn)
    exptime = float(f[1].header['EXPTIME'])
    bin_y = float(f[1].header['CCDSUM'].split()[1])
    pix_size_arcsec = bin_y*0.5
    # load the cube data
    init_obj_cube_data, init_obj_cube_var, lam_array = load_wifes_cube(cube_fn)
    inlam, iny, inx = numpy.shape(init_obj_cube_data)
    obj_cube_data = init_obj_cube_data[:,ytrim:iny-ytrim,:]
    obj_cube_var  = init_obj_cube_var[:,ytrim:iny-ytrim,:]
    nlam, ny, nx = numpy.shape(obj_cube_data)
    # get stdstar centroid
    lin_x = numpy.arange(nx,dtype='d')
    lin_y = numpy.arange(ny,dtype='d')
    full_x, full_y = numpy.meshgrid(lin_x, lin_y)
    cube_x = full_x*numpy.ones([nlam,ny,nx])
    cube_y = full_y*numpy.ones([nlam,ny,nx])
    flux = numpy.sum(numpy.sum(obj_cube_data,axis=2),axis=1)
    # make a mask for halfframe
    if halfframe:
        halfframe_mask = full_y < 12.5
    else:
        halfframe_mask = full_y < 50
    # centroid version
    #if x_ctr == None:
    #    std_x = numpy.sum(numpy.sum(
    #        obj_cube_data*cube_x,axis=1),axis=1)/flux
    #else:
    #    std_x = x_ctr*numpy.ones(nlam, dtype='d')
    #if y_ctr == None:
    #    std_y = numpy.sum(numpy.sum(
    #        obj_cube_data*cube_y,axis=2),axis=1)/flux
    #else:
    #    std_y = y_ctr*numpy.ones(nlam, dtype='d')
    if x_ctr == None or y_ctr == None:
        cube_im = numpy.sum(obj_cube_data, axis=0)
        maxind = numpy.nonzero(cube_im == cube_im.max()) 
        yc = maxind[0][0]
        xc = maxind[1][0]
        std_x = xc*numpy.ones(nlam, dtype='d')
        std_y = yc*numpy.ones(nlam, dtype='d')
    else:
        std_x = x_ctr*numpy.ones(nlam, dtype='d')
        std_y = y_ctr*numpy.ones(nlam, dtype='d')
    # fit smooth curves!
    polydeg=6
    lin_lam = numpy.arange(nlam,dtype='d')
    xpoly = numpy.polyfit(lin_lam, std_x, polydeg)
    xfvals = numpy.polyval(xpoly, lin_lam)
    ypoly = numpy.polyfit(lin_lam, std_y, polydeg)
    yfvals = numpy.polyval(ypoly, lin_lam)
    # now extract
    std_flux = numpy.zeros(nlam,dtype='d')
    sky_flux = numpy.zeros(nlam,dtype='d')
    std_var  = numpy.zeros(nlam,dtype='d')
    # THERE IS A FASTER WAY TO DO THIS WITH MASKING...
    for i in range(nlam):
        # get *distance* of each pixels from stdstar center x/y
        pix_dists = ((slice_size_arcsec*(full_x-std_x[i]))**2+
                     (pix_size_arcsec*(full_y-std_y[i]))**2)**0.5
        curr_data = obj_cube_data[i,:,:]
        curr_var  = obj_cube_var[i,:,:]
        # find mean sky spectrum outside sky radius
        sky_pix = numpy.nonzero(
            (pix_dists >= sky_radius)*halfframe_mask)
        nspix = len(sky_pix[0])
        #curr_sky_flux = numpy.sum(curr_data[sky_pix])/float(nspix)
        curr_sky_flux = numpy.median(curr_data[sky_pix])
        # subtract sky and sum up obj flux
        obj_pix = numpy.nonzero(
            (pix_dists <= extract_radius)*halfframe_mask)
        obj_flux = numpy.sum(curr_data[obj_pix]-curr_sky_flux)
        obj_var = numpy.sum(curr_var[obj_pix])
        std_flux[i] = obj_flux
        sky_flux[i] = curr_sky_flux
        std_var[i] = obj_var
    # DIVIDE FLUX BY EXPTIME AND BIN SIZE!!!
    dlam = lam_array[1]-lam_array[0]
    fscale = exptime*dlam
    std_flux /= fscale
    std_var  /= (fscale**2)
    sky_flux /= fscale
    # return flux or save!
    if save_mode == None:
        f.close()
        return lam_array, std_flux
    elif save_mode == 'ascii':
        f.close()
        save_data = numpy.zeros([nlam,3],dtype='d')
        save_data[:,0] = lam_array
        save_data[:,1] = std_flux
        save_data[:,2] = std_var
        numpy.savetxt(save_fn, save_data)
    elif save_mode == 'iraf':
        out_header = f[1].header
        out_header.update('CD1_1', f[1].header['CDELT1'], savecomment=True)
        out_header.update('CD2_2', 1)
        out_header.update('CD3_3', 1)
        out_header.update('LTM3_3', 1)
        out_data = numpy.zeros([4,1,nlam],dtype='d')
        out_data[0,0,:] = std_flux
        out_data[1,0,:] = sky_flux
        out_data[2,0,:] = std_var
        out_data[3,0,:] = std_var
        out_hdu = pyfits.PrimaryHDU(
            data=out_data,
            header=out_header)
        outfits = pyfits.HDUList([out_hdu])
        outfits[0].header.update('PYWIFES', __version__, 'PyWiFeS version')
        outfits.writeto(save_fn, overwrite=True)
        f.close()
    else:
        f.close()
        raise ValueError, 'Standard Star save format not recognized'

# here write wrapper function to save the extracted star spectrum

#------------------------------------------------------------------------
# simple function to divide a cube by some spectrum
def wifes_cube_divide(inimg, outimg,
                      corr_wave, corr_flux):
    corr_interp = scipy.interpolate.interp1d(
        corr_wave, corr_flux,
        bounds_error=False,
        fill_value=numpy.inf) # set divided flux outside bounds to zero
    f3 = pyfits.open(inimg)
    # get the wavelength array
    wave0 = f3[1].header['CRVAL1']
    dwave = f3[1].header['CDELT1']
    nlam = numpy.shape(f3[1].data)[1]
    wave_array = wave0+dwave*numpy.arange(nlam,dtype='d')
    # calculate the flux calibration array
    fcal_array = corr_interp(wave_array)
    outfits = pyfits.HDUList(f3)
    for i in range(1,26):
        curr_flux = f3[i].data
        curr_var  = f3[25+i].data
        out_flux = curr_flux / fcal_array
        out_var  = curr_var / (fcal_array**2)
        # save to data cube
        outfits[i].data = out_flux
        outfits[i+25].data = out_var
    outfits[0].header.update('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f3.close()
    return


# -------------------------- Fred's update ------------------------------
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
        """
        Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
        The Savitzky-Golay filter removes high frequency noise from data.
        It has the advantage of preserving the original shape and
        features of the signal better than other types of filtering
        approaches, such as moving averages techniques.
        Parameters
        ----------
        y : array_like, shape (N,)
           the values of the time history of the signal.
       window_size : int
           the length of the window. Must be an odd integer number.
       order : int
           the order of the polynomial used in the filtering.
           Must be less then `window_size` - 1.
       deriv: int
           the order of the derivative to compute (default = 0 means only smoothing)
       Returns
       -------
       ys : ndarray, shape (N)
           the smoothed signal (or it's n-th derivative).
       Notes
       -----
       The Savitzky-Golay is a type of low-pass filter, particularly
       suited for smoothing noisy data. The main idea behind this
       approach is to make for each point a least-square fit with a
       polynomial of high order over a odd-sized window centered at
       the point.
       Examples
       --------
       t = np.linspace(-4, 4, 500)
       y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
       ysg = savitzky_golay(y, window_size=31, order=4)
       import matplotlib.pyplot as plt
       plt.plot(t, y, label='Noisy signal')
       plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
       plt.plot(t, ysg, 'r', label='Filtered signal')
       plt.legend()
       plt.show()
       References
       ----------
       .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
          Data by Simplified Least Squares Procedures. Analytical
          Chemistry, 1964, 36 (8), pp 1627-1639.
       .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
          W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
          Cambridge University Press ISBN-13: 9780521880688
       """
        from math import factorial   
        try:
           window_size = numpy.abs(numpy.int(window_size))
           order = numpy.abs(numpy.int(order))
        except ValueError, msg:
            raise ValueError("window_size and order have to be of type int")
        if window_size % 2 != 1 or window_size < 1:
            raise TypeError("window_size size must be a positive odd number")
        if window_size < order + 2:
            raise TypeError("window_size is too small for the polynomials order")
        order_range = range(order+1)
        half_window = (window_size -1) // 2
        # precompute coefficients
        b = numpy.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = numpy.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = y[0] - numpy.abs( y[1:half_window+1][::-1] - y[0] )
        lastvals = y[-1] + numpy.abs(y[-half_window-1:-1][::-1] - y[-1])
        y = numpy.concatenate((firstvals, y, lastvals))
        return numpy.convolve( m[::-1], y, mode='valid')

#------------------------------------------------------------------------
def derive_wifes_calibration(cube_fn_list,
                             calib_out_fn,
                             stdstar_name_list=None,
                             extract_in_list=None,
                             airmass_list=None,
                             ref_dir=metadata_dir,
                             ref_fname_list=None,
                             plot_stars=False,
                             plot_sensf=False,
                             norm_stars=False,
                             method = 'poly',
                             polydeg=30,
                             excise_cut = 0.5,
                             wave_min=None,
                             wave_max=None,
                             extinction_fn=None,
                             ytrim=5):
    if plot_stars or plot_sensf:
        import pylab
    # get extinction curve
    if extinction_fn == None:
        extinct_interp = sso_extinct_interp
    else:
        ext_data = numpy.loadtxt(extinction_fn)
        extinct_interp = scipy.interpolate.interp1d(ext_data[:,0],
                                                    ext_data[:,1],
                                                    bounds_error=False,
                                                    fill_value=numpy.nan)
    # first extract stdstar spectra and compare to reference
    fratio_results = []
    for i in range(len(cube_fn_list)):
        f = pyfits.open(cube_fn_list[i])
        cube_hdr = f[1].header
        f.close()
        #------------------------------------
        # figure out which star it is
        # NEW VERSION 0.7.0: smart star name lookup!
        #
        # top priority: user forces the name
        if stdstar_name_list != None:
            star_name = stdstar_name_list[i]
            # if you forced an unknown star name, reset name to None
            if star_name not in ref_fname_lookup.keys():
                star_name = None
        else:
            star_name = None
        # try to find the nearest standard in the list
        if star_name == None:
            try:
                star_name, dist = find_nearest_stdstar(cube_fn_list[i])
                if dist > 200.0:
                    star_name = None
            except:
                # last resort: use the object name from the fits header
                # and pray it's correct
                star_name = cube_hdr['OBJECT']
        #------------------------------------
        #print star_name
        if airmass_list != None:
            secz = airmass_list[i]
        else:
            try:
                secz = cube_hdr['AIRMASS']
            except:
                print 'AIRMASS header missing for %s' % cube_fn_list[i].split('/')[-1]
                secz = 1.0
        # check if there is a calib spectrum...
        if ref_fname_list != None:
            ref_fname = ref_name_list[i]
        elif star_name in ref_fname_lookup.keys():
            ref_fname = ref_fname_lookup[star_name]
        else:
            continue
        # get observed data
        if extract_in_list == None:
            obs_wave, obs_flux = extract_wifes_stdstar(cube_fn_list[i],
                                                       ytrim=ytrim)
        else:
            ex_data = numpy.loadtxt(extract_in_list[i])
            obs_wave = ex_data[:,0]
            obs_flux = ex_data[:,1]
        if wave_min == None:
            wave_min = numpy.min(obs_wave)
        if wave_max == None:
            wave_max = numpy.max(obs_wave)
        # get reference data
        ref_data = numpy.loadtxt(ref_dir+ref_fname)
        ref_interp = scipy.interpolate.interp1d(
            ref_data[:,0], ref_data[:,1],
            bounds_error=False, fill_value=numpy.nan)
        ref_flux = ref_interp(obs_wave)
        std_ext = extinct_interp(obs_wave)
        good_inds = numpy.nonzero((ref_flux==ref_flux)*
                                  (std_ext==std_ext)*
                                  (obs_wave >= wave_min)*
                                  (obs_wave <= wave_max)*
                                  (obs_flux > 0.0))[0]
        init_flux_ratio = -2.5*numpy.log10(obs_flux[good_inds] /
                                           ref_flux[good_inds])
        flux_ratio = init_flux_ratio + (secz-1.0)*std_ext[good_inds]
        fratio_results.append([obs_wave[good_inds], init_flux_ratio])
        if plot_stars:
            scaled_flux = obs_flux[good_inds]/numpy.mean(
                10.0**(-0.4*flux_ratio))
            pylab.figure()
            pylab.plot(obs_wave, ref_flux, color='b',
                       label='Reference star flux')
            pylab.plot(obs_wave[good_inds], scaled_flux, color='r', 
                       label='Scaled observed flux')
            pylab.title(star_name)
            pylab.xlabel(r'Wavelength [$\AA$]')
            pylab.legend(loc='lower right', fancybox=True,shadow=True)
    # from all comparisons, derive a calibration solution
    # EVENTUALLY WILL FIT AN EXTINCTION TERM TOO
    if norm_stars:
        i_mid = len(fratio_results[0][0])/2
        fscale_max = min([x[1][i_mid] for x in fratio_results])
        init_full_y = numpy.concatenate(
            [x[1]-x[1][i_mid]+fscale_max for x in fratio_results])
    else:
        init_full_y = numpy.concatenate(
            [x[1] for x in fratio_results])
    init_full_x = numpy.concatenate(
        [x[0] for x in fratio_results])
    init_good_inds = numpy.nonzero(
        (init_full_y == init_full_y)*
        (init_full_y < numpy.median(init_full_y)+20.0)*
        (strong_telluric_mask(init_full_x))*
        (halpha_mask(init_full_x)))[0]
    # do a first fit
    next_full_y = init_full_y[init_good_inds]
    next_full_x = init_full_x[init_good_inds]
    sort_order = next_full_x.argsort()
    temp_full_x = next_full_x[sort_order]
    temp_full_y = next_full_y[sort_order]
    # ----------- Fred's update 3 -------------------  
    if method == 'smooth_SG':
        # Savitzky-Golay requires continuous data. ->need to fill the 'holes' 
        # It is a problem for red spectra (at this point at least)
        # Check if there are gaps (telluric, Halpha, etc ...)
        init_bad_inds = \
            numpy.nonzero(1-
                          ((init_full_y == init_full_y)*
                           (init_full_y < numpy.median(init_full_y)+20.0)*
                           (telluric_mask(init_full_x))*
                           (halpha_mask(init_full_x))))[0]
        if len(init_bad_inds) > 0 : 
            # if yes, first fit a polynomial, then use it to 'fill the gaps.
            temp_calib = numpy.polyfit(temp_full_x, temp_full_y, polydeg)
            temp_fvals = numpy.polyval(temp_calib, init_full_x)
            init_full_y[init_bad_inds] = temp_fvals[init_bad_inds]
            temp_full_y = init_full_y # to ensure this case is then compatible
            temp_full_x = init_full_x
            # Fails if multiple stars ... need to order the array !
            this_sort_order = temp_full_x.argsort()
            temp_full_x = temp_full_x[this_sort_order]
            temp_full_y = temp_full_y[this_sort_order]        
        # Then fit SG normally
        temp_fvals = savitzky_golay(temp_full_y,101,1,0)
        excise_cut = 0.003
    else :
        temp_best_calib = numpy.polyfit(temp_full_x, temp_full_y, polydeg)
        temp_fvals = numpy.polyval(temp_best_calib, temp_full_x) 
    # excise outliers
    final_good_inds = numpy.nonzero(
        numpy.abs(temp_fvals-temp_full_y)/numpy.abs(temp_fvals) < excise_cut)[0]
    full_x = temp_full_x[final_good_inds]
    full_y = temp_full_y[final_good_inds]
    # ------------ Fred's update 3 ----------------
    if  method == 'smooth_SG': # Fails if multiple stars ... need to order the array !
        this_sort_order = full_x.argsort()
        full_x = full_x[this_sort_order]
        full_y = full_y[this_sort_order]
        final_fvals = savitzky_golay(full_y,101,1,0)
        this_f = scipy.interpolate.interp1d(
            full_x,final_fvals,bounds_error=False, 
            kind='linear')
        all_final_fvals = this_f(init_full_x)
        final_x = full_x
        final_y = full_y
    else :
        best_calib = numpy.polyfit(full_x, full_y, polydeg)
        final_fvals = numpy.polyval(best_calib, full_x)
        final_x = numpy.arange(
            numpy.min(full_x),
            1.000001*numpy.max(full_x),
            0.0001*(numpy.max(full_x)-numpy.min(full_x)))
        final_y = numpy.polyval(best_calib, final_x)
    # plot if requested
    if plot_sensf:
        pylab.figure()
        # MC update - raw fit on top
        pylab.axes([0.10, 0.35, 0.85, 0.60])
        pylab.plot(temp_full_x, temp_full_y, 'r.',markerfacecolor='none',
                   markeredgecolor='r', 
                   label='Raw sensitivity (initial regions)')
        pylab.plot(full_x, full_y, color='b', 
                   label ='Raw sensitivity (valid regions)')
        pylab.plot(temp_full_x,temp_fvals, color=r'#FF6103',
                   lw=2,label='Initial fit')
        if  method == 'smooth_SG':
            pylab.plot(init_full_x, all_final_fvals, color=r'#00FF00', lw=2, 
                       label='Final fit') 
        else :
            pylab.plot(full_x, final_fvals, color=r'#00FF00', lw=2, 
                       label='Final fit')
        #pylab.hlines(-37.5,numpy.min(full_x),numpy.max(full_x), 'k')
        pylab.xlim([numpy.min(full_x),numpy.max(full_x)])
        curr_ylim = pylab.ylim()
        curr_xlim = pylab.xlim()
        pylab.ylim(curr_ylim[::-1])
        pylab.title('Derived sensitivity function')
        pylab.legend(loc='lower right',fancybox=True, shadow=True)
        # lower plot - residuals!
        pylab.axes([0.10, 0.10, 0.85, 0.25])
        pylab.plot(full_x,full_y - final_fvals,'k.', mec=r'#666666',
                   markerfacecolor='none', label='Residuals')
        pylab.axhline(0.0, color='k')
        pylab.xlim(curr_xlim)
        pylab.ylim([-0.2, 0.2])
        pylab.xlabel(r'Wavelength [$\AA$]')
        pylab.ylabel('Residuals')
    if plot_stars or plot_sensf:
        pylab.show()
    # Fred's update ... now, careful, because that's dirty ... 
    # the function does not always return the same thing !
    # SAVE IN THE PICKLE FILE THE WAVELENGTH AND CALIB FVAL ARRAYS
    save_calib = {'wave' : final_x,
                  'cal'  : final_y}
    f1 = open(calib_out_fn, 'w')
    pickle.dump(save_calib, f1)
    f1.close()
    return

#------------------------------------------------------------------------
def calibrate_wifes_cube(inimg, outimg,
                         calib_fn,
                         mode='pywifes',
                         extinction_fn=None):
    # get extinction curve
    if extinction_fn == None:
        extinct_interp = sso_extinct_interp
    else:
        ext_data = numpy.loadtxt(extinction_fn)
        extinct_interp = scipy.interpolate.interp1d(ext_data[:,0],
                                                    ext_data[:,1],
                                                    bounds_error=False,
                                                    fill_value=numpy.nan)
    # open data
    f3 = pyfits.open(inimg)
    # get the wavelength array
    wave0 = f3[1].header['CRVAL1']
    dwave = f3[1].header['CDELT1']
    exptime = f3[1].header['EXPTIME']
    try:
        secz = f3[1].header['AIRMASS']
    except:
        secz = 1.0
        print 'AIRMASS keyword not found, assuming airmass=1.0'
    nlam = numpy.shape(f3[1].data)[1]
    wave_array = wave0+dwave*numpy.arange(nlam,dtype='d')
    # calculate the flux calibration array
    if mode == 'pywifes':
        f1 = open(calib_fn, 'r')
        calib_info = pickle.load(f1)
        f1.close()
        sort_order = calib_info['wave'].argsort()
        calib_x = calib_info['wave'][sort_order]
        calib_y = calib_info['cal'][sort_order]
        # Fred's update 2 ... to account for the smooth method for B3000 ...
        # That's to fix the ugly of the previous function ...
        this_f = scipy.interpolate.interp1d(
            calib_x,calib_y,
            bounds_error=False,fill_value=-100.0,
            kind='linear')
        all_final_fvals = this_f(wave_array)
        inst_fcal_array = 10.0**(-0.4*all_final_fvals) 
        #import pdb
        #pdb.set_trace()
    elif mode == 'iraf':
        f = pyfits.open(calib_fn)
        calib_wave = (
            f[0].header['CRVAL1']
            +f[0].header['CDELT1']
            *numpy.arange(f[0].header['NAXIS1'], dtype='d'))
        calib_flux = f[0].data
        calib_interp = scipy.interpolate.interp1d(
            calib_wave, calib_flux,
            bounds_error=False, fill_value=0.0)
        inst_fcal_array = calib_interp(wave_array)
        f.close()
    else:
        raise ValueError, 'Calibration mode not defined'
    # calculate extinction curve for observed airmass
    obj_ext = 10.0**(-0.4*((secz-1.0)*extinct_interp(wave_array)))
    fcal_array = inst_fcal_array*obj_ext
    # apply flux cal to data!
    outfits = pyfits.HDUList(f3)
    for i in range(1,26):
        curr_flux = f3[i].data
        curr_var  = f3[25+i].data
        out_flux = curr_flux / (fcal_array*exptime*dwave)
        # Fred's update 2 : add dwave !
        out_var  = curr_var / ((fcal_array*exptime*dwave)**2)
        # Fred'supadte 2 : add dwave !a
        # save to data cube
        outfits[i].data = out_flux
        outfits[i+25].data = out_var
    outfits[0].header.update('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f3.close()
    return  

#------------------------------------------------------------------------
# telluric corrections!!!
def derive_wifes_telluric(cube_fn_list,
                          out_fn,
                          plot=False,
                          plot_stars=False,
                          extract_in_list=None,
                          airmass_list=None,
                          telluric_threshold=0.97,
                          fit_wmin =  5400.0,
                          fit_wmax = 10000.0,
                          H2O_power=0.72,
                          O2_power=0.40,
                          polydeg=4,
                          ytrim=3):
    if plot_stars or plot:
        import pylab
    #---------------------------------------------
    # for each star, get its airmass if not specified in input
    if airmass_list == None:
        airmass_list = []
        for fn in cube_fn_list:
            try:
                f = pyfits.open(fn)
                new_am = float(f[1].header['AIRMASS'])
                f.close()
            except:
                new_am = 1.0
                print 'AIRMASS keyword not found, assuming 1.0'
            airmass_list.append(new_am)
    #---------------------------------------------
    # now extract each star spectrum and derive telluric correction spectra
    O2_corrections  = []
    H2O_corrections = []
    for i in range(len(cube_fn_list)):
        # get extracted spectrum
        if extract_in_list == None:
            obs_wave, obs_flux = extract_wifes_stdstar(cube_fn_list[i],
                                                       ytrim=ytrim)
        else:
            ex_data = numpy.loadtxt(extract_in_list[i])
            obs_wave = ex_data[:,0]
            obs_flux = ex_data[:,1]
        # define all the telluric regions
        O2_mask  = wavelength_mask(obs_wave, O2_telluric_bands)
        H2O_mask = wavelength_mask(obs_wave, H2O_telluric_bands)
        O2_inds  = numpy.nonzero(O2_mask==0)[0]
        H2O_inds = numpy.nonzero(H2O_mask==0)[0]
        # fit smooth polynomial to non-telluric regions!
        fit_inds = numpy.nonzero(O2_mask*H2O_mask
                                 *(obs_wave >= fit_wmin)
                                 *(obs_wave <= fit_wmax))[0]
        smooth_poly = numpy.polyfit(obs_wave[fit_inds],
                                    obs_flux[fit_inds],
                                    polydeg)
        # get ratio of data to smooth continuum
        smooth_cont = numpy.polyval(smooth_poly, obs_wave)
        init_ratio = obs_flux / smooth_cont
        if plot_stars:
            pylab.figure()
            pylab.plot(obs_wave, obs_flux, 'b')
            pylab.plot(obs_wave, smooth_cont, 'g')
        # isolate desired regions, apply thresholds!
        O2_ratio = numpy.ones(len(obs_wave), dtype='d')
        O2_ratio[O2_inds] = init_ratio[O2_inds]
        O2_ratio[numpy.nonzero(O2_ratio >= telluric_threshold)[0]] = 1.0
        O2_corrections.append([obs_wave, O2_ratio])
        H2O_ratio = numpy.ones(len(obs_wave), dtype='d')
        H2O_ratio[H2O_inds] = init_ratio[H2O_inds]
        H2O_ratio[numpy.nonzero(H2O_ratio >= telluric_threshold)[0]] = 1.0
        H2O_corrections.append([obs_wave, H2O_ratio])
    #---------------------------------------------
    # now using all, derive the appropriate solutions!
    # wavelength range shouldn't change much, use the first one!
    base_wave = O2_corrections[0][0]
    O2_corr_temp = numpy.zeros([len(cube_fn_list),
                                len(base_wave)], dtype='d')
    H2O_corr_temp = numpy.zeros([len(cube_fn_list),
                                 len(base_wave)], dtype='d')
    O2_corr_temp[0,:]  = O2_corrections[0][1]
    H2O_corr_temp[0,:] = (H2O_corrections[0][1])**(
        1.0/(airmass_list[0]**0.55))
    for i in range(1, len(cube_fn_list)):
        O2_interp = scipy.interpolate.interp1d(
            O2_corrections[i][0],
            (O2_corrections[i][1])**(
            1.0/(airmass_list[i]**O2_power)),
            bounds_error=False, fill_value=1.0)
        O2_corr_temp[i,:] = O2_interp(base_wave)
        H2O_interp = scipy.interpolate.interp1d(
            H2O_corrections[i][0],
            (H2O_corrections[i][1])**(
            1.0/(airmass_list[i]**H2O_power)),
            bounds_error=False, fill_value=1.0)
        H2O_corr_temp[i,:] = H2O_interp(base_wave)
    final_O2_corr = numpy.mean(O2_corr_temp, axis=0)
    final_H2O_corr = numpy.mean(H2O_corr_temp, axis=0)
    # fix zero values
    final_O2_corr[numpy.nonzero(final_O2_corr < 0.01)[0]] = 0.01
    final_H2O_corr[numpy.nonzero(final_H2O_corr < 0.01)[0]] = 0.01
    # fix nan values
    final_O2_corr[numpy.nonzero(final_O2_corr != final_O2_corr)[0]] = 1.0
    final_H2O_corr[numpy.nonzero(final_H2O_corr != final_H2O_corr)[0]] = 1.0
    #---------------------------------------------
    # PLOT FOR INSPECTION...
    if plot:
        fig1=pylab.figure()
        sp1 = fig1.add_subplot(111)
        fig2=pylab.figure()
        sp2 = fig2.add_subplot(111)
        sp1.plot(base_wave, final_O2_corr, color='k', lw=2, label='Default O2')
        sp2.plot(base_wave, final_H2O_corr, color='k', lw=2, label='Default H2O')
        for i in range(len(cube_fn_list)):
            airmass = airmass_list[i]
            wave, O2_ratio = O2_corrections[i]
            H2O_ratio = H2O_corrections[i][1]
            sp1.plot(wave, O2_ratio**(1.0/(airmass**O2_power)), 
                     label=cube_fn_list[i].split('/')[-1])
            sp2.plot(wave, H2O_ratio**(1.0/(airmass**H2O_power)), 
                     label=cube_fn_list[i].split('/')[-1])
        sp1.set_xlabel(r'Wavelength [$\AA$]')
        sp1.set_title('Individual O2 Telluric Correction functions')
        sp1.legend(loc='lower left',fancybox=True, shadow=True)
        sp1.set_xlim([numpy.min(base_wave),numpy.max(base_wave)])
        sp2.set_xlabel(r'Wavelength [$\AA$]')
        sp2.set_title('Individual H2O Telluric Correction functions')
        sp2.legend(loc='lower left',fancybox=True, shadow=True)
        sp2.set_xlim([numpy.min(base_wave),numpy.max(base_wave)])
        pylab.figure()
        pylab.plot(base_wave, final_O2_corr, color='r', lw=2, 
                   label='Default O2 lines' )
        pylab.plot(base_wave, final_H2O_corr, color='k', lw=2, 
                   label='Default H2O lines')
        pylab.xlabel(r'Wavelength [$\AA$]')
        pylab.legend(loc='lower left',fancybox=True,shadow=True)
        pylab.xlim([numpy.min(base_wave),numpy.max(base_wave)])
        pylab.title('Telluric Correction functions')
        pylab.show()
    #---------------------------------------------
    # save to output file!
    tellcorr_info = {
        'wave' : base_wave,
        'O2'   : final_O2_corr,
        'H2O'  : final_H2O_corr,
        'O2_power' : O2_power,
        'H2O_power' : H2O_power}
    f1 = open(out_fn, 'w')
    pickle.dump(tellcorr_info, f1)
    f1.close()
    return

def apply_wifes_telluric(inimg,
                         outimg,
                         tellcorr_fn,
                         airmass=None):
    #---------------------------------------------
    # open the telluric corrction file
    f1 = open(tellcorr_fn)
    tellcorr_info = pickle.load(f1)
    O2_interp = scipy.interpolate.interp1d(
        tellcorr_info['wave'],
        tellcorr_info['O2'],
        bounds_error=False, fill_value=1.0)
    H2O_interp = scipy.interpolate.interp1d(
        tellcorr_info['wave'],
        tellcorr_info['H2O'],
        bounds_error=False, fill_value=1.0)
    try:
        O2_power = tellcorr_info['O2_power']
        H2O_power = tellcorr_info['H2O_power']
    except:
        O2_power = 0.55
        H2O_power = 1.0
    f1.close()
    #---------------------------------------------
    # apply to chosen data
    f3 = pyfits.open(inimg)
    # get airmass
    if airmass == None:
        try:
            airmass = float(f3[1].header['AIRMASS'])
        except:
            airmass = 1.0
            print 'AIRMASS keyword not found, assuming airmass=1.0'
    # get the wavelength array
    wave0 = f3[1].header['CRVAL1']
    dwave = f3[1].header['CDELT1']
    nlam = numpy.shape(f3[1].data)[1]
    wave_array = wave0+dwave*numpy.arange(nlam,dtype='d')
    # calculate the telluric correction array
    base_O2_corr = O2_interp(wave_array)
    O2_corr = base_O2_corr**(airmass**O2_power)
    base_H2O_corr = H2O_interp(wave_array)
    H2O_corr = base_H2O_corr**(airmass**H2O_power)
    fcal_array = O2_corr*H2O_corr
    # correct the data
    outfits = pyfits.HDUList(f3)
    for i in range(1,26):
        curr_flux = f3[i].data
        curr_var  = f3[25+i].data
        out_flux = curr_flux / fcal_array
        out_var  = curr_var / (fcal_array**2)
        # save to data cube
        outfits[i].data = out_flux
        outfits[i+25].data = out_var
    outfits[0].header.update('PYWIFES', __version__, 'PyWiFeS version')
    outfits.writeto(outimg, overwrite=True)
    f3.close()
    return


