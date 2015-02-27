import os
import re
import numpy
import pylab
import scipy.interpolate

#------------------------------------------------------------------------
#snid_dir = '/Users/mchildress/mjc_lib/snid-5.0/'
#snid_temp_dir = snid_dir+'all_templates/'
snid_dir = os.getcwd()
# for now just take results from best phase for each SN
temp_list = [x for x in os.listdir(snid_temp_dir)
             if x.split('.')[-1] == 'lnw']
temp_list.sort()
temp_sn_list = [x.split('.')[0] for x in temp_list]
temp_sn_array = numpy.array(temp_sn_list)
n_temp = len(temp_list)

# get all sn types
sn_types = {}
for temp_fn in temp_list:
    sn = temp_fn.split('.')[0]
    f1 = open(snid_temp_dir+temp_fn, 'r')
    sn_types[sn] = f1.readlines()[0].split()[7]
    f1.close()

temp_sn_type_array = numpy.array([sn_types[x] for x in temp_sn_list])

#------------------------------------------------------------------------
def retrieve_snid_template_data(sn, epoch):
    template_fn = '%s%s.lnw' % (snid_temp_dir, sn)
    wave, all_flux, epochs = load_snid_templates(template_fn)
    epoch_arr = numpy.array([float(x) for x in epochs])
    e_ind = numpy.nonzero(epoch_arr==epoch)[0][0]
    flux = all_flux[:,e_ind]
    return wave, flux

def clean_sn_spectrum(wave, flux,
                      fft_filter = False):
    nk = 9
    wmin = numpy.min(wave)
    wmax = numpy.max(wave)
    nw = len(wave)
    knots = numpy.array(
        [wave[int(float(nw)*float(i)/float(nk))]
         for i in numpy.arange(1,nk)])
    sn_spline = scipy.interpolate.splrep(wave, flux, k=3,
                                         t=knots)
    pseudo_cont = scipy.interpolate.splev(wave, sn_spline)
    cdiv_flux = flux / pseudo_cont
    # b - apodize at ends
    init_apod_flux = cdiv_flux-1.0
    apod_curve = numpy.ones(nw, dtype='d')
    edge_buffer = 0.05*(wmax-wmin)
    lo_w_inds = numpy.nonzero(wave < wmin+edge_buffer)[0]
    apod_curve[lo_w_inds] = (1.0+numpy.cos(
        numpy.pi*(1.0-(wave[lo_w_inds]-wmin)/edge_buffer)))/2.0
    hi_w_inds = numpy.nonzero(wave > wmax-edge_buffer)[0]
    apod_curve[hi_w_inds] = (1.0+numpy.cos(
        numpy.pi*(1.0-(wmax-wave[hi_w_inds])/edge_buffer)))/2.0
    apod_flux = init_apod_flux*apod_curve
    # final rescaling
    adj_flux = apod_flux
    return adj_flux

#------------------------------------------------------------------------
def load_snid_templates(template_fn):
    f = open(template_fn, 'r')
    template_lines = f.readlines()
    f.close()
    # parse first line, get number of epochs
    lns = template_lines[0].split()
    n_epochs = int(float(lns[0]))
    # parse line 2, figure out number of infomration lines
    lns = template_lines[1].split()
    n_info = int(float(lns[0]))
    # get list of epoch
    lns = template_lines[n_info+2].split()
    epoch_array = numpy.array(lns[1:])
    # set up wavelength and flux arrays
    data_array = numpy.loadtxt(template_fn, skiprows=(n_info+3))
    wave_array = data_array[:,0]
    flux_array = data_array[:,1:]
    return wave_array, flux_array, epoch_array

def compare_data_to_snid_template(data_wave,
                                  data_flux,
                                  template_fn,
                                  zmax = 0.1):
    #------------------------------------
    # 1 - load in template spectrum
    my_w, my_f, my_e = load_snid_templates(template_fn)
    N=len(my_w)
    # set up fourier filter
    onek_fft_filt = numpy.ones(N,dtype='d')
    kmin = 3
    ktrans = 50
    kmax = 100
    onek_fft_filt[:kmin] = 0.0
    onek_fft_filt[kmax:] = 0.0
    trans_inds = numpy.arange(ktrans, kmax)
    trans_vals = numpy.arange(ktrans, kmax, dtype='d')
    onek_fft_filt[trans_inds] = (1.0+numpy.cos(
        numpy.pi*(1.0-(kmax-trans_vals)/(kmax-ktrans))))/2.0
    # set up z array
    z_array = my_w/my_w[0]-1.0
    good_inds = numpy.nonzero(z_array <= zmax)[0]
    good_z = z_array[good_inds]
    #------------------------------------
    # 2 - clean up input data
    # a - fit spline pseudo-continuum
    nk = 9
    wmin = numpy.min(data_wave)
    wmax = numpy.max(data_wave)
    nw = len(data_wave)
    knots = numpy.array(
        [data_wave[int(float(nw)*float(i)/float(nk))]
         for i in numpy.arange(1,nk)])
    sn_spline = scipy.interpolate.splrep(data_wave, data_flux, k=3,
                                         t=knots)
    pseudo_cont = scipy.interpolate.splev(data_wave, sn_spline)
    cdiv_flux = data_flux / pseudo_cont
    # b - apodize at ends
    init_apod_flux = cdiv_flux-1.0
    apod_curve = numpy.ones(nw, dtype='d')
    edge_buffer = 0.05*(wmax-wmin)
    lo_w_inds = numpy.nonzero(data_wave < wmin+edge_buffer)[0]
    apod_curve[lo_w_inds] = (1.0+numpy.cos(
        numpy.pi*(1.0-(data_wave[lo_w_inds]-wmin)/edge_buffer)))/2.0
    hi_w_inds = numpy.nonzero(data_wave > wmax-edge_buffer)[0]
    apod_curve[hi_w_inds] = (1.0+numpy.cos(
        numpy.pi*(1.0-(wmax-data_wave[hi_w_inds])/edge_buffer)))/2.0
    apod_flux = init_apod_flux*apod_curve
    # c - bin in same log-lambda scale as template
    flux_interp = scipy.interpolate.interp1d(data_wave, apod_flux,
                                             bounds_error=False,
                                             fill_value=0.0)
    int_flux = flux_interp(my_w)
    sigma1 = numpy.mean(int_flux**2)**0.5
    # d - take fourier transform
    init_fft = numpy.fft.fft(int_flux)
    #------------------------------------
    # 3 - cross-correlate it with template
    n_epochs = len(my_e)
    zbest_results = numpy.zeros(n_epochs, dtype='d')
    dz_results    = numpy.zeros(n_epochs, dtype='d')
    r_results     = numpy.zeros(n_epochs, dtype='d')
    lap_results   = numpy.zeros(n_epochs, dtype='d')
    rlap_results  = numpy.zeros(n_epochs, dtype='d')
    for i in range(n_epochs):
        curr_flux = my_f[:,i]
        sigma2 = numpy.mean(curr_flux**2)**0.5
        curr_fft = numpy.fft.fft(curr_flux)
        curr_corr = (init_fft*curr_fft.conjugate())/(
            sigma1*sigma2*float(N))
        corr_values = numpy.fft.ifft(curr_corr)
        good_xc = corr_values[good_inds]
        # find height and redshift of best peak
        h = numpy.max(good_xc.real)
        ibest = good_xc.argmax()
        zbest = good_z[ibest]
        # get peak width and thus z error
        zhi = good_z[numpy.nonzero(good_xc > 0.5*h)[0][-1]]
        init_zlo = good_z[numpy.nonzero(good_xc > 0.5*h)[0][0]]
        if init_zlo == 0.0:
            zlo = zbest-(zhi-zbest)
        else:
            zlo = init_zlo
        w = (zhi-zlo)
        # get correlation strength r from antisymmetric part
        xc_sym = numpy.zeros(N, dtype='d')
        xc_sym[N/2-ibest:] = corr_values[:N/2+ibest]
        xc_sym[:N/2-ibest] = corr_values[N/2+ibest:]
        xc_asym = 0.5*(xc_sym - xc_sym[::-1])
        siga = numpy.mean(2.0*(xc_asym.real)**2)**0.5
        r = h/((2.0**0.5)*siga)
        # given zbest, calculate lap
        rest_wmin = wmin/(1.0+zbest)
        rest_wmax = wmax/(1.0+zbest)
        nonzero_temp_wave = my_w[numpy.nonzero(curr_flux > 0.0)[0]]
        temp_wmin = numpy.min(nonzero_temp_wave)
        temp_wmax = numpy.max(nonzero_temp_wave)
        obs_wmin = max(rest_wmin, temp_wmin)
        obs_wmax = min(rest_wmax, temp_wmax)
        lap = numpy.log(obs_wmax/obs_wmin)
        # finally z error
        #zbest_err = float(ibest)/(1.0+r*lap)
        zbest_err = 2.0*w/(1.0+r*lap)
        # SAVE RESULTS
        zbest_results[i] = zbest
        dz_results[i]    = zbest_err
        r_results[i]     = r
        lap_results[i]   = lap
        rlap_results[i]  = r*lap
    return [my_e,
            zbest_results,
            dz_results,
            r_results,
            lap_results,
            rlap_results]

def compare_all_snid(data_wave, data_flux, zmax=0.1):
    final_epochs = numpy.zeros(n_temp, dtype='d')
    final_zbests = numpy.zeros(n_temp, dtype='d')
    final_zerrs  = numpy.zeros(n_temp, dtype='d')
    final_rlaps  = numpy.zeros(n_temp, dtype='d')    
    for i in range(n_temp):
        temp_fn = temp_list[i]
        sn_name = temp_fn.split('.')[0]
        curr_results = compare_data_to_snid_template(data_wave,
                                                     data_flux,
                                                     snid_temp_dir+temp_fn,
                                                     zmax = zmax)
        curr_epochs = curr_results[0]
        curr_zbests = curr_results[1]
        curr_zerrs  = curr_results[2]
        curr_rlaps  = curr_results[5]
        # find the best epoch
        ibest = curr_rlaps.argmax()
        final_epochs[i] = curr_epochs[ibest]
        final_zbests[i] = curr_zbests[ibest]
        final_zerrs[i]  = curr_zerrs[ibest]
        final_rlaps[i]  = curr_rlaps[ibest]
    return [temp_sn_array,
            temp_sn_type_array,
            final_epochs,
            final_zbests,
            final_zerrs,
            final_rlaps]

"""
# print out results order by rlap
rlap_order = final_rlaps.argsort()[::-1]
for q in range(n_temp):
    i = rlap_order[q]
    sn = temp_list[i].split('.')[0]
    print '%s %s %5.01f %7.04f %6.04f %7.04f' % (
        sn.ljust(15),
        sn_types[sn].ljust(10),
        final_epochs[i],
        final_zbests[i],
        final_zerrs[i],
        final_rlaps[i])
          
# get final weights
weights = numpy.zeros(n_temp, dtype='d')
weights[numpy.nonzero(final_rlaps >= 4.0)[0]] = 1.0
weights[numpy.nonzero(final_rlaps >= 5.0)[0]] = 3.0
weights[numpy.nonzero(final_rlaps >= 6.0)[0]] = 5.0
wt_sum = numpy.sum(weights)
# final weighted values
weighted_zbest = numpy.sum(weights*final_zbests)/wt_sum
#weighted_zerr  = numpy.sum(weights*final_zerrs)/wt_sum
weighted_zerr  = (numpy.sum(weights/final_zerrs**2)/wt_sum)**-0.5
weighted_epoch = numpy.sum(weights*final_epochs)/wt_sum

print weighted_zbest, weighted_zerr, weighted_epoch
"""
