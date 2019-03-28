from __future__ import division, print_function
import numpy
import scipy.interpolate

#-----------------------------------------------------------------------------
def blkrep(arr, mx, my):
    # if arr is a nx by ny array, make a mx*nx by my*ny array duplicating arr
    ny, nx = numpy.shape(arr)
    newxarr = numpy.reshape(
        (arr.T*numpy.ones([mx,nx,ny],dtype='f')).T.flatten(),
        [ny,mx*nx])
    outarr = numpy.reshape(
        (newxarr*numpy.ones([my,ny,nx*mx],dtype='f')).T.flatten(),
        [mx*nx,my*ny]).T    
    return outarr

def blkavg(arr, mx, my):
    nby, nbx = numpy.shape(arr)
    ny = nby//my
    nx = nbx//my
    # average in y
    y1 = numpy.sum(numpy.reshape(arr.flatten(), [ny, my, nbx]),
                   axis=1)/float(my)
    x1 = numpy.sum(numpy.reshape(y1.T.flatten(), [nx, mx, ny]),
                   axis=1).T/float(mx)
    return x1

#-----------------------------------------------------------------------------
def transform_data(data, wave,
                   return_lambda=False,
                   out_lambda=None):
    # get dimensions...
    ny, nx = numpy.shape(data)
    xarr = numpy.arange(nx)
    yarr = numpy.arange(ny)
    full_x, full_y = numpy.meshgrid(xarr, yarr)
    # figure out wavelength coverage
    wmin = numpy.min(numpy.min(wave,axis=1))
    wmax = numpy.max(numpy.max(wave,axis=1))
    dw = numpy.abs(wave[:,1:]-wave[:,:-1])
    full_dw = numpy.zeros(numpy.shape(wave))
    full_dw[:,1:] = dw
    full_dw[:,0] = dw[:,0]
    # rectify to uniform wavelength array
    if out_lambda == None:
        wdisp = numpy.mean(dw)
        if wdisp < 0:
            out_lambda = numpy.arange(wmin, wmax-wdisp, -wdisp)[::-1]
        else:
            out_lambda = numpy.arange(wmin, wmax+wdisp, wdisp)
    else:
        # NOTE: BREAKS IF YOU HAVE NON-REGULAR WAVELENGTH SPACING
        wdisp = out_lambda[1]-out_lambda[0]
    #out_lambda_full, out_y_full = numpy.meshgrid(out_lambda, yarr)
    scaled_data = data / full_dw
    # new shortened version
    scaled_interp_data = numpy.zeros([ny, len(out_lambda)], dtype='d')
    for i in range(ny):
        init_curr_data = scaled_data[i,:]
        init_curr_wave = wave[i,:]
        sort_order = init_curr_wave.argsort()
        curr_wave = init_curr_wave[sort_order]
        curr_data = init_curr_data[sort_order]
        curr_interp = scipy.interpolate.interp1d(
            curr_wave, curr_data,
            bounds_error=False, fill_value=0.0)
        new_data = curr_interp(out_lambda)
        scaled_interp_data[i,:] = new_data
    #scaled_interp_data = scipy.interpolate.griddata(
    #    (wave.flatten(), full_y.flatten()),
    #    scaled_data.flatten(),
    #    (out_lambda_full, out_y_full),
    #    method='linear')
    interp_data = scaled_interp_data*wdisp
    bad_inds = numpy.nonzero(interp_data != interp_data)
    interp_data[bad_inds] = 0.0
    if return_lambda:
        return interp_data, out_lambda
    else:
        return interp_data

def detransform_data(new_data, orig_data, wave):
    # get dimensions...
    ny, nx = numpy.shape(orig_data)
    xarr = numpy.arange(nx)
    yarr = numpy.arange(ny)
    full_x, full_y = numpy.meshgrid(xarr, yarr)
    # figure out wavelength coverage
    wmin = numpy.min(numpy.min(wave,axis=1))
    wmax = numpy.max(numpy.max(wave,axis=1))
    dw = numpy.abs(wave[:,1:]-wave[:,:-1])
    full_dw = numpy.zeros(numpy.shape(wave))
    full_dw[:,1:] = dw
    full_dw[:,0] = dw[:,0]
    wdisp = numpy.mean(dw)
    # rectify to uniform wavelength array
    if wdisp < 0:
        out_lambda = numpy.arange(wmin, wmax-wdisp, -wdisp)[::-1]
    else:
        out_lambda = numpy.arange(wmin, wmax+wdisp, wdisp)
    #out_lambda_full, out_y_full = numpy.meshgrid(out_lambda, yarr)
    # new quick version
    scaled_interp_data = numpy.zeros([ny, len(wave[0,:])], dtype='d')
    for i in range(ny):
        curr_data = new_data[i,:]
        curr_wave = wave[i,:]
        curr_interp = scipy.interpolate.interp1d(
            out_lambda, curr_data,
            bounds_error=False, fill_value=0.0)
        fixed_data = curr_interp(curr_wave)
        scaled_interp_data[i,:] = fixed_data
    #scaled_interp_data = scipy.interpolate.griddata(
    #    (out_lambda_full.flatten(), out_y_full.flatten()),
    #    new_data.flatten(),
    #    (wave, full_y),
    #    method='linear')
    interp_data = scaled_interp_data*(full_dw/wdisp)
    bad_inds = numpy.nonzero(interp_data != interp_data)
    interp_data[bad_inds] = 0.0
    return interp_data

#-----------------------------------------------------------------------------
