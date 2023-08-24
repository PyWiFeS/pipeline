"""
Determine radial velocities
How is template determined?
"""

# PyWiFeS tools
import process_stellar as ps

ps.calc_rv_template(spect,wave,sig, template_dir,bad_intervals,smooth_distance=101, \
    gaussian_offset=1e-4,nwave_log=1e4,oversamp=1,fig_fn='',convolve_template=True,\
    starnumber=0, plotit=False, save_figures=False, save_dir='./', heliocentric_correction=0.)


    """Compute a radial velocity based on an best fitting template spectrum.
    Teff is estimated at the same time.
    
    Parameters
    ----------
    spect: array-like
        The reduced WiFeS spectrum
        
    wave: array-like
        The wavelengths corresponding to the reduced WiFeS spectrum
        
    template_conv_dir: string
        The directory containing template spectra convolved to 0.1 Angstrom resolution
        
    bad_intervals: 
        List of wavelength intervals where e.g. telluric absorption is bad.
        
    smooth_distance: float
        Distance to smooth for "continuum" correction
        
    oversamp: float
        Oversampling of the input wavelength scale. The slit is assumed 2 pixels wide.
    
    gaussian_offset: float
        Offset for the likelihood function from a Gaussian normalised to 1. 
        
    Returns
    -------
    rv: float
        Radial velocity in km/s
    rv_sig: float
        Uncertainty in radial velocity (NB assumes good model fit)
    temp: int
        Temperature of model spectrum used for cross-correlation.
    """
