import numpy

"""
Atmospheric differential refraction, based on model of Filippenko (1982PASP...94..715F).
See detailed usage in adr_x_y().
"""


def ha_degrees(ha_str):
    full_hours, mins, secs = ha_str.split(':')
    hours = abs(float(full_hours))
    # below is necessary to deal with -00:XX
    if full_hours[0] == '-':
        sign = -1.0
    else:
        sign = 1.0
    hours_total = sign * (hours
                          + float(mins) / 60.0
                          + float(secs) / 3600.0)
    return 15.0 * hours_total


def dec_dms2dd(dec_str):
    if dec_str[0] == '+':
        sign = 1.0
    else:
        sign = -1.0
    deg, min, sec = dec_str[1:].split(':')
    return sign * (float(deg) + float(min) / 60.0 + float(sec) / 3600.0)


def adr_n1(lam):
    # convert angstroms to microns
    lmic = lam * 1.0e-4
    term1 = 64.328
    term2 = 29498.1 / (146.0 - (1.0 / lmic)**2)
    term3 = 255.4 / (41.0 - (1.0 / lmic)**2)
    return 1.0e-6 * (term1 + term2 + term3)


def adr_f1(p, t):
    term1 = p / 720.883
    term2 = 1.0 + 1.0e-6 * p * (1.049 - 0.0157 * t)
    term3 = 1.0 + 0.003661 * t
    return term1 * term2 / term3


def adr_g1(lam, t, f):
    lmic = lam * 1.0e-4
    term1 = 0.0624 - 0.000680 / (lmic**2)
    term2 = 1.0 + 0.003661 * t
    return term1 * f / term2


def adr_ntot(lam, p, t, f):
    real_n = 1.0e-6 * (adr_n1(lam) * 1.0e6 - adr_g1(lam, t, f))
    return real_n * adr_f1(p, t) + 1.0


def adr_r(lam, secz, p, t, f):
    nlam = adr_ntot(lam, p, t, f)
    tanz = (secz**2 - 1.0)**0.5
    return 206265.0 * nlam * tanz


def adr_eta(ha, lat, dec):
    term1 = numpy.sin(numpy.radians(ha))
    term4 = (numpy.cos(numpy.radians(dec))
             * numpy.tan(numpy.radians(lat))
             - numpy.sin(numpy.radians(dec))
             * numpy.cos(numpy.radians(ha))
             )
    return numpy.arctan2(term1, term4)


def adr_x_y(wavelength_array,
            secz,
            objha,
            objdec,
            tellat,
            teltemp,
            telpres,
            telwvp=8.0,
            telpa=0.0,
            ref_wl=5000.0):

    """
    Atmospheric differential refraction, based on model of Filippenko (1982PASP...94..715F).

    Assumes plane-parallel atmosphere; Barrell (1951) model for dependence of index of refraction
    on temperature, air pressure, and water vapour pressure; and values for index of refraction at
    sea level from Edlen 1953 / Coleman, Bozman, & Meggers 1960.

    Default water vapour pressure of 8 mm of mercury from Filippenko corresponds to _saturation_
    pressure (not actual pressure) at air temperature of ~7 C, typical of average conditions for
    latitude of +/-30 deg and altitude of ~2000 m (Allen 1973), where the air pressure is ~600 mm
    of mercury.

    Inputs:

    - wavelength_array: (vacuum) wavelengths at which to calculate the differential refraction.
    - secz: secant of zenith angle, equivalent to the airmass for moderate zenith angles (<60 deg).
    - objha: hour angle of the target in degrees.
    - objdec: declination of the target in degrees.
    - tellat: telescope latitude in degrees.
    - telltemp: air temperature at telescope in Celsius.
    - telepres: air pressure at telescope in mm of mercury.
    - telwvp: water vapour pressure at telescope in mm of mercury.
        Default 8 mm.
    - telpa: position angle (radians East from North).
        Default 0 radians.
    - ref_wl: reference wavelength for the differential refraction in Angstroms.
        Default 5000 A.

    Outputs:

    - list of [array of x-axis offsets for each wavelength, array of y-axis offsets for each wavelength].
        Units of arcseconds.

    """

    # get adr results
    eta = adr_eta(objha, tellat, objdec)
    obj_eta = eta - telpa
    n_pix = len(wavelength_array)
    adrx = numpy.zeros(n_pix, dtype='f')
    adry = numpy.zeros(n_pix, dtype='f')
    # starting point
    r_set = adr_r(ref_wl, secz, telpres, teltemp, telwvp)
    for i in range(n_pix):
        delta_r = adr_r(wavelength_array[i],
                        secz, telpres, teltemp, telwvp) - r_set
        adry[i] = delta_r * numpy.cos(obj_eta)
        adrx[i] = delta_r * numpy.sin(obj_eta)
    return [adrx, adry]
