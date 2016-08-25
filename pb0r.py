from __future__ import division

import numpy as np
import datetime
from sunpy.time import parse_time, julian_day
from sunpy.wcs import convert_hcc_hg, convert_hg_hcc
from sunpy.sun import constants

def pb0r(date, stereo=None, soho=False, arcsec=False):
    """To calculate the solar P, B0 angles and the semi-diameter.

    Parameters
    -----------
    date: a date/time object

    stereo: { 'A' | 'B' | None }
        calculate the solar P, B0 angles and the semi-diameter from the point
        of view of either of the STEREO spacecraft.

    soho: { False | True }
        calculate the solar P, B0 angles and the semi-diameter from the point
        of view of the SOHO spacecraft.  SOHO sits at the Lagrange L1 point
        which is about 1% closer to the Sun than the Earth.  Implementation
        of this seems to require the ability to read SOHO orbit files.

    arcsec: { False | True }
        return the semi-diameter in arcseconds.

    Returns:
    -------
    A dictionary with the following keys with the following meanings:

    p  -  Solar P (position angle of pole)  (degrees)
    b0 -  latitude of point at disk centre (degrees)
    sd -  semi-diameter of the solar disk in arcminutes

    See Also
    --------
    IDL code equavalent:
        http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/pb0r.pro
    """
    if (stereo is not None) and soho:
        raise ValueError("Cannot set STEREO and SOHO simultaneously")

    # place holder for STEREO calculation
    if stereo is not None:
        raise ValueError("STEREO solar P, B0 and semi-diameter calcution" + \
                         " is not supported.")

    # number of Julian days since 2415020.0
    de = julian_day(date) - 2415020.0

    # get the longitude of the sun etc.
    sun_position = sun_pos(date)
    longmed = sun_position["longitude"]
    #ra = sun_position["ra"]
    #dec = sun_position["dec"]
    appl = sun_position["app_long"]
    oblt = sun_position["obliq"]

    # form the aberrated longitude
    Lambda = longmed - (20.50 / 3600.0)

    # form longitude of ascending node of sun's equator on ecliptic
    node = 73.6666660 + (50.250 / 3600.0) * ((de / 365.250) + 50.0)
    arg = Lambda - node

    # calculate P, the position angle of the pole
    p = np.rad2deg(\
        np.arctan(-np.tan(np.deg2rad(oblt)) * np.cos(np.deg2rad(appl))) + \
        np.arctan(-0.127220 * np.cos(np.deg2rad(arg))))

    # B0 the tilt of the axis...
    b = np.rad2deg(np.arcsin(0.12620 * np.sin(np.deg2rad(arg))))

    # ... and the semi-diameter
    # Form the mean anomalies of Venus(MV),Earth(ME),Mars(MM),Jupiter(MJ)
    # and the mean elongation of the Moon from the Sun(D).
    t = de / 36525.0
    mv = 212.60 + np.mod(58517.80 * t, 360.0)
    me = 358.4760 + np.mod(35999.04980 * t, 360.0)
    mm = 319.50 + np.mod(19139.860 * t, 360.0)
    mj = 225.30 + np.mod(3034.690 * t, 360.0)
    d = 350.70 + np.mod(445267.110 * t, 360.0)

    # Form the geocentric distance(r) and semi-diameter(sd)
    r = 1.0001410 - (0.0167480 - 0.00004180 * t) * np.cos(np.deg2rad(me)) \
        - 0.000140 * np.cos(np.deg2rad(2.0 * me)) \
        + 0.0000160 * np.cos(np.deg2rad(58.30 + 2.0 * mv - 2.0 * me)) \
        + 0.0000050 * np.cos(np.deg2rad(209.10 + mv - me)) \
        + 0.0000050 * np.cos(np.deg2rad(253.80 - 2.0 * mm + 2.0 * me)) \
        + 0.0000160 * np.cos(np.deg2rad(89.50 - mj + me)) \
        + 0.0000090 * np.cos(np.deg2rad(357.10 - 2.0 * mj + 2.0 * me)) \
        + 0.0000310 * np.cos(np.deg2rad(d))

    sd_const = constants.radius / constants.au
    sd = np.arcsin(sd_const / r) * 10800.0 / np.pi

    # place holder for SOHO correction
    if soho:
        raise ValueError("SOHO correction (on the order of 1% " + \
                        "since SOHO sets at L1) not yet supported.")

    if arcsec:
        return {"p": p, "b0": b, "sd": sd * 60.0}
    else:
        return {"p": p, "b0": b, "sd": sd, "l0": 0.0}

def sun_pos(date, is_julian=False, since_2415020=False):
    """ Calculate solar ephemeris parameters.  Allows for planetary and lunar
    perturbations in the calculation of solar longitude at date and various
    other solar positional parameters. This routine is a truncated version of
    Newcomb's Sun and is designed to give apparent angular coordinates (T.E.D)
    to a precision of one second of time.

    Parameters
    -----------
    date: a date/time object or a fractional number of days since JD 2415020.0

    is_julian: { False | True }
        notify this routine that the variable "date" is a Julian date
        (a floating point number)

    since_2415020: { False | True }
        notify this routine that the variable "date" has been corrected for
        the required time offset

    Returns:
    -------
    A dictionary with the following keys with the following meanings:

    longitude  -  Longitude of sun for mean equinox of date (degs)
    ra         -  Apparent RA for true equinox of date (degs)
    dec        -  Apparent declination for true equinox of date (degs)
    app_long   -  Apparent longitude (degs)
    obliq      -  True obliquity (degs)longditude_delta:

    See Also
    --------
    IDL code equavalent:
        http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/sun_pos.pro

    Examples
    --------
    >>> sp = sun_pos('2013-03-27')

    """
    # check the time input
    if is_julian:
        # if a Julian date is being passed in
        if since_2415020:
            dd = date
        else:
            dd = date - 2415020.0
    else:
        # parse the input time as a julian day
        if since_2415020:
            dd = julian_day(date)
        else:
            dd = julian_day(date) - 2415020.0

    # form time in Julian centuries from 1900.0
    t = dd / 36525.0

    # form sun's mean longitude
    l = (279.6966780 + np.mod(36000.7689250 * t, 360.00)) * 3600.0

    # allow for ellipticity of the orbit (equation of centre) using the Earth's
    # mean anomaly ME
    me = 358.4758440 + np.mod(35999.049750 * t, 360.0)
    ellcor = (6910.10 - 17.20 * t) * np.sin(np.deg2rad(me)) + \
    72.30 * np.sin(np.deg2rad(2.0 * me))
    l = l + ellcor

    # allow for the Venus perturbations using the mean anomaly of Venus MV
    mv = 212.603219 + np.mod(58517.8038750 * t, 360.0)
    vencorr = 4.80 * np.cos(np.deg2rad(299.10170 + mv - me)) + \
          5.50 * np.cos(np.deg2rad(148.31330 + 2.0 * mv - 2.0 * me)) + \
          2.50 * np.cos(np.deg2rad(315.94330 + 2.0 * mv - 3.0 * me)) + \
          1.60 * np.cos(np.deg2rad(345.25330 + 3.0 * mv - 4.0 * me)) + \
          1.00 * np.cos(np.deg2rad(318.150 + 3.0 * mv - 5.0 * me))
    l = l + vencorr

    # Allow for the Mars perturbations using the mean anomaly of Mars MM
    mm = 319.5294250 + np.mod(19139.858500 * t, 360.0)
    marscorr = 2.0 * np.cos(np.deg2rad(343.88830 - 2.0 * mm + 2.0 * me)) + \
            1.80 * np.cos(np.deg2rad(200.40170 - 2.0 * mm + me))
    l = l + marscorr

    # Allow for the Jupiter perturbations using the mean anomaly of Jupiter MJ
    mj = 225.3283280 + np.mod(3034.69202390 * t, 360.00)
    jupcorr = 7.20 * np.cos(np.deg2rad(179.53170 - mj + me)) + \
          2.60 * np.cos(np.deg2rad(263.21670 - mj)) + \
          2.70 * np.cos(np.deg2rad(87.14500 - 2.0 * mj + 2.0 * me)) + \
          1.60 * np.cos(np.deg2rad(109.49330 - 2.0 * mj + me))
    l = l + jupcorr

    # Allow for the Moons perturbations using the mean elongation of the Moon
    # from the Sun D
    d = 350.73768140 + np.mod(445267.114220 * t, 360.0)
    mooncorr = 6.50 * np.sin(np.deg2rad(d))
    l = l + mooncorr

    # Note the original code is
    # longterm  = + 6.4d0 * sin(( 231.19d0  +  20.20d0 * t )*!dtor)
    longterm = 6.40 * np.sin(np.deg2rad(231.190 + 20.20 * t))
    l = l + longterm
    l = np.mod(l + 2592000.0, 1296000.0)
    longmed = l / 3600.0

    # Allow for Aberration
    l = l - 20.5

    # Allow for Nutation using the longitude of the Moons mean node OMEGA
    omega = 259.1832750 - np.mod(1934.1420080 * t, 360.0)
    l = l - 17.20 * np.sin(np.deg2rad(omega))

    # Form the True Obliquity
    oblt = 23.4522940 - 0.01301250 * t + \
    (9.20 * np.cos(np.deg2rad(omega))) / 3600.0

    # Form Right Ascension and Declination
    l = l / 3600.0
    ra = np.rad2deg(np.arctan2(np.sin(np.deg2rad(l)) * \
                        np.cos(np.deg2rad(oblt)), np.cos(np.deg2rad(l))))

    if isinstance(ra, np.ndarray):
        ra[ra < 0.0] += 360.0
    elif ra < 0.0:
        ra = ra + 360.0

    dec = np.rad2deg(np.arcsin(np.sin(np.deg2rad(l)) * \
                                np.sin(np.deg2rad(oblt))))

    # convert the internal variables to those listed in the top of the
    # comment section in this code and in the original IDL code.
    return {"longitude": longmed, "ra": ra, "dec": dec, "app_long": l,
            "obliq": oblt}
