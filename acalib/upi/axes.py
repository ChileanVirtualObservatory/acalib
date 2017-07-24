import numpy as np
import astropy.units as u
from astropy.nddata import support_nddata
from astropy import log

from acalib import core

# axes_names used in formatting
from acalib.upi.data import Data


@support_nddata
def axes_names(data, wcs=None):
    """
        Get the axes's names.

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.

        Returns
        -------
        result: numpy.ndarray
            Numpy ndarray with the axes's names from the WCS.

    """
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    return np.array(wcs.axis_type_names)[::-1]


from acalib.upi.formatting import _unitize, _world_table_creator


@support_nddata
def extent(data, wcs=None, region=None):
    """
        Get the axes extent.

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.
        region :(lower : (M,N) or (M,N,Z), upper : (M,N) or (M,N,Z))
            Start and End index in data (int tuples)

        Returns
        -------
        result: (M, N) tuple of astropy.units.quantity.Quantity
            Axes extent

    """
    # TODO: These can be a decorator
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None

    if region is None:
        lower=None
        upper=None
    else:
        (lower,upper) = region

    if lower == None:
        lower = np.zeros(data.ndim)
    if upper == None:
        upper = data.shape
    idx = [lower, upper]
    idx_f = np.fliplr(idx)
    values = wcs.wcs_pix2world(idx_f, 0)
    values = np.fliplr(values)
    return (_unitize(values[0], wcs), _unitize(values[1], wcs))


@support_nddata
def center(data, wcs=None):
    """
        Get center of the data

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.

        Returns
        -------
        result: astropy.units.quantity.Quantity
            Center of the data

    """
    # TODO: These can be a decorator
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    val = wcs.wcs.crval[::-1]
    return _unitize(val, wcs)


@support_nddata
def axes_units(data, wcs=None):
    """
        Get units of the axes

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.

        Returns
        -------
        result: (M,N) or (M,N,Z) numpy.ndarray
            Vector with the units of the axes

    """
    # TODO: These can be a decorator (
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    uvec = np.array(wcs.wcs.cunit)[::-1]
    return uvec


@support_nddata
def resolution(data, wcs=None):
    """
        Get the resolution of data

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.

        Returns
        -------
        result: (M,N) or (M,N,Z) numpy.ndarray
            Resolution of the data

    """
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    val = wcs.wcs.cdelt[::-1]
    return _unitize(val, wcs)


@support_nddata
def spectral_velocities(data, wcs=None, fqs=None, fqis=None, restfrq=None):
    """
        Get the spectral velocities from frequencies fqs given a rest
        frequency (by default search for it in the WCS). If fqs is None,
        then frequencies indices (fqis) need to be given.

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.
        fqs : astropy.units.quantity.Quantity
            Array of frequencies with units.
        fqis : list of integers
            Array of frequencies indices
        restfrq : astropy.units.quantity.Quantity
            Rest frequency

        Returns
        -------
        result: astropy.units.quantity.Quantity
            Array of Spectral velocities.


    """
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    if data.ndim != 3:
        log.error("Not spectral axis found.")
        return None
    if restfrq is None:
        restfrq = wcs.wcs.restfrq * u.Hz
    if fqs is None:
        dim = wcs.wcs.spec
        if fqis is None:
            # Semi Hardconded...
            fqis = np.arange(data.shape[data.ndim - dim - 1])
        idx = np.zeros((fqis.size, data.ndim))
        idx[:, dim] = fqis
        vals = wcs.all_pix2world(idx, 0)
        fqs = vals[:, dim] * u.Hz
    eq = u.doppler_radio(restfrq)
    return fqs.to(u.km / u.s, equivalencies=eq)


@support_nddata
def features(data, wcs=None, region=None):
    """
        Creates an array with WCS axea in features format

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.
        region :(lower : (M,N) or (M,N,Z), upper : (M,N) or (M,N,Z))
            Start and End index in data (int tuples)

        Returns
        -------
        result: astropy.table.Table
            Table with WCS information of a section from the data.

    """
    if region  is None:
        lower = None
        upper = None
    else:
        (lower, upper) = region
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    ii = core.index_features(data, lower, upper)
    f = wcs.wcs_pix2world(ii.T, 0)
    return _world_table_creator(f, wcs)

def _get_axis(wcs,name):
    return wcs.naxis - wcs.axis_type_names.index(name) - 1

def _is_spectra(data,wcs):
    try:
        freq_axis=_get_axis(wcs,"FREQ")
    except ValueError:
        return False
    for i in range(wcs.naxis):
        if i == freq_axis:
            continue
        if data.shape[i]!= 1:
            return False
    return True

@support_nddata
# TODO: Consider using "box" structure rather than up and low
def opening(data, center, window, wcs=None):
    """
        Field of view (center +- window) converted to indices

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        center : astropy.units.quantity.Quantity
            Center of the field of view in WCS.
        window : astropy.units.quantity.Quantity
            Window for the field in WCS.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.
        Returns
        -------
        result: ((M1,N1,Z1),(M2,N2,Z2)) tuple of tuple of ints
    """
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    off_low = center - window
    off_low = np.array([x.value for x in off_low])
    off_up = center + window
    off_up = np.array([x.value for x in off_up])
    ld = np.rint(wcs.wcs_world2pix([off_low[::-1]], 0))
    lu = np.rint(wcs.wcs_world2pix([off_up[::-1]], 0))
    lower = np.array([ld, lu]).min(axis=0)
    upper = np.array([ld, lu]).max(axis=0)
    lower = core.fix_limits(data, lower[0][::-1])
    upper = core.fix_limits(data, upper[0][::-1])
    return lower, upper
