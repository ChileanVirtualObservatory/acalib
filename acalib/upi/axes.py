import numpy as np
import astropy.units as u
from astropy.nddata import support_nddata, NDDataRef
from astropy import log

from acalib import core

# axes_names used in formatting
@support_nddata
def axes_names(data,wcs=None):
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
def cut(data, wcs=None, mask=None, unit=None, lower=None, upper=None):
    """
        Get the axes's names.

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.
        mask : numpy.ndarray
            Mask for data.
        unit : astropy.units.Unit
            Astropy unit (http://docs.astropy.org/en/stable/units/).
        lower : tuple
            Start coordinate from where to cut.
        upper : tuple
            Coordinate to end cut.

        Returns
        -------
        result: astropy.nddata.NDDataRef.
            data cut from lower to upper.

    """
    # Check for NDDataSlicing... maybe this is already done by astropy.nddata package.
    mslab = core.slab(data, lower, upper)
    scube = data[mslab]
    newwcs = wcs.slice(mslab, numpy_order=True)
    return NDDataRef(scube, wcs=newwcs, unit=unit)

@support_nddata
def extent(data,wcs=None,lower=None,upper=None):
    """
        Get the axes extent.

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.
        lower : (M,N) or (M,N,Z) tuple of int
            Start coordinate in data
        upper : (M,N) or (M,N,Z) tuple of int
            End coordinate in data

        Returns
        -------
        result: (M, N) tuple of astropy.units.quantity.Quantity
            Axes extent

    """
    #TODO: These can be a decorator
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    if lower==None:
        lower=np.zeros(data.ndim)
    if upper==None:
        upper=data.shape
    idx=[lower,upper]
    idx_f  = np.fliplr(idx)
    values = wcs.wcs_pix2world(idx_f, 0)
    values = np.fliplr(values)
    return (_unitize(values[0],wcs),_unitize(values[1],wcs))

@support_nddata
def center(data,wcs=None):
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
    #TODO: These can be a decorator
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    val=wcs.wcs.crval[::-1]
    return _unitize(val,wcs)

@support_nddata
def axes_units(data,wcs=None):
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
    #TODO: These can be a decorator (
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    uvec=np.array(wcs.wcs.cunit)[::-1]
    return uvec

@support_nddata
def resolution(data,wcs=None):
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
    val=wcs.wcs.cdelt[::-1]
    return _unitize(val,wcs)

@support_nddata
def spectral_velocities(data,wcs=None,fqs=None,fqis=None,restfrq=None):
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
    if restfrq is None:
        restfrq=wcs.wcs.restfrq*u.Hz
    if fqs is None:
        if fqis is None:
           return None
        dim=wcs.wcs.spec
        idx=np.zeros((fqis.size,data.ndim))
        idx[:,dim]=fqis
        vals=wcs.all_pix2world(idx,0)
        fqs=vals[:,dim]*u.Hz
    eq=u.doppler_radio(restfrq)
    return fqs.to(u.km/u.s, equivalencies=eq)


@support_nddata
def features(data,wcs=None,lower=None,upper=None):
    """
        Creates an array with WCS axea in features format

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.
        lower : (M,N) or (M,N,Z) tuple of integers
            Start coordinate in data.
        upper : (M,N) or (M,N,Z) tuple of integers
            End coordinate in data.

        Returns
        -------
        result: astropy.table.Table
            Table with WCS information of a section from the data.

    """
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    ii=core.index_features(data,lower,upper)
    f=wcs.wcs_pix2world(ii.T,0)
    return _world_table_creator(f,wcs)


@support_nddata
# TODO: Consider using "box" structure rather than up and low
def opening(data,center,window,wcs=None):
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
    off_low=center-window
    off_low = np.array([x.value for x in off_low])
    off_up = center+window
    off_up = np.array([x.value for x in off_up])
    #dim = len(center.colnames)
    #off_low = np.array([center[0][i] - window[0][i] for i in range(dim)])
    #off_up  = np.array([center[0][i] + window[0][i] for i in range(dim)])
    ld=np.rint(wcs.wcs_world2pix([off_low[::-1]],0))
    lu=np.rint(wcs.wcs_world2pix([off_up[::-1]],0))
    lower=np.array([ld,lu]).min(axis=0)
    upper=np.array([ld,lu]).max(axis=0)
    lower=core.fix_limits(data,lower[0][::-1])
    upper=core.fix_limits(data,upper[0][::-1])
    return (lower,upper)
    #values=np.vstack((lower,upper))
    #return _pix_table_creator(values,wcs)
