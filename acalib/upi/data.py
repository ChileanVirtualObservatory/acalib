from astropy import nddata as ndd

from acalib import upi
from acalib import io


class Data(ndd.NDDataRef):
    """
    A generic represenation of astronomical n-dimensional data array. Extends NDData.
   
    """

    def axes_names(self):
        """  Get the axes's names.
        Returns
        -------
        result: numpy.ndarray
            Numpy ndarray with the axes's names from the WCS.
        """

        return upi.axes_names(self)

    def cut(self, region=None):
        """
        Get a cut of the cube.

        Parameters
        ----------
        region :(lower : (M,N) or (M,N,Z), upper : (M,N) or (M,N,Z))
            Start and End index in data (int tuples)
        Returns
        -------
        result: acalib.upi.AData.
            AData cut from lower to upper.
        """
        return upi.reduction.cut(self, region=region)

    def extent(self, region=None):
        """
        Get the axes extent.

        Parameters
        ----------
    region :(lower : (M,N) or (M,N,Z), upper : (M,N) or (M,N,Z))
            Start and End index in data (int tuples)

        Returns
        -------
        result: (M, N) tuple of astropy.units.quantity.Quantity
            Axes extent
        """
        return upi.axes.extent(self, region)

    def center(self):
        """
        Get center of the data

        Returns
        -------
        result: astropy.units.quantity.Quantity
            Center of the data
        """
        return upi.axes.center(self)

    def axes_units(self):
        """
        Get units of the axes

        Returns
        -------
        result: (M,N) or (M,N,Z) numpy.ndarray
            Vector with the units of the axes
        """
        return upi.axes.axes_units(self)

    def resolution(self):
        """
        Get the resolution of data

        Returns
        -------
        result: (M,N) or (M,N,Z) numpy.ndarray
            Resolution of the data
        """
        return upi.axes.resolution(self)

    def spectral_velocities(self, fqs=None, fqis=None, restfrq=None):
        """
        Get the spectral velocities from frequencies fqs given a rest
        frequency (by default search for it in the WCS). If fqs is None,
        then frequencies indices (fqis) need to be given.

        Parameters
        ----------
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
        return upi.axes.spectral_velocities(self, fqs=fqs, fqis=fqis, restfrq=restfrq)

    def features(self, region=None):
        """
        Creates an array with WCS axea in features format

        Parameters
        ----------
        region :(lower : (M,N) or (M,N,Z), upper : (M,N) or (M,N,Z))
            Start and End index in data (int tuples)

        Returns
        -------
        result: astropy.table.Table
            Table with WCS information of a section from the data.
        """
        return upi.axes.features(self, region=region)

    def opening(self, center, window):
        """
        Field of view (center +- window) converted to indices

        Parameters
        ----------
        center : astropy.units.quantity.Quantity
            Center of the field of view in WCS.
        window : astropy.units.quantity.Quantity
            Window for the field in WCS.

        Returns
        -------
        result: ((M1,N1,Z1),(M2,N2,Z2)) tuple of tuple of ints
        """
        return upi.axes.opening(self, center=center, window=window)

    def rms(self):
        """
        Compute the RMS of data.

        Returns
        -------
        rms : float
            RMS of data
        """
        return upi.flux.rms(self)

    def moment0(self):
        """
        Calculate moment 0 from a data cube.

        Returns
        -------
        result: astropy.nddata.AData
            Moment 0 of the data cube
        """
        return upi.reduction.moment0(self, restfrq=None)

    def moment1(self, restfrq=None):
        """
        Calculate moment 1 from a data cube.

        restfrq : astropy.units.quantity.Quantity
            Rest frequency

        Returns
        -------
        result: astropy.nddata.NDData
            Moment 1 of the data cube
        """
        return upi.reduction.moment1(self, restfrq=restfrq)

    def moment2(self, restfrq=None):
        """
        Calculate moment 2 from a data cube.

        restfrq : astropy.units.quantity.Quantity
            Rest frequency

        Returns
        -------
        result: astropy.nddata.AData
            Moment 2 of the data cube
        """
        return upi.reduction.moment2(self, restfrq=restfrq)

    def visualize(self):
        """
        Generic function to visualize data, line-plot for 1D and image for 2D.
        """
        io.visualize(self)

    def visualize_image(self, contour=False,cmap = None):
        io.visualize_image(self, contour=contour,cmap=cmap)

    def visualize_spectra(self, velocities=False):
        io.visualize_spectra(self, velocities=velocities)

    def select_region(self, ra_1=None, dec_1=None, ra_2=None, dec_2=None, interactive=False):
        upi.reduction.select_region(self, ra_1=ra_1, dec_1=dec_1, ra_2=ra_2, dec_2=dec_2, interactive=interactive)

    def select_band(self,freq_1=None,freq_2=None,interactive=False):
        upi.reduction.select_band(self, freq_1=freq_1,freq_2=freq_2, interactive=interactive)
