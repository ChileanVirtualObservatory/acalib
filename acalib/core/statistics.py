def rms(data, mask=None):
    """
    Compute the RMS of data. If mask != None, then we use that mask.

    Parameters
    ----------
    data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData or or astropy.nddata.NDDataRef
        Astronomical data cube.

    mask : numpy.ndarray (default = None)

    Returns
    -------
    RMS of the data (float)
    """
    # TODO: check photutils background estimation for using that if possible
    if mask is not None:
        data = fix_mask(data, mask)
    mm = data * data
    rms = np.sqrt(mm.sum() * 1.0 / mm.size)
    return rms

def snr_estimation(data, mask=None, noise=None, points=1000, full_output=False):
    """
    Heurustic that uses the inflexion point of the thresholded RMS to estimate where signal is dominant w.r.t. noise

    Parameters
    ----------
    data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
        Astronomical data cube.

    mask : numpy.ndarray (default = None)

    noise : float (default=None)
        Noise level, if not given will use rms of the data.

    points : (default=1000)

    full_output : boolean (default=False)
        Gives verbose results if True

    Returns
    --------

    "Signal to Noise Radio" value

    """
    if noise is None:
        noise = rms(data, mask)
    x = []
    y = []
    n = []
    sdata = data[data > noise]
    for i in range(1, int(points)):
        val = 1.0 + 2.0 * i / points
        sdata = sdata[sdata > val * noise]
        if sdata.size < 2:
            break
        n.append(sdata.size)
        yval = sdata.mean() / noise
        x.append(val)
        y.append(yval)
    y = np.array(y)
    v = y[1:] - y[0:-1]
    p = v.argmax() + 1
    snrlimit = x[p]
    if full_output == True:
        return snrlimit, noise, x, y, v, n, p
    return snrlimit

