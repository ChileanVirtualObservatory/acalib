import numpy as np


def gaussian_function(mu, P, feat, peak):
    """
    Generates an N-dimensional Gaussian using the feature matrix feat,
    centered at mu, with precision matrix P and with intensity peak.

    Parameters
    ----------
    mu : numpy.ndarray
        Centers of gaussians array.
    P : numpy.ndarray
        Precision matrix.
    feat : numpy.ndarray.
        Features matrix.
    peak : float
        Peak value of the resulting evaluation.

    Returns
    -------
    result: 2D numpy.ndarray
        Returns the gaussian function evaluated at the value on feat. 
    """

    cent_feat = np.empty_like(feat)
    for i in range(len(mu)):
        cent_feat[i] = feat[i] - mu[i]
    qform = (P.dot(cent_feat)) * cent_feat
    quad = qform.sum(axis=0)
    res = np.exp(-quad / 2.0)
    res = peak * (res / res.max())
    return res


def create_mould(P, delta):
    """
    Creates a Gaussian mould with precision matrix P, using the already computed values of delta.

    Parameters
    ----------
    P : numpy.ndarray
        Precision matrix.
    delta : list or numpy.ndarray
        Delta values used to generate the mould.

    Returns
    -------
    result : numpy.ndarray
        Mould matrix.
    """
    n = len(delta)
    ax = []
    elms = []
    for i in range(n):
        lin = np.linspace(-delta[i], delta[i], delta[i] * 2 + 1)
        elms.append(len(lin))
        ax.append(lin)
    grid = np.meshgrid(*ax, indexing='ij')
    feat = np.empty((n, np.product(elms)))
    for i in range(n):
        feat[i] = grid[i].ravel()
    mould = gaussian_function(np.zeros(n), P, feat, 1)
    mould = mould.reshape(*elms)
    return mould
