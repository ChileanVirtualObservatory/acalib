from astropy.nddata import support_nddata, NDData
from astropy.table import Table
from skimage.measure import regionprops
from astropy import log
import numpy as np
from .utils import *
from acalib.core._morph import *


@support_nddata
def rms(data,mask=None):
    """Compute the RMS of data. If mask != None, then 
       we use that mask.
    """
    #TODO: check photutils background estimation for using that if possible
    if mask is not None:
        data=fix_mask(data,mask)
    mm=data*data
    rms=np.sqrt(mm.sum()*1.0/mm.size)
    return rms

@support_nddata
def snr_estimation(data,mask=None,noise=None,points=1000,full_output=False):
    """Heurustic that uses the inflexion point of the thresholded RMS to estimate 
       where signal is dominant w.r.t. noise
    """ 
    if noise is None:
       noise=rms(data,mask)
    x=[]
    y=[]
    n=[]
    sdata=data[data>noise]
    for i in range(1,int(points)):
        val=1.0 + 2.0*i/points
        sdata=sdata[sdata>val*noise]
        if sdata.size < 2:
            break
        n.append(sdata.size)
        yval=sdata.mean()/noise
        x.append(val)
        y.append(yval)
    y=np.array(y)
    v=y[1:]-y[0:-1]
    p=v.argmax() + 1
    snrlimit=x[p]
    if full_output==True:
       return snrlimit,noise,x,y,v,n,p 
    return snrlimit


def _moment(data,order,wcs=None,mask=None,unit=None,restfrq=None): 
    if wcs is None: 
        log.error("A world coordinate system (WCS) is needed") 
        return None 
    data=fix_mask(data,mask) 
    dim=wcs.wcs.spec 
    rdim=data.ndim - 1 - dim 
    v=spectral_velocities(data,wcs,np.arange(data.shape[rdim]),restfrq) 
    v=v.value 
    #delta=np.mean(np.abs(v[:v.size-1] - v[1:v.size])) 
    #newdata=data.sum(axis=rdim)*delta 
    m0=data.sum(axis=rdim) 
    if order==0: 
        mywcs=wcs.dropaxis(dim) 
        return NDData(m0.data, uncertainty=None, mask=m0.mask,wcs=mywcs, meta=None, unit=unit) 
    #mu,alpha=np.average(data,axis=rdim,weights=v,returned=True) 
    mu,alpha=np.ma.average(data,axis=rdim,weights=v,returned=True) 
    m1=alpha*mu/m0 
    if order==1: 
        mywcs=wcs.dropaxis(dim) 
        return NDData(m1.data, uncertainty=None, mask=m1.mask,wcs=mywcs, meta=None, unit=u.km/u.s) 
    v2=v*v 
    var,beta=np.ma.average(data,axis=rdim,weights=v2,returned=True) 
    #var,beta=data.average(axis=rdim,weights=v2,returned=True) 
    m2=np.sqrt(beta*var/m0 - m1*m1) 
    if order==2: 
        mywcs=wcs.dropaxis(dim) 
        return NDData(m2.data, uncertainty=None, mask=m2.mask,wcs=mywcs, meta=None, unit=u.km*u.km/u.s/u.s) 
    log.error("Order not supported") 
    return None 

@support_nddata
def integrate(data, wcs=None, mask=None, unit=None, axis=(0)):
    """ Returns a numpy array (no WCS!) with the integration results. For updated WCS use spectra, moments, or similar. """
    if mask is not None:
        data=fix_mask(data,mask)
    newdata = np.sum(data, axis=axis)
    mask = np.isnan(newdata)
    return newdata

        
# Should return a NDData 
@support_nddata 
def moment0(data,wcs=None,mask=None,unit=None,restfrq=None): 
    return _moment(data,0,wcs,mask,unit,restfrq) 
 
@support_nddata 
def moment1(data,wcs=None,mask=None,unit=None,restfrq=None): 
    return _moment(data,1,wcs,mask,unit,restfrq) 
 
@support_nddata 
def moment2(data,wcs=None,mask=None,unit=None,restfrq=None): 
    return _moment(data,2,wcs,mask,unit,restfrq) 


@support_nddata
def measure_shape(data, labeled_images, min_freq = None, max_freq = None, wcs=None):
    """ Measure a few statistics from labeled images """
    #TODO: Document this function
    objects = list()
    intensity_image = data
    for image in labeled_images:
        objs_properties = _get_shape(image, intensity_image)
        objects.extend(objs_properties)


    if len(objects) == 0:
        return Table()

    names  = ["CentroidRa", "CentroidDec", "MajorAxisLength","MinorAxisLength",
              "Area", "Eccentricity", "Solidity", "FilledPercentaje", "MaxIntensity", "MinIntensity", "AverageIntensity" ]

    meta = {"name": "Object Shapes"}

    if min_freq is not None:
        meta["min_freq_hz"] = min_freq

    if max_freq is not None:
        meta["max_freq_hz"] = max_freq
            


    t = Table(rows = objects, names=names, meta=meta)
    return t

@support_nddata
def _get_shape(data, intensity_image, wcs=None):
    objs_properties = []
    fts = regionprops(data, intensity_image = intensity_image)
    
    for obj in fts:
        if wcs:
            matrix = wcs.pixel_scale_matrix
            deg_per_pix_x = matrix[0,0]
            deg_per_pix_y = matrix[1,1]

        centroid = wcs.celestial.all_pix2world(obj.centroid[0],obj.centroid[1],1) if wcs else obj.centroid
        centroid_ra = centroid[0] if wcs else centroid[0]
        centroid_dec = centroid[1] if wcs else centroid[1]
        major_axis = abs(deg_per_pix_x) * obj.major_axis_length if wcs else obj.major_axis_length 
        minor_axis = abs(deg_per_pix_x) * obj.minor_axis_length if wcs else obj.minor_axis_length
        area = obj.area * abs(deg_per_pix_x)  if wcs else obj.area
        eccentricity = obj.eccentricity
        solidity = obj.solidity
        filled = obj.area / obj.filled_area

        objs_properties.append((centroid_ra, centroid_dec, major_axis, minor_axis, area,
                               eccentricity , solidity, filled, obj.max_intensity, obj.min_intensity, obj.mean_intensity))        

    return objs_properties



# TODO: generalize this function... is not very generic :S
# TODO: update the WCS!
@support_nddata
def spectra(data,wcs=None,mask=None,unit=None,position=None,aperture=None):
    if position is None:
        # Get celestial center
        position=wcs.celestial.wcs.crval*u.deg
    if aperture is None:
        # Get 1 pixel aperture
        aperture=np.abs(wcs.celestial.wcs.cdelt[0])*u.deg
    if position.unit == u.pix and aperture.unit == u.pix:
        # TODO:  Here is the nasty part
        lb=np.array([0,            position[1].value - aperture.value, position[0].value - aperture.value])
        ub=np.array([data.shape[2],position[1].value + aperture.value, position[0].value + aperture.value])
    else:
        log.error("Not Implemented Yet!")
    specview=data[slab(data,lb,ub)]
    return specview.sum(axis=(1,2))


@support_nddata
def spectra_sketch(data,samples, random_state = None):
    """
    Create the sketch spectra using pixel samples.
    
    :param samples: Number of pixel samples used for the sketch.
    :type samples: int
    :returns: ( spectra (array), slices  (list)).
    """

    if random_state is not None:
        np.random.seed(random_state)

    dims = data.shape
    P_x = dims[2]
    P_x_range = range(P_x)
    P_y = dims[1]
    P_y_range = range(P_y)
    frec = dims[0]

    spectra = np.zeros(frec)

    for i in xrange(samples):
        x_ = np.random.choice(P_x_range,1)
        y_ = np.random.choice(P_y_range,1)
        pixel = data[:, y_, x_] 
        pixel_masked = _pixel_processing(pixel)
        spectra += pixel_masked
    spectra = _pixel_processing(spectra)

    slices = []
    min_slice = -1
    max_slice = -1
    for i in range(frec-1):
        if spectra[i] != 0:
            if min_slice == -1:
                min_slice = i
            else:
                if spectra[i+1] == 0:
                    max_slice = i+1
                    slices.append(slice(min_slice,max_slice))
                    min_slice = -1
                else:
                    if i == frec-2:
                        max_slice = i+1
                        slices.append(slice(min_slice,max_slice))

    return spectra,slices

def _pixel_processing(pixels):
    pixels = pixels.astype(np.float64)
    acum = _accumulating(pixels)
    diff = _differenting(acum)
    boxing = _segmenting(diff)
    boxing = _erosing(boxing)
    return _masking(boxing,pixels)


def _accumulating(pixels):
    return np.cumsum(pixels)

def _differenting(cumPixels):
    d = diff(cumPixels)
    return d

def _segmenting(diff):
    
    boxing = seg(diff)
    
    return boxing

def _erosing(boxing):
    
    boxing = eros(boxing)

    return boxing

def _masking(boxing, pixels):
    return boxing.reshape(-1)*pixels.reshape(-1)



#@support_nddata
#def vel_stacking(data,data_slice, wcs=None, mask=None,uncertainty=None, meta=None, unit=None):
#    """
#        Create an image collapsing the frecuency axis
#        
#        :param data_slice: Sector to be collapsed
#        :type data_slice: slice 
#        :returns: image (NDData): 2D-Array with the stacked cube.

#    """
#    if len(data.shape) != 3:
#        log.error("Cube needs to be a 3D array")
#        raise ValueError("Cube needs to be a 3D array")
#
#    dims = data.shape
#    subcube = data[data_slice, :,:]
#    stacked = np.sum(subcube,axis=0)
#    wcs = wcs.dropaxis(2)

#    return NDData(stacked, uncertainty=uncertainty, mask=mask,wcs=wcs, meta=meta, unit=unit)
