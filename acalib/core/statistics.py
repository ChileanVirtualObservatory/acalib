from astropy.nddata import support_nddata, NDData
from astropy.table import Table
from skimage.measure import regionprops
from astropy import log
from .utils import *

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
    v=get_velocities(data,wcs,np.arange(data.shape[rdim]),restfrq) 
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
