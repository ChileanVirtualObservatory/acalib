import glob
import scipy.ndimage as scnd
from astropy.nddata import support_nddata, NDData
from .axes import *

@support_nddata
def cut(data,wcs=None,mask=None,unit=None,lower=None,upper=None):
    mslab=slab(data,lower,upper)
    scube=data[mslab]
    newwcs=wcs.slice(mslab,numpy_order=True)
    return scube,wcs=newwcs,unit=unit


def scale(inputCont, majorAxisTemplate):
	scaledData = []

	for i in np.arange(len(inputCont.images)):
	    prop = pr.fits_props(inputCont.images[i].data)
	    scale = majorAxisTemplate/prop['major']
	    scaledData.append(scnd.zoom(prop['orig'],scale))
	return scaledData

@support_nddata
def rotate(data, templateAngle):
	rotatedData = []
	angles = []

	for i in np.arange(len(data)):
	    prop = pr.img_props(data[i])
	    angles.append(templateAngle - prop['angle'])
	    rotatedData.append(scnd.rotate(data[i], angles[-1], reshape = True))
	return rotatedData, angles	

def limits(img,angle):
    if angle > 0:
        cx, cy = np.nonzero(np.array(img.T))
    else:
        cx, cy = np.nonzero(np.array(img))
        
    upper = (cx[0], cy[0])
    lower = (cx[-1], cy[-1])
    
    return upper,lower
    
def cropAndAlign(data,angles):
	alignedData = []
	shapes = []

	for i in np.arange(len(data)):
		upper, lower = limits(data[i],angles[i])
		crop = data[i][upper[1]:lower[1], upper[1]:lower[1]]
		shapes.append(list(crop.shape))

		alignedData.append(crop)

	minShape = tuple(np.amin(shapes,axis = 0))
	

	for i in np.arange(len(alignedData)):
		dxl = (alignedData[i].shape[0] - minShape[0])/2
		dxr = (alignedData[i].shape[0] + minShape[0])/2
		dyu = (alignedData[i].shape[1] - minShape[1])/2
		dyd = (alignedData[i].shape[1] + minShape[1])/2

		alignedData[i] = alignedData[i][dxl:dxr, dyu:dyd]

	return alignedData

