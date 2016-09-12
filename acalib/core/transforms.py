
import props as pr
import numpy as np

import glob
import scipy.ndimage as scnd



def scale(inputDir, majorAxisTemplate):
	data = glob.glob(inputDir+'/*.fits')
	scaledData = []

	for i in np.arange(len(data)):
	    prop = pr.fits_props(data[i])
	    print prop
	    scale = majorAxisTemplate/prop['major']
	    scaledData.append(scnd.zoom(prop['orig'],scale))
	

	return scaledData



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

	print shapes
	minShape = tuple(np.amin(shapes,axis = 0))
	

	for i in np.arange(len(alignedData)):
		dxl = (alignedData[i].shape[0] - minShape[0])/2
		dxr = (alignedData[i].shape[0] + minShape[0])/2
		dyu = (alignedData[i].shape[1] - minShape[1])/2
		dyd = (alignedData[i].shape[1] + minShape[1])/2

		alignedData[i] = alignedData[i][dxl:dxr, dyu:dyd]

	return alignedData

