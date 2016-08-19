
import props as pr
import numpy as np

import glob
import scipy.ndimage as scnd



def scale(inputDir, majorAxisTemplate):

	inputDir = 'test'
	data = glob.glob(inputDir+'/*.fits')
	scaledData = []

	for i in np.arange(len(data)):
	    prop = pr.fits_props(data[i])
	    scale = majorAxisTemplate/prop['major']
	    scaledData.append(scnd.zoom(prop['orig'],scale))
	    
	return scaledData



def rotate(data, templateAngle):
	rotatedData = []

	for i in np.arange(len(data)):
	    prop = pr.img_props(data[i])
	    angle  = templateAngle - prop['angle']
	    rotatedData.append(scnd.rotate(arrays[i], angle, reshape = True))
	
	return rotatedData	



        
def limits(img,angle):
    if angle > 0:
        cx, cy = np.nonzero(np.array(img.T))
    else:
        cx, cy = np.nonzero(np.array(img))
        
    upper = (cx[0], cy[0])
    lower = (cx[-1], cy[-1])
    
    return upper,lower
    
    
def cropAndAlign(data,angle):
	alignedData = []
	shapes = []

	for i in np.arange(len(data)):
		upper, lower = limits(data[i],angle)
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