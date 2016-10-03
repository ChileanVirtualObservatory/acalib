import numpy as np

from astropy.nddata import support_nddata, NDData
from astropy.table import Table
from astropy import units as u

from skimage.filter import threshold_adaptive
from skimage.morphology import binary_opening
from skimage.morphology import disk

from skimage.measure import label
from skimage.measure import regionprops

from skimage.segmentation import clear_border
from astropy import log

@support_nddata
def gaussian_mix(data,prob = 0.05, precision=0.02, wcs=None):
    """
    Using a mixture of gaussians make an multiscale segmentation to get the region of interest of a 2D astronomical image.
    
    :param image: Velocity collapsed image
    :returns: list of skimage.measure.regionprops Objects, with detected regions properties
    """
    if len(data.shape) > 2:
        log.error("Only 2D images supported")        
        raise ValueError("Only 2D images supported")


    image_list = []
 
    
    image = data
    image[np.isnan(image)] = 0

    prob = prob
    dims = image.shape
    rows = dims[0]
    cols = dims[1]
    size = np.min([rows,cols])
    precision = size * precision 
    
    image = image.astype('float64')

    w_max = _optimal_w(image,prob)
    diff = (image - np.min(image)) / (np.max(image)- np.min(image))

    tt=w_max*w_max
    if tt%2==0:
       tt+=1
    g = threshold_adaptive(diff, tt,method='mean',offset=0)

    r = w_max/2
    rMin = 2*np.round(precision)


    while (r > rMin):
        background = np.zeros((rows,cols))
        selem = disk(r)
        sub = binary_opening(g,selem)
        sub = clear_border(sub)
        sub = label(sub)
        fts = regionprops(sub)
        
        image_list.append(NDData(sub,wcs=wcs))

        if len(fts) > 0:
            for props in fts:
                C_x, C_y = props.centroid

                radius = props.equivalent_diameter / 2.
                kern = 0.01 * np.ones((2*radius, 2*radius))
                krn = _kernelsmooth( x = np.ones((2*radius, 2*radius)), kern = kern)
                krn = np.exp(np.exp(krn))
                if np.max(krn) > 0:
                    krn = (krn-np.min(krn))/(np.max(krn)-np.min(krn))
                    background = _kernel_shift(background,krn, C_x, C_y)
        if np.max(background) > 0:
            background = (background-np.min(background))/(np.max(background)-np.min(background))
            diff = diff - background
        diff = (diff-np.min(diff))/(np.max(diff)-np.min(diff))
        tt=r*r
        if tt%2==0:
           tt+=1
        g = threshold_adaptive(diff,tt,method='mean',offset=0)
        r = np.round(r/2.)
    
    return image_list



def _optimal_w(image, p = 0.05):
    #radiusMin, radius Max and inc in percentages of the image size, p as [0,1] value, image is the original version
    radiusMin = 5
    radiusMax = 40
    inc = 1

    f = (image-np.min(image))/(np.max(image)-np.min(image))
    dims = f.shape
    rows = dims[0]
    cols = dims[1]

    maxsize = np.max([rows,cols])
    imagesize = cols*rows
    radius_thresh = np.round(np.min([rows,cols])/4.)
    unit = np.round(maxsize/100.)

    radiusMin = radiusMin*unit
    radiusMax = radiusMax*unit
    radiusMax = int(np.min([radiusMax,radius_thresh]))
    radius = radiusMin
    inc = inc*unit

    bg = np.percentile(f , p * 100)        
    fg = np.percentile(f, (1-p) * 100)
    min_ov = imagesize

    while(radius <= radiusMax):
        tt=radius*radius
        if tt%2==0:
           tt+=1
        
        g = threshold_adaptive(f,tt,method='mean',offset=0)
        ov = _bg_fg(f,g,bg,fg)
        if(ov < min_ov):
            w = radius
            min_ov = ov
        
        radius += inc
    return w

def _bg_fg(f,g,bg,fg):
    dims = f.shape
    rows = dims[0]
    cols = dims[1]
    fp = 0
    fn = 0
    for rowID in xrange(rows):
        for colId in xrange(cols):
            if g[rowID][colId] == True:
                if (np.abs(f[rowID][colId]- bg) < np.abs(f[rowID][colId]-fg)):
                    fp += 1
            elif g[rowID][colId] == False:
                if (np.abs(f[rowID][colId]- bg) > np.abs(f[rowID][colId]-fg)):
                    fn += 1
    overall = fp + fn
    return overall

def _kernelsmooth(x,kern, norm = True):
    # how many rows/cols of zeroes are used to pad.
    width = kern.shape[0]
    pad = int(width/2.)

    #record the width and height the input data matrix
    x_w = x.shape[0]
    x_h = x.shape[1]

    if norm:
        k = kern / np.sum(abs(kern))
    else:
        k = kern

    #Padding with zeros
    x_pad = np.lib.pad(x, ((pad,pad),(pad,pad)), 'constant')

    # Pre-allocate the final (smoothed) data matrix
    s = np.zeros( (x_h, x_w) )

    # Pre-allocate a temporary matrix for the iterative calculations
    temp = np.zeros((width,width))

    # Loop through the data to apply the kernel.
    for col in xrange(x_w):
        for row in xrange(x_h):
            temp = x_pad[row:(row+width),col:(col+width)]
            s[row][col] = np.sum(k*temp)

    return s

def _kernel_shift(back,kernel, x,y):
    rows_back = back.shape[0]
    cols_back = back.shape[1]
    rowsKernel = kernel.shape[0]
    colsKernel = kernel.shape[1]
    rowInit = int(x-rowsKernel/2)
    colInit = int(y-colsKernel/2)

    for row in xrange(rowsKernel):
        for col in xrange(colsKernel):
            if (rowInit + row < rows_back-1) and (colInit+col < cols_back-1):
                back[rowInit+row][colInit+col] = kernel[row][col]

    return back

@support_nddata
def measure_shape(data, labeled_images, wcs=None):
    objects = list()
    intensity_image = data
    for image in labeled_images:
        objs_properties = _get_shape(image, intensity_image)
        objects.extend(objs_properties)


    if len(objects) == 0:
        return Table()

    names  = ["CentroidRa", "CentroidDec", "MajorAxisLength","MinorAxisLength",
              "Area", "Eccentricity", "Solidity", "FilledPercentaje"]

    t = Table(rows = objects, names=names, meta={"name": "Object Shapes"})
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
                               eccentricity , solidity, filled))        

    return objs_properties