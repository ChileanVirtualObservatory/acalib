from astropy.nddata import support_nddata, NDData
from astropy.table import Table

from skimage.measure import regionprops

from astropy import log

@support_nddata
def measure_shape(data, labeled_images, min_freq = None, max_freq = None, wcs=None):
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