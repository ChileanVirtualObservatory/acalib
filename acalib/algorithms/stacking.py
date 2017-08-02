import acalib
from .algorithm import Algorithm

import numpy
from numpy import mean
from scipy.stats import signaltonoise
import scipy.ndimage as scnd
from astropy.nddata import NDData,NDDataRef

import skimage

import matplotlib.pyplot as plt

class Stacking(Algorithm):
    """
    Create a stacked image using a template image and a set of different images from same object.
    """

    def default_params(self):
        pass

    def run(self, template_data, images):
        """
            Run the stacking algorithm given a template image and a container of images.

            Parameters
            ----------
            template_data : (M,N) numpy.ndarray
                Astronomical image.
            images : list of (M,N) numpy.ndarray
                A list of images.

            Returns
            -------
            result : (M,N) numpy.ndarray
                Image that adds all images after they are transformed so that their main object fits the position of the template_data.
        """
        # Transform input:
        if type(template_data) is NDData or type(template_data) is NDDataRef:
            template_data = template_data.data

        for i in range(len(images)):
            if type(images[i]) is not NDData or type(images[i]) is not NDDataRef:
                images[i] = NDDataRef(images[i])

        # Create a copy of template_data
        template_data = numpy.copy(template_data)
        # Get main object properties on it
        tprops = acalib.core.transform.fits_props(template_data)

        ORDER = 1 # Interpolation order for image operations

        for i in range(len(images)):
            data1 = images[i].data
            # Get main object properties
            props = acalib.core.transform.fits_props(data1)

            # Extend matrix so that the centroid is at the center of the image.
            max_img_dim = max(data1.shape[0]*2,data1.shape[1]*2)
            data2 = numpy.zeros((max_img_dim,max_img_dim))
            cy = data2.shape[0]//2-int(numpy.round(props['centroid'][0]))
            cx = data2.shape[1]//2-int(numpy.round(props['centroid'][1]))
            data2[cy:cy+data1.shape[0],cx:cx+data1.shape[1]] = data1

            # Rotate the image so that the larger radious is at 0 degrees:
            data3 = skimage.transform.rotate(data2,
                numpy.rad2deg(-props['angle']),cval=0,order=ORDER)

            # Rescale the image so that it meets the shape of the template:
            ysize = int(numpy.round(
                data3.shape[0]*tprops['ratio']/props['ratio']))
            data4 = skimage.transform.resize(data3,(ysize,data3.shape[1]),
                order=ORDER)
            # Extend matrix so that the centroid is at the center of the image.
            data5 = numpy.zeros(data2.shape)
            cy = (data3.shape[0]-data4.shape[0])//2
            data5 = acalib.core.transform.blit_add(data5,data4,(cy,0))

            # Rotate image again, now to match the angle of the template:
            data6 = skimage.transform.rotate(data5,
                numpy.rad2deg(tprops['angle']),cval=0,order=ORDER)

            # Scale the image:
            main_ratio = tprops['major']/props['major']
            data7 = skimage.transform.rescale(data6,main_ratio,
                order=ORDER,cval=0, clip=False)

            delta_y = int(numpy.round(tprops['centroid'][0]-data7.shape[0]//2))
            delta_x = int(numpy.round(tprops['centroid'][1]-data7.shape[1]//2))
            template_data = acalib.core.transform.blit_add(
                template_data,data7,(delta_y,delta_x))

        return template_data
