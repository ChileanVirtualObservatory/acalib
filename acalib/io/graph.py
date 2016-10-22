from astropy import log
import numpy as np
from astropy.wcs import wcs
#from mayavi import mlab
from astropy.nddata import support_nddata
from ..core.analysis import *
import matplotlib.pyplot as plt

#TODO: complete the nddata support (i.e. data, meta...)
#TODO: make animation possible again

@support_nddata
def visualize(data,wcs=None,unit=None,contour=False):
    if data.ndim == 1:
        return visualize_plot(data,wcs,unit)
    elif data.ndim == 2:
        return visualize_image(data,wcs,unit,contour)
    elif data.ndim == 3:
        if contour:
            visualize_contour3D(data,wcs,unit)
        else:
            visualize_volume(data,wcs,unit)
    else:
        log.error("Data dimensions must be between 1 and 3")

@support_nddata
def visualize_plot(data,wcs=None,unit=None):
     if wcs is None:
         plt.plot(data)
         plt.ylabel(unit)
     else:
         #TODO: Implement x vector, but check why the wcs cannot be onedimensional!
         plt.plot(data)
         plt.ylabel(unit)
         plt.xlabel(wcs.axis_type_names[0])
     plt.show()
         
@support_nddata
def visualize_image(data,wcs=None,unit=None,contour=False):
     if wcs is None:
         plt.imshow(data, origin='lower', cmap=plt.cm.gist_heat)
         cb=plt.colorbar()
         cb.ax.set_ylabel(unit)
     else:

         gax=plt.subplot(111,projection=wcs)
         plt.imshow(data, origin='lower', cmap=plt.cm.gist_heat)
         g0=gax.coords[0]
         g1=gax.coords[1]
         g0.set_axislabel(wcs.axis_type_names[0])
         g1.set_axislabel(wcs.axis_type_names[1])
         g0.grid(color='yellow', alpha=0.5, linestyle='solid')
         g1.grid(color='yellow', alpha=0.5, linestyle='solid')
         cb=plt.colorbar()
         cb.ax.set_ylabel(unit)
     if contour:
         rms=estimate_rms(data)
         dmax=data.max()
         crs=np.arange(1,dmax/rms)
         plt.contour(data,levels=rms*crs,alpha=0.5)
     plt.show()

# TODO: Remove hardocded stuff
@support_nddata
def visualize_volume(data,wcs=None,unit=None):
     if wcs is None:
        log.error("WCS is needed by this function")
     figure = mlab.figure('Volume Plot')
     mesh=get_mesh(data)
     xi,yi,zi=mesh
     ranges=axes_ranges(data,wcs)
     grid = mlab.pipeline.scalar_field(xi, yi, zi, data)
     mmin = data.min()
     mmax = data.max()
     mlab.pipeline.volume(grid)#,vmin=mmin, vmax=mmin)
     ax=mlab.axes(xlabel="VEL [km/s] ",ylabel="DEC [deg]",zlabel="RA [deg]",ranges=ranges,nb_labels=5)
     ax.axes.label_format='%.3f'
     mlab.colorbar(title='flux', orientation='vertical', nb_labels=5)
     mlab.show()

@support_nddata
def visualize_contour3D(data,wcs=None,unit=None):
     if wcs is None:
        log.error("WCS is needed by this function")
     figure = mlab.figure('Contour Plot')
     mesh=get_mesh(data)
     xi,yi,zi=mesh
     ranges=axes_ranges(data,wcs)
     mmin = data.min()
     mmax = data.max()
     mlab.contour3d(xi,yi,zi,data,transparent=True,contours=10,opacity=0.5)
     ax=mlab.axes(xlabel="VEL [km/s] ",ylabel="DEC [deg]",zlabel="RA [deg]",ranges=ranges,nb_labels=5)
     ax.axes.label_format='%.3f'
     mlab.colorbar(title='flux', orientation='vertical', nb_labels=5)
     mlab.show()




#    def animate(self, inte, rep=True):
#               #TODO: this is not ported to the new wcs usage: maybe we must use wcsaxes to plot the wcs information...
#               """ Simple animation of the data.
#                               - inte       : time interval between frames
#                               - rep[=True] : boolean to repeat the animation
#                       """
#               fig = plt.figure()
#               self.im = plt.imshow(self.data[0, :, :], cmap=plt.get_cmap('jet'), vmin=self.data.min(), vmax=self.data.max(), \
#                                                                                                extent=(
#                                                                                                                self.alpha_border[0], self.alpha_border[1], self.delta_border[0],
#                                                                                                                self.delta_border[1]))
#               ani = animation.FuncAnimation(fig, self._updatefig, frames=range(len(self.freq_axis)), interval=inte, blit=True,
#                                                                                                                                       repeat=rep)
#               plt.show()

