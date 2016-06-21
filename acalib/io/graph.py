from astropy import log
import numpy as np
from astropy.wcs import wcs
from mayavi import mlab
from ..core.utils import get_mesh,get_ranges

def volume(cube,wcs):
     figure = mlab.figure('Volume Plot')
     xi,yi,zi=get_mesh(cube)
     ranges=get_ranges(cube,wcs)
     grid = mlab.pipeline.scalar_field(xi, yi, zi, cube)
     mmin = cube.min()
     mmax = cube.max()
     mlab.pipeline.volume(grid)#,vmin=mmin, vmax=mmin)
     ax=mlab.axes(xlabel="VEL [km/s] ",ylabel="DEC [deg]",zlabel="RA [deg]",ranges=ranges,nb_labels=5)
     ax.axes.label_format='%.3f'
     mlab.colorbar(title='flux', orientation='vertical', nb_labels=5)
     mlab.show()


def contour(cube,wcs):
     figure = mlab.figure('Contour Plot')
     xi,yi,zi=get_mesh(cube)
     ranges=get_ranges(cube,wcs)
     mmin = cube.min()
     mmax = cube.max()
     mlab.contour3d(xi,yi,zi,cube,transparent=True,contours=10,opacity=0.5)
     ax=mlab.axes(xlabel="VEL [km/s] ",ylabel="DEC [deg]",zlabel="RA [deg]",ranges=ranges,nb_labels=5)
     ax.axes.label_format='%.3f'
     mlab.colorbar(title='flux', orientation='vertical', nb_labels=5)
     mlab.show()

#def stacked(cube):
#     figure = mlab.figure('Stacked Plot')
#     ranges=cube.get_ranges()
#     ranges=[ranges[2],ranges[3],ranges[4],ranges[5],ranges[0],ranges[1]]
#     img=cube.stack()
#     mlab.imshow(img)
#     ax=mlab.axes(xlabel="DEC [deg]",ylabel="RA [deg]",zlabel="VEL [km/s] ",ranges=ranges,nb_labels=5)
#     ax.axes.label_format='%.3f'
#     mlab.colorbar(title='flux', orientation='vertical', nb_labels=5)
#     mlab.show()


#def velocity(cube):
#     figure = mlab.figure('Velocity Map')
#     ranges=cube.get_ranges()
#     nn=cube.data.shape[0]
#     ranges=[ranges[2],ranges[3],ranges[4],ranges[5],-nn/2,nn/2]
#     #nn=cube.data.shape[0]
#     #vect=np.linspace(0.0,1.0,nn)
#     #vfield=np.average(cube.data,axis=0,weights=vect)
#     rms=cube.estimate_rms()
#     afield=np.argmax(cube.data,axis=0) - nn/2
#     vfield=np.max(cube.data,axis=0)
#     afield[vfield<1.5*rms]=0
#     afield[afield==-20]=0
#     mlab.surf(afield,warp_scale="auto")
#     ax=mlab.axes(xlabel="DEC [deg]",ylabel="RA [deg]",zlabel="PIX",ranges=ranges,nb_labels=5)
#     ax.axes.label_format='%.3f'
#     mlab.colorbar(title='Pix', orientation='vertical', nb_labels=5)
#     mlab.show()



#    def animate(self, inte, rep=True):
#               #TODO: this is not ported to the new wcs usage: maybe we must use wcsaxes to plot the wcs information...
#               """ Simple animation of the cube.
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


