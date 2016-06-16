import os
os.environ['ETS_TOOLKIT'] = 'qt4'
from pyface.qt import QtGui, QtCore
from traits.api import HasTraits, Instance, on_trait_change
#from traitsui.api import View, Item, HSplit, Group
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor
import sys
sys.path.append('../../')
from acalib.core import *
from acalib.io import graph as gp
from astropy import log
from astropy.wcs import wcs
from matplotlib import pyplot as plt
from mayavi import mlab
from qrangeslider import QRangeSlider
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from tvtk.api import tvtk

############### The mayavi visualization. ###############
#########################################################
class MayaviVisualization(HasTraits):
	scene = Instance(MlabSceneModel, ())
	mesh_grid = 0
	wireframe = 0
	@on_trait_change('scene.activated')
	def update_volume(self):
		# Clear the current visualization.
		mlab.clf()
		# Copy the cube to avoid modifying the original, and then eliminate lower values.
		updated_cube = cube.copy()
		updated_cube.data = np.nan_to_num(updated_cube.data)
		updated_cube.data[updated_cube.data < min_threshold_factor*rms] = 0
		updated_cube.data[updated_cube.data > max_threshold_factor*rms] = 0
		# Create the volume visualization and axes.
		xi, yi, zi = updated_cube.get_mesh()
		ranges=updated_cube.get_ranges()
		if(volume_type == 'volume'):
			volume = mlab.pipeline.scalar_field(xi*shape[1]/shape[0], yi, zi, updated_cube.data)
			mlab.pipeline.volume(volume)
		if(volume_type == 'contour'): 
			mlab.contour3d(xi*shape[1]/shape[0], yi, zi, updated_cube.data, transparent=True, contours=10, opacity=0.5)
		# Adjust the axes to make it a cube.
		extent = [0, shape[1], 0, shape[1], 0, shape[2]]
		ax=mlab.axes(xlabel="VEL [km/s] ",ylabel="DEC [deg]",zlabel="RA [deg]",ranges=ranges,nb_labels=5, extent=extent)
		ax.axes.label_format='%.3f'
	
	@on_trait_change('scene.activated')
	def create_wireframe(self):
		# Create a wireframe around the selected area for the stacked view.
		mesh_data = np.random.random((2, 2, 2))
		self.mesh_grid = tvtk.RectilinearGrid()
		self.mesh_grid.point_data.scalars = mesh_data.ravel()
		self.mesh_grid.point_data.scalars.name = 'scalars'
		self.mesh_grid.dimensions = mesh_data.shape
		self.mesh_grid.x_coordinates = np.array((lower_stack_limit[0]*shape[1]/shape[0], upper_stack_limit[0]*shape[1]/shape[0]), dtype=np.int32)
		self.mesh_grid.y_coordinates = np.array((lower_stack_limit[1], upper_stack_limit[1]), dtype=np.int32)
		self.mesh_grid.z_coordinates = np.array((lower_stack_limit[2], upper_stack_limit[2]), dtype=np.int32)
		self.wireframe = mlab.pipeline.surface(self.mesh_grid, opacity=0)
		mlab.pipeline.surface(mlab.pipeline.extract_edges(self.wireframe), color=(0, 0, 0))
		
	def update_wireframe(self):
		self.mesh_grid.x_coordinates = np.array((lower_stack_limit[0]*shape[1]/shape[0], upper_stack_limit[0]*shape[1]/shape[0]), dtype=np.int32)
		self.mesh_grid.y_coordinates = np.array((lower_stack_limit[1], upper_stack_limit[1]), dtype=np.int32)
		self.mesh_grid.z_coordinates = np.array((lower_stack_limit[2], upper_stack_limit[2]), dtype=np.int32)
		self.wireframe.render()

	view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=300, width=300, show_label=False),
                resizable=True)
             
############### The QWidget containing the visualization for mayavi. ############### 
####################################################################################
class MayaviQWidget(QtGui.QWidget):
	def __init__(self, parent=None):
		QtGui.QWidget.__init__(self, parent)
		layout = QtGui.QVBoxLayout(self)
		layout.setContentsMargins(0,0,0,0)
		layout.setSpacing(0)
		self.visualization = MayaviVisualization()
		self.ui = self.visualization.edit_traits(parent=self,kind='subpanel').control
		layout.addWidget(self.ui)
		self.ui.setParent(self)

############### The matplotlib visualization. ###############
#############################################################
class MatplotlibVisualization(FigureCanvas):
	def __init__(self, type, parent=None, width=5, height=4, dpi=100):
		# Create a new figure.
		self.fig = Figure(figsize=(width, height), dpi=dpi)
		self.axes = self.fig.add_subplot(111)
		self.axes.hold(False)
		if(type == 'spectre'):
			# Add the RA and DEC axes.
			ranges = cube.get_ranges(lower_stack_limit, upper_stack_limit)
			vel = np.linspace(ranges[0], ranges[1], upper_stack_limit[0]-lower_stack_limit[0])
			spectre = cube.stack(lower_stack_limit, upper_stack_limit, (1,2))
			self.axes.plot(vel, spectre)
			self.axes.set_xlabel("VEL [km/s]")
			plt.subplots_adjust(bottom=0.15)
			if(ranges[0] > ranges[1]): self.axes.invert_xaxis()
		if(type == 'stacked'):
			# Calculate a new cube adjusted to limits.
			stacked_cube = cube.stack(lower_stack_limit, upper_stack_limit)
			ranges = cube.get_ranges(lower_stack_limit, upper_stack_limit)
			extent = ranges[4:6] + ranges[2:4]
			img = self.axes.imshow(stacked_cube, origin='lower', extent=extent)
			# Reduce the number of ticks.
			xticks = np.linspace(ranges[4], ranges[5], 4)
			yticks = np.linspace(ranges[2], ranges[3], 4)
			self.axes.set_xticks(xticks)
			self.axes.set_yticks(yticks)
			# Add the colorbar.
			plt.colorbar(img, ax=self.axes)
		FigureCanvas.__init__(self, self.fig)
		self.setParent(parent)
		FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)
		
	def update_spectre(self):
		ranges = cube.get_ranges(lower_stack_limit, upper_stack_limit)
		vel = np.linspace(ranges[0], ranges[1], upper_stack_limit[0]-lower_stack_limit[0])
		spectre = cube.stack(lower_stack_limit, upper_stack_limit, (1,2))
		self.axes.plot(vel, spectre)
		self.axes.set_xlabel("VEL [km/s]")
		self.fig.canvas.draw()
		
	def update_stacked(self):
		self.fig.clf()
		self.axes = self.fig.add_subplot(111)
		stacked_cube = cube.stack(lower_stack_limit, upper_stack_limit)
		ranges = cube.get_ranges(lower_stack_limit, upper_stack_limit)
		extent = ranges[4:6] + ranges[2:4]
		img = self.axes.imshow(stacked_cube, origin='lower', extent=extent)
		xticks = np.linspace(ranges[4], ranges[5], 4)
		yticks = np.linspace(ranges[2], ranges[3], 4)
		self.axes.set_xticks(xticks)
		self.axes.set_yticks(yticks)
		plt.colorbar(img, ax=self.axes) 
		self.fig.canvas.draw()

############### The QWidget containing the visualization for matplotlib. ############### 
########################################################################################
class MatplotlibQWidget(QtGui.QWidget):
	def __init__(self, type, parent=None):
		QtGui.QWidget.__init__(self, parent)
		layout = QtGui.QVBoxLayout(self)
		layout.setContentsMargins(0,0,0,0)
		layout.setSpacing(0)
		self.main_widget = QtGui.QWidget(self)
		self.visualization = MatplotlibVisualization(type, self.main_widget, width=5, height=4, dpi=100)
		layout.addWidget(self.visualization)
	
############### Load the file and the data inside. ############### 
##################################################################
folder = '../../../../fits/'
fits_file = folder+'M100line.image.fits'
#fits_file = folder+'Boom.cm.cln.fits'
#fits_file = folder+'Orion.methanol.cbc.contsub.image.fits'
#fits_file = '~/calibrated.ms.contsub.bin4.line.fits'

# Load from container.
c = AContainer()
c.load(fits_file)
cube=c.adata[0]

# Metadata.
header = cube.meta

# This variable toggles between a volume and a contour surface in the mayavi visualization.
volume_type = 'volume'

# Simple Checks.
shape=cube.shape()
pixels=cube.count()
tflux=cube.flux()
vflux=cube.variance()
mmin = cube.min()
mmax = cube.max()

# Determine the threshold.
rms=cube.estimate_rms()
max_threshold = int(mmax[0]/rms)
min_threshold_factor = 5
max_threshold_factor = max_threshold

# Create the limits for the stacked view.
lower_stack_limit = (0,0,0)
upper_stack_limit = shape

# WCS Limits.
l1=cube.wcs_limits(0)
l2=cube.wcs_limits(1)
l3=cube.wcs_limits(2)
iw=cube.index_to_wcs((20,300,300))

############### Create the windows and add the container. ###############
#########################################################################
app = QtGui.QApplication.instance()
container = QtGui.QWidget()
container.setWindowTitle("Embedding Mayavi in a PyQt4 Application")
layout = QtGui.QGridLayout(container)

# Initialize the window and the main widgets.
label = QtGui.QLabel(container)
label.setText("Vista de volumen")
label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
layout.addWidget(label, 0, 3)
mayavi_widget = MayaviQWidget()
layout.addWidget(mayavi_widget, 1, 2, 1, 3)

label = QtGui.QLabel(container)
label.setText("Espectro de VEL")
label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
layout.addWidget(label, 3, 0)
matplotlib_spectre = MatplotlibQWidget('spectre')
layout.addWidget(matplotlib_spectre, 4, 0, 7, 5)

label = QtGui.QLabel(container)
label.setText("Vista acumulada")
label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
layout.addWidget(label, 0, 6)
matplotlib_stacked = MatplotlibQWidget('stacked')
layout.addWidget(matplotlib_stacked, 1, 5, 3, 3)

container.show()
window = QtGui.QMainWindow()
screen_resolution = app.desktop().screenGeometry()
window.resize(screen_resolution.width(), screen_resolution.height())
window.setCentralWidget(container)
window.show()

###############  GUI elements and related functions. ############### 
####################################################################

# Functions to update the view when the threshold is changed.
def threshold_slider_change(extreme_changed):
	global min_threshold_factor, max_threshold_factor	
	if(extreme_changed == 'start' and int(min_threshold_factor) != threshold_slider.start()):
		min_threshold_factor = threshold_slider.start()
		min_threshold_text.setText(str(min_threshold_factor))
		mayavi_widget.visualization.update_volume()
		mayavi_widget.visualization.create_wireframe()
	if(extreme_changed == 'end' and int(max_threshold_factor) != threshold_slider.end()):
		max_threshold_factor = threshold_slider.end()
		max_threshold_text.setText(str(threshold_slider.end()))
		mayavi_widget.visualization.update_volume()
		mayavi_widget.visualization.create_wireframe()
	
def threshold_text_change(extreme_changed):
	global min_threshold_factor, max_threshold_factor	
	if(extreme_changed == 'start'):
		if(float(min_threshold_text.text()) < 1):
			min_threshold_text.setText('1')
		if(float(min_threshold_text.text()) >= float(max_threshold_text.text())):
			min_threshold_text.setText(str(float(max_threshold_text.text())-0.01))
		min_threshold_factor = float(min_threshold_text.text())
		threshold_slider.setStart(int(float(min_threshold_text.text())))
	if(extreme_changed == 'end'):
		if(float(max_threshold_text.text()) > max_threshold):
			max_threshold_text.setText(str(max_threshold))
		if(float(max_threshold_text.text()) <= float(min_threshold_text.text())):
			max_threshold_text.setText(str(float(min_threshold_text.text())+0.01))
		max_threshold_factor = float(max_threshold_text.text())
		threshold_slider.setEnd(int(float(max_threshold_text.text())))
	mayavi_widget.visualization.update_volume()
	mayavi_widget.visualization.create_wireframe()
		
# Label for the threshold slider.
label = QtGui.QLabel(container)
label.setText("Threshold")
label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
layout.addWidget(label, 2, 3)
# Slider used to change the threshold range of the volume.
threshold_slider = QRangeSlider(window)
threshold_slider.setMin(1)
threshold_slider.setMax(max_threshold)
threshold_slider.setRange(min_threshold_factor, max_threshold)
threshold_slider.startValueChanged.connect(lambda: threshold_slider_change('start'))
threshold_slider.endValueChanged.connect(lambda: threshold_slider_change('end'))
layout.addWidget(threshold_slider, 3, 3)
# Textbox for display and change directly the start limit value.
min_threshold_text = QtGui.QLineEdit(window)
min_threshold_text.setFixedWidth(50)
min_threshold_text.setValidator(QtGui.QDoubleValidator(0., float(max_threshold), 2, min_threshold_text))
min_threshold_text.validator().setNotation(0)
min_threshold_text.returnPressed.connect(lambda: threshold_text_change('start'))
min_threshold_text.setText(str(min_threshold_factor))
layout.addWidget(min_threshold_text, 3, 2)
# Textbox for display and change directly the end limit value.
max_threshold_text = QtGui.QLineEdit(window)
max_threshold_text.setFixedWidth(50)
max_threshold_text.setValidator(QtGui.QDoubleValidator(0., 99., 2, max_threshold_text))
max_threshold_text.validator().setNotation(0)
max_threshold_text.returnPressed.connect(lambda: threshold_text_change('end'))
max_threshold_text.setText(str(max_threshold))
layout.addWidget(max_threshold_text, 3, 4)

# Functions to update the view when the stack range is changed.
user_control = True
def stack_slider_change(slider, textbox, extreme_changed):
	global lower_stack_limit, upper_stack_limit
	if(user_control and extreme_changed == 'start' and (vel_slider.start(), dec_slider.start(), ra_slider.start()) != lower_stack_limit):
		lower_stack_limit = (vel_slider.start(), dec_slider.start(), ra_slider.start())
		textbox.setText(str(slider.start()))
		matplotlib_stacked.visualization.update_stacked()
		matplotlib_spectre.visualization.update_spectre()
		mayavi_widget.visualization.update_wireframe()
	if(user_control and extreme_changed == 'end' and (vel_slider.end(), dec_slider.end(), ra_slider.end()) != upper_stack_limit):
		upper_stack_limit = (vel_slider.end(), dec_slider.end(), ra_slider.end())
		textbox.setText(str(slider.end()))
		matplotlib_stacked.visualization.update_stacked()
		matplotlib_spectre.visualization.update_spectre()
		mayavi_widget.visualization.update_wireframe()

def stack_text_change(slider, textbox, extreme_changed):
	global lower_stack_limit, upper_stack_limit
	if(extreme_changed == 'start'):
		if(int(textbox.text()) >= slider.end()):
			textbox.setText(str(slider.end()-1))
		lower_stack_limit = (int(min_vel_text.text()), int(min_dec_text.text()), int(min_ra_text.text()))
		slider.setStart(int(textbox.text()))
	if(extreme_changed == 'end'):
		if(int(textbox.text()) <= slider.start()):
			textbox.setText(str(slider.start()+1))
		upper_stack_limit = (int(max_vel_text.text()), int(max_dec_text.text()), int(max_ra_text.text()))
		slider.setEnd(int(textbox.text()))
	matplotlib_stacked.visualization.update_stacked()
	matplotlib_spectre.visualization.update_spectre()
	mayavi_widget.visualization.update_wireframe()
	
# Label for the ra range slider.
label = QtGui.QLabel(container)
label.setText("Rango de RA")
label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
layout.addWidget(label, 4, 6)
# Slider used to change the ra range of the stacked view.
ra_slider = QRangeSlider(window)
ra_slider.setMax(shape[2])
ra_slider.setRange(0, shape[2])
ra_slider.startValueChanged.connect(lambda: stack_slider_change(ra_slider, min_ra_text, 'start'))
ra_slider.endValueChanged.connect(lambda: stack_slider_change(ra_slider, max_ra_text, 'end'))
layout.addWidget(ra_slider, 5, 6)
# Textbox for display and change directly the start limit value.
min_ra_text = QtGui.QLineEdit(window)
min_ra_text.setFixedWidth(50)
min_ra_text.setValidator(QtGui.QIntValidator(0, shape[2], min_ra_text))
min_ra_text.returnPressed.connect(lambda: stack_text_change(ra_slider, min_ra_text, 'start'))
min_ra_text.setText('0')
layout.addWidget(min_ra_text, 5, 5)
# Textbox for display and change directly the end limit value.
max_ra_text = QtGui.QLineEdit(window)
max_ra_text.setFixedWidth(50)
max_ra_text.setValidator(QtGui.QIntValidator(0, shape[2], max_ra_text))
max_ra_text.returnPressed.connect(lambda: stack_text_change(ra_slider, max_ra_text, 'end'))
max_ra_text.setText(str(shape[2]))
layout.addWidget(max_ra_text, 5, 7)
	
# Label for the dec range slider.
label = QtGui.QLabel(container)
label.setText("Rango de DEC")
label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
layout.addWidget(label, 6, 6)
# Slider used to change the dec range of the stacked view.
dec_slider = QRangeSlider(window)
dec_slider.setMax(shape[1])
dec_slider.setRange(0, shape[1])
dec_slider.startValueChanged.connect(lambda: stack_slider_change(dec_slider, min_dec_text, 'start'))
dec_slider.endValueChanged.connect(lambda: stack_slider_change(dec_slider, max_dec_text, 'end'))
layout.addWidget(dec_slider, 7, 6)
# Textbox for display and change directly the start limit value.
min_dec_text = QtGui.QLineEdit(window)
min_dec_text.setFixedWidth(50)
min_dec_text.setValidator(QtGui.QIntValidator(0, shape[1], min_dec_text))
min_dec_text.returnPressed.connect(lambda: stack_text_change(dec_slider, min_dec_text, 'start'))
min_dec_text.setText('0')
layout.addWidget(min_dec_text, 7, 5)
# Textbox for display and change directly the end limit value.
max_dec_text = QtGui.QLineEdit(window)
max_dec_text.setFixedWidth(50)
max_dec_text.setValidator(QtGui.QIntValidator(0, shape[1], max_dec_text))
max_dec_text.returnPressed.connect(lambda: stack_text_change(dec_slider, max_dec_text, 'end'))
max_dec_text.setText(str(shape[1]))
layout.addWidget(max_dec_text, 7, 7)

# Label for the vel range slider.
label = QtGui.QLabel(container)
label.setText("Rango de VEL")
label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
layout.addWidget(label, 8, 6)
# Slider used to change the vel range of the stacked view.
vel_slider = QRangeSlider(window)
vel_slider.setMax(shape[0])
vel_slider.setRange(0, shape[0])
vel_slider.startValueChanged.connect(lambda: stack_slider_change(vel_slider, min_vel_text, 'start'))
vel_slider.endValueChanged.connect(lambda: stack_slider_change(vel_slider, max_vel_text, 'end'))
layout.addWidget(vel_slider, 9, 6)
# Textbox for display and change directly the start limit value.
min_vel_text = QtGui.QLineEdit(window)
min_vel_text.setFixedWidth(50)
min_vel_text.setValidator(QtGui.QIntValidator(0, shape[0], min_vel_text))
min_vel_text.returnPressed.connect(lambda: stack_text_change(vel_slider, min_vel_text, 'start'))
min_vel_text.setText('0')
layout.addWidget(min_vel_text, 9, 5)
# Textbox for display and change directly the end limit value.
max_vel_text = QtGui.QLineEdit(window)
max_vel_text.setFixedWidth(50)
max_vel_text.setValidator(QtGui.QIntValidator(0, shape[0], max_vel_text))
max_vel_text.returnPressed.connect(lambda: stack_text_change(vel_slider, max_vel_text, 'end'))
max_vel_text.setText(str(shape[0]))
layout.addWidget(max_vel_text, 9, 7)

# Label for the header data.
label = QtGui.QLabel(container)
label.setText("Metadatos")
label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
layout.addWidget(label, 0, 0)
# Add the header data under the view.
header_text = QtGui.QTextEdit(window)
header_text.setText(str(header))
header_text.setReadOnly(True)
layout.addWidget(header_text, 1, 0, 1, 2)

# Buttons to change between a volume and a contour.
def change_volume_type(type):
	global volume_type
	if(type != volume_type):
		volume_type = type
		mayavi_widget.visualization.update_volume()
		mayavi_widget.visualization.create_wireframe()

volume_button = QtGui.QPushButton("Volumen", window)
volume_button.clicked.connect(lambda: change_volume_type('volume'))
layout.addWidget(volume_button, 2, 1)
contour_button = QtGui.QPushButton("Contorno", window)
contour_button.clicked.connect(lambda: change_volume_type('contour'))
layout.addWidget(contour_button, 3, 1)

# Button to reset the wireframe to its original boundaries.
def reset_stack():
	global user_control, lower_stack_limit, upper_stack_limit
	user_control = False
	lower_stack_limit = (0,0,0)
	upper_stack_limit = shape
	ra_slider.setRange(0, shape[2])
	min_ra_text.setText('0')
	max_ra_text.setText(str(shape[2]))
	dec_slider.setRange(0, shape[1])
	min_dec_text.setText('0')
	max_dec_text.setText(str(shape[1]))
	vel_slider.setRange(0, shape[0])
	min_vel_text.setText('0')
	max_vel_text.setText(str(shape[0]))
	matplotlib_stacked.visualization.update_stacked()
	matplotlib_spectre.visualization.update_spectre()
	mayavi_widget.visualization.update_wireframe()
	user_control = True
	
reset_button = QtGui.QPushButton("Restaurar los valores originales de los rangos", window)
reset_button.clicked.connect(reset_stack)
layout.addWidget(reset_button, 10, 6)

# Start the main event loop.
app.exec_()


	
