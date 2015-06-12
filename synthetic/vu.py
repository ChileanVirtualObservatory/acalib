from core.cube import *
from astropy import log
import numpy as np
import astropy.units as u
import astropy.constants as const

### Core ASYDO Classes ###

class Universe:
    """
    A synthetic universe where to put synthetic objects.
    """

    def __init__(self):
        self.sources = dict()

    def create_source(self, name, alpha, delta):
        """
        A source needs a name and a spatial position (alpha,delta).
        """
        self.sources[name] = Source(name, alpha, delta)

    def add_component(self, source_name, model):
        """
        To add a component a Component object must be instantiated (model), and added to 
        a source called source_name.
        """
        self.sources[source_name].add_component(model)

    def _gen_sources_table(self):
        #TODO: generate table
    
    def gen_cube(self, name, alpha, delta, freq, ang_res, ang_fov, spe_res, spe_bw):
        """
        Returns a Cube object where all the sources within the FOV and BW are projected, and 
        a dictionary with a sources astropy Table and all the parameters of the components
        in the succesive Tables

        This function needs the following parameters:
        - name    : name of the cube
        - alpha   : right-ascension center
        - delta   : declination center
        - freq    : spectral center (frequency)
        - ang_res : angular resolution
        - ang_fov : angular field of view
        - spe_res : spectral resolution
        - spe_bw  : spectral bandwidth
        """
        # TODO: add a proper constructor to the Cube class in core
        #cube = Cube(name, alpha, delta, freq, ang_res, ang_fov, spe_res, spe_bw)
        
        tables=dict()
        tables['sources']=self._gen_sources_table():              
 
        for src in self.sources:
            log.info('[Synthetic] Projecting source ' + src)
            dsource=self.sources[src].project(cube)
            tables.update(dsource)
        return cube,tables

    def save_cube(self, cube, filename):
        """
        Wrapper function that saves a cube into a FITS (filename).
        """
        cube.save_fits(self.sources, filename)

class Source:
    """
    A generic source of electromagnetic waves with several components.
    """

    def __init__(self, name, alpha, delta):
        """ Parameters:
               * name: a name of the source
               * alpha: right ascension 
               * delta: declination
        """
        log.info('[Synthetic] Source \'' + name + '\' added\n')
        self.alpha = alpha
        self.delta = delta
        self.name = name
        self.comp = list()

    def add_component(self, model):
        """ Defines a new component from a model.
        """
        code = self.name + '::' + str(len(self.comp) + 1)
        self.comp.append(model)
        model.register(code, self.alpha, self.delta)
        log.info('Added component ' + code + ' with model '+ model.info())

    def project(self, cube):
        """
        Projects all components in the source to a cube.
        """
        ctable=dict()
        for component in self.comp:
            log.info('Projecting ' + component.comp_name)
            table=component.project(cube)
            if table!=None:
               ctable[component.comp_name]=table
        return ctable

class Component:
    """Abstract component model"""

    def __init__(self):
        """
        Assume object in rest velocity/redshift
        """
        self.z = 0*u.Unit("")

    def set_velocity(self, rvel):
        """Set radial velocity rvel. If rvel has no units, we assume km/s"""
        c=const.c.to('km/s')
        if not isinstance(rvel,u.Quantity):
            rvel=rvel*u.km/u.s
        self.z = np.sqrt((1 + rvel/c) / (1 - rvel/c)) - 1

    def set_redshift(self, z):
        """Set the redshift"""
        self.z = z
 
    def get_velocity(self):
        """Get radial velocity rvel"""
        z=self.z
        c=const.c.to('km/s')
        rv = c*(2*z + np.square(z))/(2*z + np.square(z) + 2)
        return rv

    def set_redshift(self, z):
        """Get the redshift"""
        return(self.z)

    def info(self):
        """Print relevant information of the component"""
        return "(none)"

    def register(self, comp_name, alpha, delta):
        """Register the component name and angular position (alpha,delta)"""
        self.comp_name = comp_name
        self.alpha = alpha
        self.delta = delta

    def project(self, cube):
        """Project the component in the cube and return the component astropy Table"""
        pass

