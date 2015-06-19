from astropy.table import table
from core.cube import *
from astropy import log
import numpy as np
import astropy.units as u
import astropy.constants as const


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
        """
        Will generate a table with the following columns:
            - Component ID
            - Source name
            - Model
            - Alpha
            - Delta
            - Red shift (z)
            - Radial velocity

        :return: an astropy table.
        """

        # create columns for all fields.
        col_comp_id = []
        col_source_name = []
        col_model = []
        col_alpha = []
        col_delta = []
        col_redshift = []
        col_radial_velocity = []

        # run through all sources and components, and add those values to the above lists.
        for source in self.sources:
            for component in source.comp:
                col_comp_id.append(component.comp_name)
                col_source_name.append(source.name)
                col_model.append("Not yet.")
                col_alpha.append(component.alpha)
                col_delta.append(component.delta)
                col_redshift.append(component.get_redshift())
                col_radial_velocity.append(component.get_velocity())

        # create two lists for the table, this is to not pollute the table construction line.
        col_values = [col_comp_id, col_source_name, col_model, col_alpha, col_delta, col_redshift, col_radial_velocity]
        col_names = ["Comp ID", "Source name", "Model", "Alpha", "Delta", "Redshift", "Radial Vel"]

        return table.Table(col_values, col_names)
    
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

        tables = dict()
        tables['sources'] = self._gen_sources_table()

        for source in self.sources:
            log.info('[Synthetic] Projecting source ' + source)
            dsource = source.project(cube)
            tables.update(dsource)

        return cube, tables

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
        """
        :param name:    a name of the source
        :param alpha:   right ascension
        :param delta:   declination
        """

        self.alpha = alpha
        self.delta = delta
        self.name = name
        self.comp = list()

        log.info('[Synthetic] Source \'' + name + '\' added\n')

    def add_component(self, model):
        """
        Defines a new component from a model.
        """

        code = self.name + '::' + str(len(self.comp) + 1)

        # create a deep copy of the model.
        model_cpy = copy.deepcopy(model)
        model_cpy.register(code, self.alpha, self.delta)
        self.comp.append(model_cpy)

        log.info('Added component ' + code + ' with model ' + model_cpy.info())

    def project(self, cube):
        """
        Projects all components in the source to a cube.
        """

        component_tables = dict()

        for component in self.comp:
            log.info('Projecting ' + component.comp_name)
            table = component.project(cube)

            if table is not None:
                component_tables[component.comp_name] = table

        return component_tables

class Component:
    """Abstract component model"""

    def __init__(self):
        """
        Assume object in rest velocity/redshift
        """

        self.z = 0 * u.Unit("")

    def set_velocity(self, rvel):
        """Set radial velocity rvel. If rvel has no units, we assume km/s"""
        c = const.c.to('km/s')

        if not isinstance(rvel, u.Quantity):
            rvel = rvel*u.km/u.s

        self.z = np.sqrt((1 + rvel/c) / (1 - rvel/c)) - 1

    def set_redshift(self, z):
        """Set the redshift"""
        self.z = z
 
    def get_velocity(self):
        """
        Get radial velocity rvel
        """
        z = self.z
        c = const.c.to('km/s')

        rv = c * (2 * z + np.square(z)) / (2 * z + np.square(z) + 2)

        return rv

    def get_redshift(self, z):
        """
        Get the redshift
        """

        return self.z

    def info(self):
        """
        Print relevant information of the component
        """

        return "(none)"

    def register(self, comp_name, alpha, delta):
        """
        Register the component name and angular position (alpha, delta)
        """

        self.comp_name = comp_name
        self.alpha = alpha
        self.delta = delta

    def project(self, cube):
        """
        Project the component in the cube and return the component astropy Table
        """
        pass
