from astropy import log
import numpy as np
import astropy.constants as const

import astropy.wcs as wcs
import astropy.units as u
from astropy.table.table import Table
from astropy.nddata import NDDataRef
import datetime
import copy

from acalib import *
from .convert import *


class Universe:
    """
    A synthetic universe where to put synthetic objects.
    """

    def __init__(self):
        self.sources = dict()

    def create_source(self, name, pos):
        """
        A source needs a name and a spatial position (alpha,delta).
        """
        self.sources[name] = Source(name, pos)

    def add_component(self, source_name, model):
        """
        To add a component a Component object must be instantiated (model), and added to
        a source called source_name.
        """
        self.sources[source_name].add_component(model)

    def _gen_sources_table(self):
        """
        Will generate a resumed table of each component.
        :return: an Astropy Table.
        """

        table = Table(names=("Source Name", "Comp ID", "Model", "Alpha", "Delta", "Redshift", "Radial Vel"),dtype=('S80','S80','S40','f8','f8','f8','f8'))
        for source in self.sources:
            for component in self.sources[source].comp:
                table.add_row((self.sources[source].name, component.comp_name, component.get_model_name(),
                          component.pos[0].value + component.offset[1].value, component.pos[1].value + component.offset[0].value, component.get_redshift().value,
                          component.get_velocity().value))

        return table

    def gen_cube(self, pos, ang_res, fov, freq, spe_res, bw, noise, cutlev):
        """
        Returns a container object where all the sources within the FOV and BW are projected in the
        primary object (NDDataRef), and with a sources astropy Table and all the parameters of the components
        in the successive tables as extensions

        This function needs the following parameters:
        - name    : name of the cube
        - pos     : right-ascension and declination center
        - ang_res : angular resolution x2
        - fov     : angular field of view x 2
        - freq    : spectral center (frequency)
        - spe_res : spectral resolution
        - bw      : spectral bandwidth
        - noise   : noise level
        - cutlev  : cut level
        """

        # Create a new WCS object.
        pos = to_deg(pos)
        ang_res = to_deg(ang_res)
        fov = to_deg(fov)
        #print pos,ang_res,fov
        freq = to_hz(freq)
        spe_res = to_hz(spe_res)
        bw = to_hz(bw)
        #print freq,spe_res,bw
        #w = wcs.WCS(naxis=3)
        #w.wcs.crval = np.array([pos[0].value, pos[1].value, freq.value])
        #w.wcs.restfrq = freq.value
        #w.wcs.cdelt = np.array([ang_res[0].value, ang_res[1].value, spe_res.value])
        mm = np.array([int(abs(fov[0]/ang_res[0])), int(abs(fov[1]/ang_res[1])), int(abs(bw/spe_res))])
        #w.wcs.crpix = mm / 2.0
        #w.wcs.ctype = ["RA---SIN", "DEC--SIN","FREQ"]
        data=np.zeros((mm[2],mm[1],mm[0]))
        meta=dict()
        meta['SIMPLE'] = True
        meta['BITPIX'] = -32
        meta['NAXIS'] = 3
        meta['NAXIS1'] = mm[0]
        meta['NAXIS2'] = mm[1]
        meta['NAXIS3'] = mm[2]
        #meta['NAXIS4'] = 1
        #meta['BSCALE'] = 1.0
        #meta['BZERO'] = 0.0
        meta['BMAJ'] = 3*ang_res[0].value
        meta['BMIN'] = 3*ang_res[0].value
        #meta['BPAN'] = 164.862686157
        meta['BTYPE'] = "Intensity"
        meta['OBJECT'] = "ACALIB Synthetic Object"
        meta['BUNITT'] = "JY/BEAM"
        meta['EQUINOX'] = 2000.0
        meta['RADESYS'] = "FK5"
        #meta['LONPOLE'] = 180.0
        #meta['LATPOLE'] = -5.37638888889
        #meta['PC001001'] = 1.0
        #meta['PC002001'] = 0.0
        #meta['PC003001'] = 0.0
        #meta['PC004001'] = 0.0
        #meta['PC001002'] = 0.0
        #meta['PC002002'] = 1.0
        #meta['PC003002'] = 0.0
        #meta['PC004002'] = 0.0
        #meta['PC001003'] = 0.0
        #meta['PC002003'] = 0.0
        #meta['PC003003'] = 1.0
        #meta['PC004003'] = 0.0
        #meta['PC001004'] = 0.0
        #meta['PC002004'] = 0.0
        #meta['PC003004'] = 0.0
       # meta['PC004004'] = 1.0

        meta['CTYPE1'] = "RA---SIN"
        meta['CRVAL1'] = pos[0].value
        meta['CDELT1'] = ang_res[0].value
        meta['CRPIX1'] = int(mm[0]/2)
        meta['CUNIT1'] = "deg"
        meta['CTYPE2'] = "DEC--SIN"
        meta['CRVAL2'] = pos[1].value
        meta['CDELT2'] = ang_res[1].value
        meta['CRPIX2'] = int(mm[1]/2)
        meta['CUNIT2'] = "deg"
        meta['CTYPE3'] = "FREQ"
        meta['CRVAL3'] = freq.value
        meta['CDELT3'] = spe_res.value
        meta['CRPIX3'] = int(mm[2]/2)
        meta['CUNIT3'] = "Hz"
        #meta['CTYPE4'] = "STOKES"
        #meta['CRVAL4'] = 1.0
        #meta['CDELT4'] = 1.0
        #meta['CRPIX4'] = 1.0
        #meta['CUNIT4'] = ""
        meta['PV2_1']= 0.0
        meta['PV2_2'] =0.0
        meta['RESTFRQ'] = freq.value
        meta['SPECSYS'] = "LSRK"
        #meta['ALTRVAL'] = -5429.74494828
        #meta['ALTRPIX'] = 1.0
        meta['VELREF ']= 257
        meta['COMMENT'] = "ACALIB Synthetic Object"
        meta['TELESCOP'] = "ACALB"
        meta['OBSERVER'] = "ACALIB"
        now=a=datetime.datetime.utcnow()
        meta['DATE-OBS'] = now.isoformat()
        meta['TIMESYS'] = "UTC"
        #meta['OBSRA'] = 83.8097916667
        #meta['OBSDEC'] = -5.37638888889
        #meta['OBSGEO-X'] = 2225142.18027
        #meta['OBSGEO-Y'] = -5440307.37035
        #meta['OBSGEO-Z'] = -2481029.85187
        meta['DATE'] = meta['DATE-OBS']
        meta['ORIGIN'] = "ACALIB"
        mywcs=wcs.WCS(meta)
        tab = [self._gen_sources_table()]
        cube = NDDataRef(data, wcs=mywcs,meta=meta,unit=u.Jy / u.beam)
        for source in self.sources:
            log.info('Projecting source ' + source)
            gen_tables = self.sources[source].project(cube, cutlev)
            #print gen_tables
            #gen_tables = self.sources[source].project(cube, noise)

            # add all tables generated from each component of the source
            # on the complete components dictionary.
            tab+=gen_tables
        cube.add_flux(2 * noise * (np.random.random(data.shape) - 0.5))
        cont=Container()
        cont.primary=cube

        cont.tables=tab

        return cont

    #def save_cube(self, cube, filename):
    #    """
    #    Wrapper function that saves a cube into a FITS (filename).
    #    """
    #    cube.save_fits(self.sources, filename)


class Source:
    """
    A generic source of electromagnetic waves with several components.
    """

    def __init__(self, name, pos):
        """
        :param name:    a name of the source
        """

        self.pos=to_deg(pos)
        self.name = name
        self.comp = list()

        log.info('Source \'' + name + '\' added\n')

    def add_component(self, model):
        """
        Defines a new component from a model.
        """

        code = self.name + '::' + str(len(self.comp) + 1)

        # create a deep copy of the model.
        model_cpy = copy.deepcopy(model)
        model_cpy.register(code, self.pos)
        self.comp.append(model_cpy)

        log.info('Added component ' + code + ' with model ' + model_cpy.info())

    def project(self, cube, limit):
        """
        Projects all components in the source to a cube.
        """

        component_tables = []
        log.info('Projecting Source at ' + str(self.pos))

        for component in self.comp:
            log.info('Projecting ' + component.comp_name)
            comp_table = component.project(cube, limit)

            if comp_table is not None:
                component_tables += [comp_table]
                meta_data = component.get_meta_data()

                if isinstance(meta_data, dict):
                    comp_table.meta = component.get_meta_data()
                else:
                    raise ValueError("get_meta_data return expected a dict, got %s instead" % type(meta_data))

        return component_tables

class Component:
    """Abstract component model"""

    def __init__(self):
        """
        Assume object in rest velocity/redshift
        """

        self.z = 0 * u.Unit("")

        self.comp_name = None
        self.pos = None

    def set_velocity(self, rvel):
        """Set radial velocity rvel. If rvel has no units, we assume km/s"""
        c = const.c.to('m/s')
        rvel = to_m_s(rvel)

        self.z = np.sqrt((1 + rvel/c) / (1 - rvel/c)) - 1

    def set_redshift(self, z):
        """Set the redshift"""
        self.z = z

    def get_velocity(self):
        """
        Get radial velocity rvel
        """
        z = self.z
        c = const.c.to('m/s')

        rv = c * (2 * z + np.square(z)) / (2 * z + np.square(z) + 2)

        return rv

    def get_redshift(self):
        """
        Get the redshift
        """

        return self.z

    def info(self):
        """
        Print relevant information of the component
        """

        return "(none)"

    def register(self, comp_name, pos):
        """
        Register the component name and angular position (alpha, delta)
        """

        self.comp_name = comp_name
        self.pos = pos

    def project(self, cube, limit):
        """
        Project the component in the cube and return the component astropy Table
        """
        pass

    def get_model_name(self):
        raise NotImplementedError("Extending class did not override the get_model_name method.")

    def get_meta_data(self):
        raise NotImplementedError("Get meta data not implemented.")
