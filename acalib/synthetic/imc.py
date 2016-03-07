from numpy import random
from . import db
from .vu import Component
from ..core import flux as flx
import astropy.units as u
import numpy as np
from astropy import log
import urllib
import shutil
import os.path
from astropy.io import fits
from ..core import parameter as par
from ..core import atable
from ..core.atable import ATable

#INTEN_GROUP = [('default'), ('COv=0'), ('13COv=0'), ('HCO+, HC3N, CS, C18O, CH3OH, N2H, HDO')]
#INTEN_VALUES = [[0.1, 2], [20, 60], [5, 20], [1, 10]]
#DEFAULT_ISO_ABUND = {'13C': 1. / 30., '18O': 1. / 60., '17O': 1. / 120., '34S': 1. / 30., '33S': 1. / 120.,
#                         '13N': 1. / 30., 'D': 1./40.}
#DEFAULT_ABUND_RANGE=[10**-5,10**-6]
#GAUSS_STRINGS = ["Gaussian", "gaussian", "Gauss", "gauss", "normal", "Normal"]
#DEFAULT_CO_ABUND=1.0



DEFAULT_DBPATH= '../../bindata/db/ASYDO'


class IMC(Component):
    """ Interstellar Molecular Core """
    def __init__(self, mol_list, temp,dbpath=DEFAULT_DBPATH,equiv=u.doppler_radio):
        Component.__init__(self)
        self.equiv=equiv
        self.dbpath = dbpath
        self.temp = temp
        self.mol_list = mol_list

    def _draw(self,cube,flux,freq,cutoff):
        raise NotImplementedError("Draw not implemented.")

    def info(self):
        return ""

    def project(self, cube, cutoff):
        # TODO Make all this with astropy units from the call functions

        table = ATable(names=("Line Code", "Mol", "Ch Name", "Rest Freq", "Obs Freq", "Intensity"),dtype=('S80','S40','S40','f8','f8','f8'))

        dba = db.lineDB(self.dbpath)  # Maybe we can have an always open DB
        dba.connect()
        fwin = cube.wcs_limits(axis=2)
        # print "fwin",fwin
        cor_fwin = np.array(fwin/(1 + self.z))*u.Hz
        cor_fwin = cor_fwin.to(u.MHz).value
        # print "cor_fwin",cor_fwin
        counter = 0
        used = False
        for mol in self.mol_list:
            # For each molecule specified in the dictionary
            # load its spectral lines
            linlist = dba.getSpeciesLines(mol, cor_fwin[0], cor_fwin[1])  # Selected spectral lines for this molecule
            # rinte = INTEN_VALUES[0]
            # for j in range(len(INTEN_GROUP)):  # TODO: baaad python, try a more pythonic way..
            #     if mol in INTEN_GROUP[j]:
            #         rinte = INTEN_VALUES[j]
            abun = random.uniform(self.mol_list[mol][0], self.mol_list[mol][1])*u.Jy/u.beam

            for lin in linlist:
                counter += 1
                trans_temp = lin[5]*u.K
                inten = lin[4]

                if inten != 0:
                    inten = 10 ** inten

                flux = np.exp(-abs(trans_temp - self.temp) / trans_temp) * inten * abun
                flux = flux.value * u.Jy/u.beam
                # print trans_temp, self.temp, flux
                freq = (1 + self.z) * lin[3]*u.MHz  # TODO: astropy
                # print flux, cutoff
                if flux < cutoff: # TODO: astropy units!
                    log.info('    - Discarding ' + str(lin[1]) + ' at freq=' + str(freq) + '('+str(lin[3]*u.MHz)+') because I='+str(flux)+' < '+str(cutoff))
                    continue

                log.info('   - Projecting ' + str(lin[2]) + ' (' + str(lin[1]) + ') at freq=' + str(freq) + '('+str(lin[3]*u.MHz)+') intens='+ str(flux))
                self._draw(cube,flux,freq,cutoff)
                used = True

                # add line to the table.
                # TODO: modificar ultimo valor, que corresponde a la intensidad.
                table += (self.comp_name + "-l" + str(counter), mol, str(lin[2]), str(lin[3]),freq, flux)

        dba.disconnect()
        if not used:
            return None

        return table
        #return []

    def get_model_name(self):
        return "IMC"

    def get_meta_data(self):
        raise NotImplementedError("Get meta data not implemented.")

class GaussianIMC(IMC):
    def get_model_name(self):
        return "Gaussian IMC"

    def __init__(self, mol_list, temp, offset, std, angle, fwhm, gradient, dbpath=DEFAULT_DBPATH, equiv=u.doppler_radio):
        IMC.__init__(self,mol_list, temp, dbpath, equiv)
        self.offset = par.to_deg(offset)
        self.std = std
        self.angle = angle
        self.fwhm = fwhm
        self.gradient = gradient

    def _draw(self, cube, flux, freq, cutoff):
        new_pos = self.pos + self.offset
        mu, p = flx.clump_to_gauss(new_pos, self.std, self.angle, freq, self.fwhm, self.gradient)
        mcub, lower, upper = flx.create_gauss_flux(cube, mu, p, flux, cutoff)
        cube.add_flux(mcub, lower, upper)

    def info(self):
        return "species = " + str(self.mol_list) + " temp = "+str(self.temp)+" offset = "+str(self.offset)+" std = "+str(self.std)+" angle = "+str(self.angle)+" fwhm = "+str(self.fwhm)+"  gradient = "+str(self.gradient)

    def get_meta_data(self):
        return {
            # TODO: make the keys constants, NOT HARDCODED!!!
            'CRVAL1': self.pos[0],
            'CRVAL2': self.pos[1],
            '__STD1': self.std[0],
            '__STD2': self.std[1],
            '__GRD1': self.gradient[0],
            '__GRD2': self.gradient[1],
            '__FWHM': self.fwhm,
            '__TEMP': self.temp,
            '__ANGL': self.angle,
        }

################ ATTIC #######################

        #for mol in mol_list.split(','):
        #    abun = random.uniform(abun_range[1], abun_range[0])
        #    if mol in ('COv=0', '13COv=0', 'C18O', 'C17O', '13C18O'):
        #        abun += abun_CO
        #    for iso in iso_abun:
        #        if iso in mol:
        #            abun *= iso_abun[iso]
        #    self.intens[mol] = abun
#    def loader(self, uri, destination=""):
#        """
#        Gets the image/Fits to use as model
#        Supports Fits only.
#        """
#
#        fileName = uri.split("/")[-1]
#        self.fullDestination = os.path.join(destination, fileName)
#        s = urllib.urlretrieve(uri)
#
#        # Copy the file to the selected destination
#        shutil.copy2(s[0], self.fullDestination)
#        self.extension = os.path.splitext(self.fileName)[-1]
#        if (self.extension == "fits"):
#            return fits.open(self.fullDestination)[0].data
#        else:
#            raise ValueError("Wrong File Type (Only .fits Currently supported)")

 #   def change_intensities(self, intens):
 #       '''User defined dictionary in the form {molecule: intensity}'''
 #       self.intens = intens

        #   self._draw_func=self._draw_gauss
        #else:
           # Assuming an image template URI (fits format)
           # Download the URI
           # Maybe rotate it? (not sure)
           # Load the URI and put it in _image
        #   self._image=self.loader(URI)
        #   self._draw_func=self._draw_image
        #for mol in mol_list.split(','):
        #    abun = random.uniform(abun_range[1], abun_range[0])
        #    if mol in ('COv=0', '13COv=0', 'C18O', 'C17O', '13C18O'):
        #        abun += abun_CO
        #    for iso in iso_abun:
        #        if iso in mol:
        #            abun *= iso_abun[iso]
        #    self.intens[mol] = abun
