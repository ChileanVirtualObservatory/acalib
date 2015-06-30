from numpy import random
from synthetic import db
from synthetic.vu import Component
import core.flux as flx
import astropy.units as u

import urllib
import shutil
import os.path
from astropy.io import fits

INTEN_GROUP = [('default'), ('COv=0'), ('13COv=0'), ('HCO+, HC3N, CS, C18O, CH3OH, N2H')]
INTEN_VALUES = [[0.1, 2], [20, 60], [5, 20], [1, 10]]
DEFAULT_ISO_ABUND = {'13C': 1.0 / 30, '18O': 1.0 / 60, '17O': 1.0 / 120, '34S': 1.0 / 30, '33S': 1.0 / 120,
                         '13N': 1.0 / 30, 'D': 1.0 / 30}
DEFAULT_ABUND_RANGE=[10**-5,10**-6]
GAUSS_STRINGS = ["Gaussian", "gaussian", "Gauss", "gauss", "normal", "Normal"]
DEFAULT_CO_ABUND=1.0


DEFAULT_DBPATH= 'ASYDO'

class IMC(Component):
    """ Interstellar Molecular Core """

    def __init__(self, template,std, angle, mol_list, temp, fwhm, gradient, 
                 abun_range=DEFAULT_ABUND_RANGE, abun_CO=DEFAULT_CO_ABUND, iso_abun=DEFAULT_ISO_ABUND,dbpath=DEFAULT_DBPATH,equiv=u.doppler_radio, URI=""):
        Component.__init__(self)
        self.equiv=equiv
        self.dbpath = dbpath
        self.temp = temp
        self.std= std
        self.angle = angle
        self.fwhm = fwhm
        self.gradient = gradient
        self.intens = dict()
        if template in GAUSS_STRINGS:
           self._draw_func=self._draw_gauss
        else:
           # Assuming an image template URI (fits format)
           # Download the URI
           # Maybe rotate it? (not sure)
           # Load the URI and put it in _image
           self._image=self.loader(URI)
           self._draw_func=self._draw_image
        for mol in mol_list.split(','):
            abun = random.uniform(abun_range[1], abun_range[0])
            if mol in ('COv=0', '13COv=0', 'C18O', 'C17O', '13C18O'):
                abun += abun_CO
            for iso in iso_abun:
                if iso in mol:
                    abun *= iso_abun[iso]
            self.intens[mol] = abun


    def loader(self, uri, destination=""):
        """
        Gets the image/Fits to use as model
        Supports Fits only.
        """

        fileName = uri.split("/")[-1]
        self.fullDestination = os.path.join(destination, fileName)
        s = urllib.urlretrieve(uri)

        # Copy the file to the selected destination
        shutil.copy2(s[0], self.fullDestination)
        self.extension = os.path.splitext(self.fileName)[-1]
        if (self.extension == "fits"):
            return fits.open(self.fullDestination)[0].data
        else:
            raise ValueError("Wrong File Type (Only .fits Currently supported)")





    def change_intensities(self, intens):
        '''User defined dictionary in the form {molecule: intensity}'''
        self.intens = intens

    def _draw_gauss(self,cube,flux,freq,cutoff):
       pos=np.array([self.alpha,self.delta])
       (mu,P)=flx.clump_to_gauss(pos,self.std,self.angle,freq,self.fwhm,self.gradient,equiv)
       (mcub,lower,upper)=flx.create_gauss_flux(cube,mu,P,flux,cutoff)
       cube.add_flux(mcub,lower,upper)
       return retval

    def _draw_image(self,cube,flux,freq,cutoff):
       # TODO      
       pass

    def info(self):
       # TODO: implement
       return "IMC"
       #return "mol_list = " + str(self.intens.keys()) + " @ spa_form=" + str(self.spa_form) + ", spe_form=" + str(
       #     self.spe_form) + ", z=" + str(self.z) + ", grad=" + str(self.z_grad)

    def project(self, cube, cutoff):
        #arr_code = []
        #arr_mol = []
        #arr_chname = []
        #arr_rest_freq = []
        #arr_rad_vel = []
        #arr_fwhm = []
        #arr_temp = []
        dba = db.lineDB(self.dbpath) # Maybe we can have an always open DB
        dba.connect()
        fwin= cube.get_wcs_limits(axis=2)
        cor_fwin =  fwin / (1 + self.z)
        counter = 0
        used = False
        for mol in self.intens:
            # For each molecule specified in the dictionary
            # load its spectral lines

            linlist = dba.getSpeciesLines(mol, cor_fwin[0],
                                          cor_fwin[1])  # Selected spectral lines for this molecule
            rinte = INTEN_VALUES[0]
            for j in range(len(INTEN_GROUP)):  # TODO: baaad python, try a more pythonic way..
                if mol in INTEN_GROUP[j]:
                    rinte = INTEN_VALUES[j]
            rinte = random.uniform(rinte[0], rinte[1])

            for lin in linlist:
                counter += 1
                trans_temp = lin[5]
                flux = np.exp(-abs(trans_temp - self.temp) / self.temp) * rinte
                if flux < cutoff: # TODO: astropy units!
                    continue
                freq = (1 + self.z) * lin[3]  # TODO: astropy unit... Catalog in Mhz
                #self.log.write('      |- Projecting ' + str(lin[2]) + ' (' + str(lin[1]) + ') around ' + str(
                #    freq) + ' Mhz, at ' + str(temp) + ' K\n')
                self._draw_func(cube,flux,freq,cutoff)
                used = True
                # TODO: generate a table: example:All the next commented lines were for generating a table: 
                #arr_code.append(self.comp_name + '-r' + str(self.alpha) + '-d' + str(self.delta) + "-l" + str(counter))
                #arr_mol.append(mol)
                #arr_temp.append(temp)
                #arr_chname.append(str(lin[2]))
                #arr_rest_freq.append(str(lin[3]))
                #arr_rad_vel.append(self.rv)
                #arr_fwhm.append(self.spe_form[1])
        dba.disconnect()
        if not used:
            return
        #hduT = fits.PrimaryHDU()
        #hduT.data = T;
        #hduG = fits.PrimaryHDU()
        #hduG.data = G;
        #tbhdu = fits.new_table(fits.ColDefs([
        #    fits.Column(name='line_code', format='60A', array=arr_code),
        #    fits.Column(name='mol', format='20A', array=arr_mol), \
        #    fits.Column(name='chname', format='40A', array=arr_chname), \
        #    fits.Column(name='rest_freq', format='D', array=arr_rest_freq), \
        #    fits.Column(name='rad_vel', format='D', array=arr_rad_vel), \
        #    fits.Column(name='fwhm', format='D', array=arr_fwhm), \
        #    fits.Column(name='temp', format='D', array=arr_temp)]))
        #cube._add_HDU(hduT)
        #cube._add_HDU(hduG)
        #cube._add_HDU(tbhdu)



