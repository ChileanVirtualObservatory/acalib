from numpy import random
from synthetic import db
from synthetic.vu import Component
import core.flux as flx
import astropy.units as u
import numpy as np
from astropy import log
import urllib
import shutil
import os.path
from astropy.io import fits
import core.parameter as par

#INTEN_GROUP = [('default'), ('COv=0'), ('13COv=0'), ('HCO+, HC3N, CS, C18O, CH3OH, N2H, HDO')]
#INTEN_VALUES = [[0.1, 2], [20, 60], [5, 20], [1, 10]]
#DEFAULT_ISO_ABUND = {'13C': 1. / 30., '18O': 1. / 60., '17O': 1. / 120., '34S': 1. / 30., '33S': 1. / 120.,
#                         '13N': 1. / 30., 'D': 1./40.}
#DEFAULT_ABUND_RANGE=[10**-5,10**-6]
#GAUSS_STRINGS = ["Gaussian", "gaussian", "Gauss", "gauss", "normal", "Normal"]
#DEFAULT_CO_ABUND=1.0


DEFAULT_DBPATH= 'ASYDO'


class IMC(Component):
    """ Interstellar Molecular Core """

    def __init__(self, mol_list, temp,dbpath=DEFAULT_DBPATH,equiv=u.doppler_radio):
        Component.__init__(self)
        self.equiv=equiv
        self.dbpath = dbpath
        self.temp = temp
        self.mol_list = mol_list


    def _draw(self,cube,flux,freq,cutoff):
       # TODO      
       pass

    def info(self):
       return ""

    def project(self, cube, cutoff):
        #TODO Make all this with astropy units from the call functions
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
        #print "fwin",fwin
        cor_fwin =  np.array(fwin/(1 + self.z))*u.Hz
        cor_fwin = cor_fwin.to(u.MHz).value
        #print "cor_fwin",cor_fwin
        counter = 0
        used = False
        for mol in self.mol_list:
            # For each molecule specified in the dictionary
            # load its spectral lines
            linlist = dba.getSpeciesLines(mol, cor_fwin[0], cor_fwin[1])  # Selected spectral lines for this molecule
            #rinte = INTEN_VALUES[0]
            #for j in range(len(INTEN_GROUP)):  # TODO: baaad python, try a more pythonic way..
            #    if mol in INTEN_GROUP[j]:
            #        rinte = INTEN_VALUES[j]
            rinte = random.uniform(self.mol_list[mol][0], self.mol_list[mol][1])*u.Jy/u.beam
             
            for lin in linlist:
                counter += 1
                trans_temp = lin[5]*u.K
                flux = np.exp(-abs(trans_temp - self.temp) / trans_temp) * rinte
                #print trans_temp, self.temp, flux, rinte
                freq = (1 + self.z) * lin[3]*u.MHz  # TODO: astropy 
                #print flux, cutoff
                if flux < cutoff: # TODO: astropy units!
                    log.info('    - Discarding ' + str(lin[1]) + ' at freq=' + str(freq) + '('+str(lin[3]*u.MHz)+') because I='+str(flux)+' < '+str(cutoff))
                    continue
                log.info('   - Projecting ' + str(lin[2]) + ' (' + str(lin[1]) + ') at freq=' + str(freq) + '('+str(lin[3]*u.MHz)+') intens='+ str(flux)+ ' '+cube.unit.to_string())
                self._draw(cube,flux,freq,cutoff)
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


class GaussianIMC(IMC):
    def __init__(self,mol_list, temp, offset,std, angle, fwhm, gradient, dbpath=DEFAULT_DBPATH,equiv=u.doppler_radio):
        IMC.__init__(self,mol_list, temp,dbpath,equiv)
        self.offset=par.to_deg(offset)
        self.std= std
        self.angle = angle
        self.fwhm = fwhm
        self.gradient = gradient
   
    def _draw(self,cube,flux,freq,cutoff):
       new_pos=self.pos + self.offset
       (mu,P)=flx.clump_to_gauss(new_pos,self.std,self.angle,freq,self.fwhm,self.gradient)
       #print "mu",mu
       (mcub,lower,upper)=flx.create_gauss_flux(cube,mu,P,flux,cutoff)
       cube.add_flux(mcub,lower,upper)

    def info(self):
       return "species = " + str(self.mol_list) + " temp = "+str(self.temp)+" offset"+str(self.offset)+" std = "+str(self.std)+" angle ="+str(self.angle)+" fwhm ="+str(self.fwhm)+"  gradient ="+str(self.gradient)


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
