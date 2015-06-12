### Core ASYDO Classes ###

class IMCM(Component):
    """ Interstellar Molecular Cloud Model """

    def __init__(self, log, dbpath, mol_list, temp, spa_form, spe_form, z_grad, z_base=0.0, abun_max=10 ** -5,
                 abun_min=10 ** -6, abun_CO=1.0, iso_abun=default_iso_abundance):
        Component.__init__(self, log, z_base)
        self.spa_form = spa_form
        self.spe_form = spe_form
        self.z_grad = z_grad
        self.dbpath = dbpath
        self.temp = temp
        self.intens = dict()
        for mol in mol_list.split(','):
            abun = random.uniform(abun_min, abun_max)
            if mol in ('COv=0', '13COv=0', 'C18O', 'C17O', '13C18O'):
                abun += abun_CO
            for iso in iso_abun:
                if iso in mol:
                    abun *= iso_abun[iso]
            self.intens[mol] = abun

    def change_intensities(self, intens):
        '''User defined dictionary in the form {molecule: intensity}'''
        self.intens = intens;

    def info(self):
        return "mol_list = " + str(self.intens.keys()) + " @ spa_form=" + str(self.spa_form) + ", spe_form=" + str(
            self.spe_form) + ", z=" + str(self.z) + ", grad=" + str(self.z_grad)

    def project(self, cube):
        arr_code = []
        arr_mol = []
        arr_chname = []
        arr_rest_freq = []
        arr_rad_vel = []
        arr_fwhm = []
        arr_temp = []
        self.log.write('   --+ Generating template image\n')  # TODO More info
        T, (ybord, xbord) = gen_surface(self.spa_form, self.alpha, self.delta, cube.alpha_axis, cube.delta_axis)
        if isinstance(T, bool):
            return
        G = gen_gradient(self.z_grad, self.alpha, self.delta, cube.alpha_axis, cube.delta_axis, ybord, xbord)
        self.log.write('   --+ Generating template line distribution\n')  #TODO More info
        self.log.write('   --+ Loading and correcting lines\n')
        dba = db.lineDB(self.dbpath)
        dba.connect()
        freq_init_corr = cube.freq_border[0] / (1 + self.z)
        freq_end_corr = cube.freq_border[1] / (1 + self.z)
        counter = 0
        used = False
        for mol in self.intens:
            # For each molecule specified in the dictionary
            # load its spectral lines

            linlist = dba.getSpeciesLines(mol, freq_init_corr,
                                          freq_end_corr)  # Selected spectral lines for this molecule
            rinte = inten_values[0]
            for j in range(len(inten_group)):  # TODO baaad python...
                if mol in inten_group[j]:
                    rinte = inten_values[j]
            rinte = random.uniform(rinte[0], rinte[1])

            for lin in linlist:
                counter += 1
                trans_temp = lin[5]
                temp = np.exp(-abs(trans_temp - self.temp) / self.temp) * rinte
                if temp < 3 * cube.noise:
                    continue
                freq = (1 + self.z) * lin[3]  # Catalogs must be in Mhz
                self.log.write('      |- Projecting ' + str(lin[2]) + ' (' + str(lin[1]) + ') around ' + str(
                    freq) + ' Mhz, at ' + str(temp) + ' K\n')
                for xp in range(xbord[0], xbord[1]):
                    for yp in range(ybord[0], ybord[1]):
                        freq = freq_correct(lin[3],self.rv + G[yp - ybord[0], xp - xbord[0]])
                        L, Lbord = gen_line(self.spe_form, freq, cube.freq_axis)
                        if isinstance(L, bool):
                            continue
                        cube.data[Lbord[0]:Lbord[1] + 1, yp, xp] = cube.data[Lbord[0]:Lbord[1] + 1, yp, xp] + T[yp - ybord[0], xp - xbord[0]] * temp* L
                        used = True

                arr_code.append(self.comp_name + '-r' + str(self.alpha) + '-d' + str(self.delta) + "-l" + str(counter))
                arr_mol.append(mol)
                arr_temp.append(temp)
                arr_chname.append(str(lin[2]))
                arr_rest_freq.append(str(lin[3]))
                arr_rad_vel.append(self.rv)
                arr_fwhm.append(self.spe_form[1])
        dba.disconnect()
        if not used:
            return
        hduT = fits.PrimaryHDU()
        hduT.data = T;
        hduG = fits.PrimaryHDU()
        hduG.data = G;
        tbhdu = fits.new_table(fits.ColDefs([
            fits.Column(name='line_code', format='60A', array=arr_code),
            fits.Column(name='mol', format='20A', array=arr_mol), \
            fits.Column(name='chname', format='40A', array=arr_chname), \
            fits.Column(name='rest_freq', format='D', array=arr_rest_freq), \
            fits.Column(name='rad_vel', format='D', array=arr_rad_vel), \
            fits.Column(name='fwhm', format='D', array=arr_fwhm), \
            fits.Column(name='temp', format='D', array=arr_temp)]))
        cube._add_HDU(hduT)
        cube._add_HDU(hduG)
        cube._add_HDU(tbhdu)



