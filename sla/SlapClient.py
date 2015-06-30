from astropy import units as u
from astropy.io.votable import parse_single_table
from astropy.io.votable.tree import Table as pTable
from astropy.io.votable.tree import Field as pField
import urllib2
import sys
from os import remove

astropy_equivalent = {
    "m"                     :   "metros",
    "MHz"                   :   u.MHz,
    "molecule type"        :   None,
    "sijmu2"                :   None,
    "aij"                   :   None,
    "lowerstateenergy"      :   None,
    "lowerstateenergyK"     :   u.K,
    "upperstateenergy"      :   None,
    "upperstateenergyK"     :   u.K,
    "frequency recommended" :   None
}


class LineDatabase:
    def __init__(self, line_source):
        self.line_source = line_source
        self.__generate_params_cache()

    def __is_astropy(self, variable):
        if not isinstance(variable,u.quantity.Quantity):
            raise ValueError("Must be a Astropy Quantity")

    def __get_specie_lines(self, table):
        pass
    def __execute_query(self, query):
        curl = self.line_source + query.encode('utf-8')
        print curl
        req = urllib2.Request(curl)
        response = urllib2.urlopen(req)
        return response.read()

    def __generate_params_cache(self):
        query = '?REQUEST=queryData&WAVELENGTH=1'
        data = self.__execute_query(query)

        location = 'fields.xml'
        f = open(location, 'w')
        f.write(data)
        f.close()

        table = parse_single_table(location)

        remove(location)

        self.fieldsCache = {}

        for field in table.fields:
            self.fieldsCache[field.name] = self.__assign_unit(field.name,field.description,field.datatype, field.unit)

        print self.fieldsCache


    def __assign_unit(self, name, description, datatype, unit):
        new_unit = None
        if datatype == "int" or datatype == "double":
            if unit:
                new_unit = astropy_equivalent[str(unit)]
            else:
                new_unit = astropy_equivalent[str(name)]

        return (str(description), str(datatype), new_unit)





    def __save_n_parse(self, data):
        location ='temp.xml'
        f = open(location, 'w')
        f.write(data)
        f.close()

        table = parse_single_table(location)

        if isinstance(table, pTable):
            # print type(table.array)
            # for i in table.fields:
            #     if isinstance(i, pField):
            #         print(i.name)
            #         print(i.description)
            #         print(i.datatype)
            #         print(i.unit)
            #         print("\n")
            print type(table.array[1])
           # print table.array[1]
            print table.array[1]['wavelength']








    def getSpecieLines(self, specie, lower_frequency, upper_frequency):
        """
        Method that accepts range of frequencies and a molecule, and return all available spectral lines for that molecule in the defined range
        This is the "loose" Version, meaning that the molecule name match is not strict

        :param specie: The Molecule to search, for example "O2", "CO2"
        :type specie: String
        :param lowerFrequency: Lower boundary for the frequency range
        :type lowerFrequency: Astropy frequency, for example 495 * astropy.units.HZ
        :param upperFrequency: Upper boundary for the frequency range
        :type upperFrequency: Astropy frequency, for example 495 * astropy.units.HZ
        :return: An arrary containing all the lines for the given molecule. The data for each lines is in astropy Unit format
        :rtype: Array containig astropy units
        """
        self.__is_astropy(lower_frequency)
        self.__is_astropy(upper_frequency)

        lower_wavelength = upper_frequency.to(u.m,equivalencies=u.spectral())
        upper_wavelength = lower_frequency.to(u.m,equivalencies=u.spectral())
        data = '?REQUEST=queryData&WAVELENGTH=' + str(lower_wavelength.value) + '/' + str(upper_wavelength.value)

        votable = self.__execute_query(data)
        self.__save_n_parse(votable)





    def getLots(self, upper_wavelength):
        data = '?REQUEST=queryData&WAVELENGTH=/' + str(upper_wavelength.to(u.m).value)
        curl = self.line_source + data.encode('utf-8')
        print curl

        req = urllib2.Request(curl)
        response = urllib2.urlopen(req)
        votable = response.read()
        location = "all-0-to-" + str(upper_wavelength) + '.xml'
        f = open(location, 'w')
        f.write(votable)
        f.close()

    def getMoleculeList(self,lowerFrequency, upperFrequency):
        pass

    def getSpeciesList(self,chem_name,lowerFrequency, upperFrequency):
        pass



test = LineDatabase("http://find.nrao.edu/splata-slap/slap")
#test.getSpecieLines("a",1 * u.GHz,1.11 * u.GHz)

