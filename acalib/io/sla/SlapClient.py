import requests
from astropy.io import votable
from astropy import units as u
import re
import StringIO

__author__ = "teohoch"


def _check(par, default):
    if not isinstance(par, u.Quantity):
        par = par * default
    return par

class SlapClient(object):
    def __init__(self, slap_service, slap_version=1.0):
        """

        :type slap_version: Double
        :type slap_service: String

        """
        self.__slap_service = slap_service
        self.__slap_version = slap_version

    @classmethod
    def query(cls, service, output_format="txt", slap_version=1.0, **kwargs):
        """
        Class method used to query a SLAP Service, with the contrains following the following scheme:
        for equality contraint, do parameter=value
        for range constraint, do parameter={"gte" : minimum_value, "lte" : maximum_value}.
        Note: for open ranges (i.e 5 to infinity), do not include the key for the open limit ( parameter={"gte" : 5} )
        for a list of constrints to one parameter, just store all the constraints in a list
        (i.e. parameter=[1,2,{"gte" : 5}]
        :param service: String. The Slap service to use
        :param output_format: String. Format in which to output the data. Current available choices are raw text ("txt")
         or votable ("votable"). Default is "txt"
        :param slap_version: Slap version of the service. Default's 1.0
        :param kwargs: The constrains are included here, as explained in the description
        :return: Depending on the output format selected, the output will change, but the contents will always be the
        reply to the query.
        """
        query_params = SlapClient.__parser(slap_version, kwargs)
        if query_params:
            request = requests.get(service, params=query_params)
            return cls.__output_formater(request.text, output_format)
        else:
            raise ValueError("The query is not valid according to the SLA Protocol")

    @classmethod
    def chemical_query(cls, service, chemical_element, output_format="txt", minimum=None, maximum=None,
                       slap_version=1.0):
        """
        Class method used to get all the spectral lines concerning to a certain chemical, within a spectral range
        Take into consideration that constraining results by Chemical Element is not a requirement of the SLAP protocol,
        so use discretion. The splatalogue and Slap-a-Chivo services nevertheless support this feature.

        :param service: String. The Slap service to use
        :param chemical_element: String. The chemical/s in question to search for.
        :param output_format: Format in which to output the data. Current available choices are raw text ("txt") or votable ("votable"). Default is "txt"
        :param minimum: Number. Lower limit for the spectral range. If astropy units are not used, it will assume the limit is in Hz
        :param maximum: Number. Upper limit for the spectral range. If astropy units are not used, it will assume the limit is in Hz
        :param slap_version: Number. The SLAP version to use. Defaults to 1.0
        :return: Depending on the output format selected, the output will change, but the contents will always be the
        reply to the query.
        """
        wave = dict()
        min_temporal = min_final = None
        max_temporal = max_final = None

        if minimum:
            min_temporal = ((_check(minimum, u.Hz)).to(u.m, equivalencies=u.spectral())).value
        if maximum:
            max_temporal = ((_check(maximum, u.Hz)).to(u.m, equivalencies=u.spectral())).value
        if max_temporal and min_temporal:
            max_final = max_temporal if (max_temporal >= min_temporal) else min_temporal
            min_final = max_temporal if (max_temporal >= min_temporal) else min_temporal
            wave["lte"] = str(max_final)
            wave["gte"] = str(min_final)
        elif max_temporal:
            max_final = max_temporal if (max_temporal >= min_temporal) else min_temporal
            wave["lte"] = str(max_final)
        elif min_temporal:
            min_final = max_temporal if (max_temporal >= min_temporal) else min_temporal
            wave["gte"] = str(min_final)
        else:
            raise ValueError("Either the minimum, the maximum or both must be valid for the query to be valid")

        data = cls.query(service, chemical_element, slap_version=slap_version, wavelenght=wave,
                         CHEMICAL_ELEMENT=chemical_element)
        return cls.__output_formater(data.text, output_format)

    @classmethod
    def __parser(cls, slap_version, parameters):

        query_params = {"VERSION": slap_version,
                        "REQUEST": "queryData"}

        for constrain_name, constrain in parameters.iteritems():
            if constrain_name.lower() == "wavelength":
                valid = True
            if isinstance(constrain, dict):  # range
                query_params[constrain_name.upper()] = cls.__range(constrain)
            elif isinstance(constrain, list):  # list of constrains
                query_params[constrain_name.upper()] = cls.__list(constrain)
            else:  # equality
                query_params[constrain_name.upper()] = cls.__equality(constrain)

        return query_params if valid else None

    @classmethod
    def __output_formater(cls, data, format):
        if format.lower() == "txt":
            return data
        elif format.lower() == "votable":
            return cls.__to_votable(data)
        else:
            raise ValueError("The Output value is not Valid")

    @classmethod
    def __list(cls, params):
        output = []
        for constrain in params:
            if isinstance(constrain, dict):  # range
                output.append(cls.__range(constrain))
            elif isinstance(constrain, list):  # list of constrains
                raise ValueError("Invalid list within list detected")
            else:  # equality
                output.append(cls.__equality(constrain))
        return ",".join(output)

    @staticmethod
    def __range(params):
        gte = params.has_key("gte")
        lte = params.has_key("lte")
        if gte and lte:
            return str(params["gte"]) + "/" + str(params["lte"])
        elif gte:
            return str(params["gte"]) + "/"
        elif lte:
            return "/" + str(params["lte"])
        else:
            raise ValueError("The Range dictionary must contain eithe 'gte' (greater than equal),"
                             " 'lte' (lower than equal), or both")

    @staticmethod
    def __equality(params):
        return str(params)

    @staticmethod
    def __to_votable(data):
        temp = StringIO.StringIO()
        temp.write(data)
        return votable.parse(temp)

    def query_service(self, output_format="txt", **kwargs):
        """
        Method used to query a SLAP Service, with the contrains following the following scheme:
        for equality contraint, do parameter=value
        for range constraint, do parameter={"gte" : minimum_value, "lte" : maximum_value}.
        Note: for open ranges (i.e 5 to infinity), do not include the key for the open limit ( parameter={"gte" : 5} )
        for a list of constrints to one parameter, just store all the constraints in a list
        (i.e. parameter=[1,2,{"gte" : 5}]
        :param output_format: String. Format in which to output the data. Current available choices are raw text ("txt")
         or votable ("votable"). Default is "txt"
        :param kwargs: The constrains are included here, as explained in the description
        :return: Depending on the output format selected, the output will change, but the contents will always be the
        reply to the query.
        :param kwargs:
        :return:
        """
        return self.query(self.__slap_service, slap_version=self.__slap_version, **kwargs)

    def query_fields(self):
        """
        Returns the fields (columns) returned by the SLAP Service for an standard query
        :return: A list of dictionaries, each one containing the information of one field.
        Common keys of these dictionaries are "name", description, "ID"
        """
        basic = SlapClient.query(self.__slap_service, self.__slap_version,
                                 wavelength=1.0)  # Splatalogue can't deal with wavelenght == 0 XD!
        # Use Regex <FIELD (.*?)>.*?<DESCRIPTION>(.*?)<\/DESCRIPTION>.*?<\/FIELD>
        # https://regex101.com/r/pL7pY1/1
        p = re.compile(ur'<FIELD (.*?)>.*?<DESCRIPTION>(.*?)<\/DESCRIPTION>.*?<\/FIELD>', re.UNICODE | re.DOTALL)

        # Use Regex (.*?)="(.*?)".*?
        # https://regex101.com/r/qW8iX8/1
        p2 = re.compile(ur'(.*?)="(.*?)".*?', re.UNICODE | re.DOTALL)

        result = re.findall(p, basic)
        fields = map(lambda s: [re.findall(p2, s[0]), s[1]], result)
        a = list()
        for i in fields:
            dict_field = (dict((x[0], x[1]) for x in i[0]))
            dict_field.update({"description": i[1]})
            a.append(dict_field)
        return a


if __name__ == "__main__":
    service = "https://find.nrao.edu/splata-slap/slap"
    client = SlapClient(service)
    data = client.query_service(wavelength={"gte": 0.00260075, "lte": 0.00260080})
    print data
