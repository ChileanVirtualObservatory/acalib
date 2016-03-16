from astropy import units as u
import urllib2
import json
from numbers import Number

astropy_equivalent = {
    "m": "metros",
    "MHz": u.MHz,
    "molecule type": None,
    "sijmu2": u.D,
    "aij": u.s ** -1,
    "lowerstateenergy": u.cm ** -1,
    "lowerstateenergyK": u.K,
    "upperstateenergy": u.cm ** -1,
    "upperstateenergyK": u.K,
    "frequency recommended": None
}


class Client():
    def __init__(self, host):
        self.__host = host

    def parse_query(self, query):
        query_params = ["REQUEST=queryData"]

        for key, item in query.items():
            value = ""
            if item["type"] == "range":
                range_value = item["value"]
                # Check that the values in the range are Number types
                if (range_value.has_key("min") and not isinstance(range_value["min"], Number)) or (
                            range_value.has_key("max") and not isinstance(range_value["max"], Number)):
                    raise TypeError("The values in a range MUST be numbers")

                min = str(range_value["min"]) if range_value.has_key("min") else ""
                max = str(range_value["max"]) if range_value.has_key("max") else ""
                value = min + "/" + max

            elif item["type"] == "equality":
                value = item["value"]

            elif item["type"] == "list":
                list_constrains = []
                for constrain in item["value"]:
                    constrain_value = ""
                    if constrain["type"] == "range":
                        range_value = constrain["value"]
                        # Check that the values in the range are Number types
                        if (range_value.has_key("min") and not isinstance(range_value["min"], Number)) or (
                                    range_value.has_key("max") and not isinstance(range_value["max"], Number)):
                            raise TypeError("The values in a range MUST be numbers")

                        min = str(range_value["min"]) if range_value.has_key("min") else ""
                        max = str(range_value["max"]) if range_value.has_key("max") else ""
                        constrain_value = min + "/" + max
                    elif constrain["type"] == "equality":
                        constrain_value = str(constrain["value"])
                    else:
                        raise ValueError("Unrecognized parameter type within parameter list")

                    list_constrains.append(constrain_value)
                value = ",".join(list_constrains)
            else:
                raise ValueError("Unrecognized parameter type")

            query_params.append(key.upper().replace(" ", "_") + "=" + str(value))
        print "&".join(query_params)
        return "&".join(query_params)

    def execute_query(self, query):
        url_with_query = self.__host + "?" + self.parse_query(query)
        req = urllib2.Request(url_with_query)
        response = urllib2.urlopen(req)
        return json.loads(response.read())


if __name__ == "__main__":
    params = \
        {
            "wavelenght":
                {
                    "type": "range",
                    "value":
                        {
                            "min": 1,
                            "max": 1.2
                        }
                }
        }
    params2 = \
        {
            "wavelenght":
                {
                    "type": "list",
                    "value":
                        [
                            {
                                "type": "range",
                                "value":
                                    {
                                        "min": 1.0,
                                        "max": 1.2
                                    }
                            },
                            {
                                "type": "equality",
                                "value": 3
                            }
                        ]
                }
        }

    test = Client("http://127.0.0.1:5000/slap/")
    test.execute_query(params2)
