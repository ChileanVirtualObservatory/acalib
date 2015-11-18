from unittest import TestCase
from Client import Client

class TestClient(TestCase):

    def setUp(self):
        self.client = Client("http://127.0.0.1:5000/slap/")

    def test_parse_query(self):
        query = \
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
                    },
                "NRAO Recommended":
                    {
                        "type":     "equality",
                        "value":    True
                    }
            }
        expected_output = "REQUEST=queryData&WAVELENGHT=1.0/1.2,3&NRAO_RECOMMENDED=True"
        assert self.client.parse_query(query) == expected_output




