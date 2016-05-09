from unittest import TestCase
import sys

sys.path.append("../../acalib/io/sla/")
from SlapClient import SlapClient


class TestSlapClient(TestCase):
    def setUp(self):
        self.Client = SlapClient("https://find.nrao.edu/splata-slap/slap",slap_version=1.0)

    def test_query(self):
        input_basic = {
            "wavelenght" :
                     {
                         "gte": 0.1,
                         "lte": 0.5
                     }
                 }

        self.fail()

    def test_query_service(self):
        self.fail()


    def test_query_fields(self):
        self.fail()

    def test_equality(self):
        input_1 = 0.01
        input_2 = "Dihydrogen oxide"
        self.assertEqual(self.Client._SlapClient__equality(input_1), str(input_1))
        self.assertEqual(self.Client._SlapClient__equality(input_2), str(input_2))

    def test_range(self):
        input_1 = {"gte": 0.1,
                   "lte": 0.5}
        input_2 = {"gte": 0.1}
        input_3 = {"lte": 0.5}

        expected_output_1 = "0.1/0.5"
        expected_output_2 = "0.1/"
        expected_output_3 = "/0.5"

        self.assertEqual(self.Client._SlapClient__range(input_1), expected_output_1)
        self.assertEqual(self.Client._SlapClient__range(input_2), expected_output_2)
        self.assertEqual(self.Client._SlapClient__range(input_3), expected_output_3)

    def test_list(self):
        input_1 = [0.1, 0.4, 5]
        input_2 = [{"gte": 0.1, "lte": 0.5}, {"gte": 0.1}, {"lte": 0.5}]
        input_3 = [0.01, {"gte": 0.1, "lte": 0.5}, {"lte": 0.5}, 5]

        expected_output_1 = "0.1,0.4,5"
        expected_output_2 = "0.1/0.5,0.1/,/0.5"
        expected_output_3 = "0.01,0.1/0.5,/0.5,5"

        self.assertEqual(self.Client._SlapClient__list(input_1), expected_output_1)
        self.assertEqual(self.Client._SlapClient__list(input_2), expected_output_2)
        self.assertEqual(self.Client._SlapClient__list(input_3), expected_output_3)

    def test_chemical_query(self):
        self.fail()
