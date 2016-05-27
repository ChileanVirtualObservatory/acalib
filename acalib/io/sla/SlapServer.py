import sqlite3 as lite
from SlapClient import SlapClient

SqlEquivalent = {
    "char": "TEXT",
    "double": "DOUBLE",
    "int": "INT",
    "boolean": "BOOLEAN"
}


class SlapServer:
    # slap_serv = 'https://find.nrao.edu/splata-slap/slap'


    def __init__(self, filename):
        self.connected = False
        self.name = filename
        try:
            self.pointer = lite.connect(self.name + ".sqlite")
            self.connected = True
        except lite.Error, e:
            print "Error %s:" % e.args[0]

    def __del__(self):
        if self.pointer:
            self.pointer.close()
            self.connected = False

    def __execute_sql(self, sentence):
        resp = self.pointer.execute(sentence)
        return resp.fetchall()

    def load_from_slap_service(self, service, slap_version=1.0, **kwargs):
        Client = SlapClient(service, slap_version)
        self.__fields = tuple(Client.query_fields())

        # Create a table for the lines
        base_command = "CREATE TABLE Lines (ID INT PRIMARY KEY NOT NULL"
        for field in self.__fields:
            base_command = " ".join(
                (base_command, ",", field["ID"].replace(" ", "_"), SqlEquivalent[field["datatype"]]))
        base_command += ")"
        self.__execute_sql(base_command)

        # insert the Data into the SQL  Database.




if __name__ == "__main__":
    service = "https://find.nrao.edu/splata-slap/slap"
    server = SlapServer("Test.SQL")
    server.load_from_slap_service(service)
