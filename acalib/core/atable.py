import collections
from astropy.table.table import Table as AstropyTable
from astropy.table.table import Column
import json
from copy import deepcopy
from os.path import isfile
from astropy.io import fits
import numpy as np


class ATable(AstropyTable):
    def __init__(self,names,dtype,meta=None):
        AstropyTable.__init__(self, names=names,dtype=dtype)

    def __iadd__(self, other):
        self.add_row(other)
        return self

    @staticmethod
    def import_from(importer, column_types=None):
        table = ATable()

        # load the data to the importer.
        importer.load()

        # read all columns from the importer.
        columns = importer.read_colnames()
        importer.columns = columns

        # return None if the columns aren't valid.
        if not columns:
            return

        # add every column to the table
        for i, column in enumerate(columns):
            # the default column type.
            dtype = 'string'

            # if the user specified a column type, apply it.
            if column_types is not None:
                dtype = column_types[i]

            column = Column([], name=column, dtype=dtype, length=10)
            table.add_column(column)

        # iterable over the rows and add them to the table.
        for row in importer.read_rows():
            table += tuple(row)

        return table

    def filled(self, fill_value=None):
        """Return a copy of self, with masked values filled.

        If input ``fill_value`` supplied then that value is used for all
        masked entries in the table.  Otherwise the individual
        ``fill_value`` defined for each table column is used.

        Parameters
        ----------
        fill_value : str
            If supplied, this ``fill_value`` is used for all masked entries
            in the entire table.

        Returns
        -------
        filled_table : Table
            New table with masked values filled
        """
        if self.masked:
            data = [col.filled(fill_value) for col in six.itervalues(self.columns)]
        else:
            data = self
        return AstropyTable(data, meta=deepcopy(self.meta))


    def get_hdu(self,primary=False):
         if primary:
            raise NotImplementedError("FITS Format do now support tables as primary HDU! You can set primary = None")        
         else:
            hdu=fits.BinTableHDU.from_columns(np.array(self.filled()))
         if self.meta is not None:
            for k, v in self.meta.iteritems():
               hdu.header[k] = v
         return hdu


class AbstractImporter(object):
    def __init__(self, source):
        if source is None:
            raise ValueError("Source cannot be None")

        self.columns = []
        self.source = source

    def load(self):
        """
        Loads the data needed to work with.
        """
        raise NotImplementedError("Not implemented!")

    def read_colnames(self):
        """
        Read column names to construct the table.
        """
        raise NotImplementedError("Not implemented!")

    def read_rows(self):
        """
        Read each row. Should return an iterable (ideally a generator).
        """
        raise NotImplementedError("Not implemented!")


class JsonImporter(AbstractImporter):
    def load(self):
        if isfile(self.source):
            self.data = json.load(open(self.source))
        else:
            self.data = json.loads(self.source)

    def walk(self, indices, data=None):
        sub_data = data if data is not None else self.data
        for c in indices:
            sub_data = sub_data[c]
        return sub_data

    def keys_from(self, indices=[], minus=[]):
        """
        Returns the keys in a given map. As an example:

            [{a: 1, b: 2}, {a: 3, b: 4}]

        You can get the keys a, b, with:

            keys_from([0])

        Which goes to the first element and retrieves all the keys.
        """
        json_keys = self.walk(indices).keys()
        for m in minus:
            json_keys.remove(m)
        return json_keys

    def values_from(self, indices, vdata=None):
        """
        Similar as keys_from, retrieves all rows that are inside the indices given. Example:

            [{items: [{a: 1, b: 2}, {a: 3, b: 4}]}, ...]

        You can get the values from 'items' as:

            values_from([0, 'items'])

        Which retrieves a generator with: [(1, 2), (3, 4)]
        """

        if vdata is None:
            vdata = self.data

        data = self.walk(indices[0], data=vdata)

        values = []
        if len(indices) > 1:
            for d in data:
                d_values = self.values_from(indices[1:], d)

                for val in d_values:
                    values.append(val)

        else:
            for row in data:
                col_values = []
                for column in self.columns:
                    col_values.append(row[column])

                values.append(col_values)

        return values



# the client of the ATable importer should extend the AbstractImporter class
# in order to tell the ATable how to read values. In this case we use a predefined
# JsonImporter that helps a little.
class TestJsonImporter(JsonImporter):
    def read_colnames(self):
        return self.keys_from(["table", 0, "elements", 1], minus=['electrons', 'molar'])

    def read_rows(self):
        return self.values_from((["table"], ['elements']))


if __name__ == "__main__":
    table = ATable.import_from(TestJsonImporter('test.json'))
    print table
