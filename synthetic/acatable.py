import collections
from astropy.table.table import Table as AstropyTable

class _AcaTableColumn(object):
    def __init__(self, name):
        self.values = []
        self.name = name

    def __iadd__(self, other):
        if type(other) in (list, tuple):
            for o in other:
                self.values.append(o)
        else:
            self.values.append(other)

        return self

    def clear(self):
        self.values = []

    def as_list(self):
        return list(self.values)

    def get_values(self):
        return self.values

    def __str__(self):
        return "%s: %s" % (self.name, ", ".join(self.values))


class AcaTable_Old(object):
    def __init__(self, *rows):
        self.cols = {}
        self.order = []

        map(self.add_col, rows)

    def add_col(self, row):
        self.cols[row] = _AcaTableColumn(row)
        self.order.append(row)

    def __getitem__(self, item):
        if item not in self.cols:
            self.add_col(item)

        return self.cols[item]

    def __setitem__(self, key, value):
        pass

    def __iadd__(self, other):
        if isinstance(other, collections.Iterable):
            if len(other) == len(self.order):
                for column, value in self.iter_cols(other):
                    column += value
            else:
                raise ValueError("The length of the iterable must be the same of the number of columns.")

        return self

    def iter_cols(self, other=None):
        if other is None:
            for col in self.order:
                yield self.cols[col]
        else:
            for i in range(len(self.order)):
                yield self.cols[self.order[i]], other[i]

    def __iter__(self):
        return self.iter_cols()

    def render_astropy(self):
        """
        Render this table to an AstroPy table.
        :return:
        """
        col_values = []

        for col in self.order:
            col_values.append(self.cols[col].get_values())

        try:
            table = AstropyTable(col_values, names=self.order)
            return table
        except ValueError as e:
            raise ValueError(e.message)

    def __str__(self):
        return self.render_astropy().__str__()

# TODO: THIS CLASS WILL GROW.
# for now, this is just a wrapper with some handy cases of use.
class AcaTable(AstropyTable):
    def __init__(self, *names):
        AstropyTable.__init__(self, names=names)

    def __iadd__(self, other):
        self.add_row(vals=other)
        return self

if __name__ == "__main__":
    table = AcaTable("Col1", "Col2", "Col3")
    table += (1, 2, 3)
    table += (4, 5, 6)

    print table
