import collections
from astropy.table.table import Table as AstropyTable

class ATable(AstropyTable):
    def __init__(self, *names):
        AstropyTable.__init__(self, names=names)

    def __iadd__(self, other):
        self.add_row(vals=other)
        return self


if __name__ == "__main__":
    table = ATable("Col1", "Col2", "Col3")
    table += (1, 2, 3)
    table += (4, 5, 6)

    print table
