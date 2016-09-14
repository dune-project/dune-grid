from __future__ import absolute_import, division, print_function, unicode_literals

from .create import module

def create(constructor, **parameters):

    dimgrid = parameters.pop["dimgrid"]
    typeName = "Dune::AlbertaGrid< " + str(dimgrid) + " >"
    includes = [ "dune/grid/albertagrid.hh", "dune/grid/albertagrid/dgfparser.hh" ]
    gridModule = module(includes, typeName)
    return gridModule.LeafGrid( gridModule.reader(constructor) )

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
