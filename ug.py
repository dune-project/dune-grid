from __future__ import absolute_import, division, print_function, unicode_literals

from .grid_generator import module

def create(constructor, **parameters):

    dimgrid = parameters.pop["dimgrid"]
    typeName = "Dune::UGGrid< " + str(dimgrid) + " >"
    includes = [ "dune/grid/uggrid.hh", "dune/grid/io/file/dgfparser/dgfug.hh" ]
    gridModule = module(includes, typeName)
    return gridModule.LeafGrid( gridModule.reader(constructor) )

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
