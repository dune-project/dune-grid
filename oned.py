from __future__ import absolute_import, division, print_function, unicode_literals

from .grid_generator import module

def create(constructor, **parameters):

    typeName = "Dune::OneDGrid"
    includes = [ "dune/grid/onedgrid.hh", "dune/grid/io/file/dgfparser/dgfoned.hh" ]
    gridModule = module(includes, typeName)
    return gridModule.LeafGrid( gridModule.reader(constructor) )

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
