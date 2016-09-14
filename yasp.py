from __future__ import absolute_import, division, print_function, unicode_literals

from .create import module

def create(constructor, dimgrid,
        coordinates="Dune::EquidistantCoordinates", ctype="double",
        **parameters):

    if coordinates == "equidistant":
        coordinates = "=Dune::EquidistantCoordinates"
    elif coordinates == "tensorproduct":
        coordinates = "=Dune::TensorProductCoordinates"
    typeName = "Dune::YaspGrid< " + str(dimgrid) + ", " +\
               coordinates + "< " + ctype + ", " + str(dimgrid) + " > >"
    includes = ["dune/grid/yaspgrid.hh", "dune/grid/io/file/dgfparser/dgfyasp.hh"]
    gridModule = module(includes, typeName)
    return gridModule.LeafGrid( gridModule.reader(constructor) )

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
