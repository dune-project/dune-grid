from __future__ import absolute_import, division, print_function, unicode_literals

from .create import module

def create(constructor, dimgrid,
        coordinates="Dune::EquidistantCoordinates", ctype="double", refinement="Dune::SPIsotropicRefinement",
        **parameters):

    if refinement == "isotropic":
        refinement = "Dune::SPIsotropicRefinement"
    elif refinement == "anisotropic":
        refinement = "Dune::SPAnisotropicRefinement"
    elif refinement == "bisection":
        refinement = "Dune::SPBisectionRefinement"
    elif refinement == "arbitrary":
        refinement = "Dune::SPArbitraryRefinement"
    typeName = "Dune::SPGrid< " + ctype + ", " + str(dimgrid) + ", " + refinement + " >"
    includes = [ "dune/grid/spgrid.hh", "dune/grid/spgrid/dgfparser.hh" ]
    gridModule = module(includes, typeName)
    return gridModule.LeafGrid( gridModule.reader(constructor) )

if __name__ == "__main__":
    import doctest
