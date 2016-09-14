from __future__ import absolute_import, division, print_function, unicode_literals

from .create import module

def create(constructor, dimgrid, dimworld=None, elementType=None,
        **parameters):

    refinement = parameters["refinement"]
    if not dimworld:
        dimworld = dimgrid
    if not elementType:
        elementType = parameters.pop("type")
    if refinement == "conforming":
        refinement="Dune::conforming"
    elif refinement == "nonconforming":
        refinement="Dune::nonconforming"

    if not (2<=dimworld and dimworld<=3):
        raise KeyError(\
            "Parameter error in ALUGrid with "+
            "dimworld=" + str(dimworld) + ": " +\
            "dimworld has to be either 2 or 3")
    if not (2<=dimgrid and dimgrid<=3):
        raise KeyError(\
            "Parameter error in ALUGrid with "+
            "dimgrid=" + str(dimgrid) + ": " +\
            "dimgrid has to be either 2 or 3")
    if not (dimgrid==dimworld or dimgrid+1==dimworld):
        raise KeyError(\
            "Parameter error in ALUGrid with "+
            "dimworld=" + str(dimworld) + " and dimgrid=" + str(dimgrid) + ": " +\
            "dimworld has to be either equal to dimgrid or "+
            "dimworld has to be equal to dimgrid+1")
    if refinement=="Dune::conforming" and elementType=="Dune::cube":
        raise KeyError(\
            "Parameter error in ALUGrid with "+
            "refinement=" + refinement + " and type=" + elementType + ": " +\
            "conforming refinement is only available with simplex element type")

    typeName = "Dune::ALUGrid< " + str(dimgrid) + ", " + str(dimworld) + ", " +\
               elementType + ", " + refinement + " >"
    includes = [ "dune/alugrid/grid.hh", "dune/alugrid/dgf.hh" ]
    gridModule = module(includes, typeName)
    return gridModule.LeafGrid( gridModule.reader(constructor) )

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
