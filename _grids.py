from __future__ import absolute_import, division, print_function, unicode_literals

from dune.common.checkconfiguration import assertHave, ConfigurationError
from dune.common import generateTypeName

def onedGrid(constructor):
    from .grid_generator import module, getDimgrid

    typeName = "Dune::OneDGrid"
    includes = ["dune/grid/onedgrid.hh", "dune/grid/io/file/dgfparser/dgfoned.hh"]
    gridModule = module(includes, typeName)

    return gridModule.reader(constructor).leafView

def yaspGrid(constructor, dimgrid=None, coordinates="Dune::EquidistantCoordinates", ctype="double"):
    from .grid_generator import module, getDimgrid

    if not dimgrid:
        dimgrid = getDimgrid(constructor)
    if coordinates == "equidistant":
        coordinates = "Dune::EquidistantCoordinates"
    elif coordinates == "tensorproduct":
        coordinates = "Dune::TensorProductCoordinates"
    coordinates, _ = generateTypeName(coordinates, ctype, dimgrid)
    typeName, includes = generateTypeName("Dune::YaspGrid", dimgrid, coordinates)
    includes += ["dune/grid/yaspgrid.hh", "dune/grid/io/file/dgfparser/dgfyasp.hh"]
    # typeName = "Dune::YaspGrid< " + str(dimgrid) + ", " + coordinates + "< " + ctype + ", " + str(dimgrid) + " > >"
    includes = ["dune/grid/yaspgrid.hh", "dune/grid/io/file/dgfparser/dgfyasp.hh"]
    gridModule = module(includes, typeName)

    return gridModule.reader(constructor).leafView

grid_registry = {
        "OneD"       : onedGrid,
        "Yasp"       : yaspGrid,
    }

try:
    assertHave("ALBERTA")
    def albertaGrid(constructor, dimgrid=None):
        from .grid_generator import module, getDimgrid

        if not dimgrid:
            dimgrid = getDimgrid(constructor)
        typeName = "Dune::AlbertaGrid< " + str(dimgrid) + " >"
        includes = ["dune/grid/albertagrid.hh", "dune/grid/albertagrid/dgfparser.hh"]
        gridModule = module(includes, typeName)

        return gridModule.reader(constructor).leafView
    grid_registry["Alberta"] = albertaGrid
except ConfigurationError:
    pass


try:
    assertHave("HAVE_DUNE_UGGRID")
    def ugGrid(constructor, dimgrid=None, **parameters):
        from .grid_generator import module, getDimgrid

        if not dimgrid:
            dimgrid = getDimgrid(constructor)
        typeName = "Dune::UGGrid< " + str(dimgrid) + " >"
        includes = ["dune/grid/uggrid.hh", "dune/grid/io/file/dgfparser/dgfug.hh"]
        gridModule = module(includes, typeName)

        return gridModule.reader(constructor).leafView
    grid_registry["UG"] = ugGrid
except ConfigurationError:
    pass

if __name__ == "__main__":
    import doctest
