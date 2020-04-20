from dune.common.checkconfiguration import assertHave, ConfigurationError
from dune.typeregistry import generateTypeName

def onedGrid(constructor):
    from .grid_generator import module, getDimgrid

    typeName = "Dune::OneDGrid"
    includes = ["dune/grid/onedgrid.hh", "dune/grid/io/file/dgfparser/dgfoned.hh"]
    gridModule = module(includes, typeName)

    return gridModule.reader(constructor).leafView

def yaspGrid(constructor, dimgrid=None, coordinates="equidistant", ctype="double"):
    from dune.generator import Constructor
    from .grid_generator import module, getDimgrid

    if not dimgrid:
        dimgrid = getDimgrid(constructor)
    if coordinates == "equidistant":
        coordinates = "Dune::EquidistantOffsetCoordinates"
    elif coordinates == "tensorproduct":
        coordinates = "Dune::TensorProductCoordinates"
    coordinates, _ = generateTypeName(coordinates, ctype, dimgrid)
    typeName, includes = generateTypeName("Dune::YaspGrid", dimgrid, coordinates)
    includes += ["dune/grid/yaspgrid.hh", "dune/grid/io/file/dgfparser/dgfyasp.hh"]

    ctor = Constructor([
              'Dune::FieldVector<'+ctype+', '+str(dimgrid)+'> lowerleft',
              'Dune::FieldVector<'+ctype+', '+str(dimgrid)+'> upperright',
              'std::array<int, '+str(dimgrid)+'> elements',
              'std::array<bool, '+str(dimgrid)+'> periodic',
              'int overlap'],
              ['std::bitset<'+str(dimgrid)+'> periodic_;',
               'for (int i=0;i<'+str(dimgrid)+';++i) periodic_.set(i,periodic[i]);',
               'return new DuneType(lowerleft,upperright,elements,periodic_,overlap);'],
              ['"lowerleft"_a', '"upperright"_a', '"elements"_a', '"periodic"_a', '"overlap"_a'])

    gridModule = module(includes, typeName, ctor)

    try:
        lowerleft  = constructor.lower
        upperright = constructor.upper
        elements   = constructor.division
        periodic   = constructor.param.get("periodic", [False,]*dimgrid)
        overlap    = constructor.param.get("overlap", 0)
        return gridModule.HierarchicalGrid(lowerleft,upperright,elements,periodic,overlap).leafView
    except AttributeError:
        return gridModule.reader(constructor).leafView

grid_registry = {
        "OneD"       : onedGrid,
        "Yasp"       : yaspGrid,
    }

try:
    assertHave("HAVE_ALBERTA")
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
