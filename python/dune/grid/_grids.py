from dune.common.checkconfiguration import assertCMakeHave, ConfigurationError
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
        coordinates_name = "Dune::EquidistantOffsetCoordinates"
    elif coordinates == "tensorproduct":
        coordinates_name = "Dune::TensorProductCoordinates"
    coordinates_type, _ = generateTypeName(coordinates_name, ctype, dimgrid)
    typeName, includes = generateTypeName("Dune::YaspGrid", dimgrid, coordinates_type)
    includes += ["dune/grid/yaspgrid.hh", "dune/grid/io/file/dgfparser/dgfyasp.hh"]

    if coordinates == "equidistant":
        ctor = Constructor(
            [ 'Dune::FieldVector<'+ctype+', '+str(dimgrid)+'> lowerleft',
              'Dune::FieldVector<'+ctype+', '+str(dimgrid)+'> upperright',
              'std::array<int, '+str(dimgrid)+'> elements',
              'std::array<bool, '+str(dimgrid)+'> periodic',
              'int overlap' ],
            [ 'std::bitset<'+str(dimgrid)+'> periodic_;',
              'for (int i=0;i<'+str(dimgrid)+';++i) periodic_.set(i,periodic[i]);',
              'return new DuneType(lowerleft,upperright,elements,periodic_,overlap);' ],
            [ '"lowerleft"_a', '"upperright"_a', '"elements"_a', '"periodic"_a', '"overlap"_a' ])
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
    elif coordinates == "tensorproduct":
        ctor = Constructor(
            [ 'std::array<std::vector<'+ctype+'>, '+str(dimgrid)+'> coords',
              'std::array<bool, '+str(dimgrid)+'> periodic',
              'int overlap' ],
            [ 'std::bitset<'+str(dimgrid)+'> periodic_;',
              'for (int i=0;i<'+str(dimgrid)+';++i) periodic_.set(i,periodic[i]);',
              'return new DuneType(coords,periodic_,overlap);' ],
            [ '"coords"_a', '"periodic"_a', '"overlap"_a' ])
        gridModule = module(includes, typeName, ctor)
        try:
            coords     = constructor.coords
            periodic   = constructor.param.get("periodic", [False,]*dimgrid)
            overlap    = constructor.param.get("overlap", 0)
        except AttributeError:
            coords     = constructor
            periodic   = [False,]*dimgrid
            overlap    = 0
        return gridModule.HierarchicalGrid(coords,periodic,overlap).leafView
    assert(False)

grid_registry = {
        "OneD"       : onedGrid,
        "Yasp"       : yaspGrid,
    }

from dune.packagemetadata import getCMakeFlags
try:
    if not getCMakeFlags()["HAVE_ALBERTA"]:
        raise KeyError
    def albertaGrid(constructor, dimgrid=None, dimworld=None):
        from .grid_generator import module, getDimgrid

        if not dimgrid:
            dimgrid = getDimgrid(constructor)
        if dimworld is None:
            dimworld = dimgrid
        typeName = "Dune::AlbertaGrid< " + str(dimgrid) + " >"
        includes = ["dune/grid/albertagrid.hh", "dune/grid/albertagrid/dgfparser.hh"]
        extraCMake = ["add_dune_alberta_flags(TARGET WORLDDIM {})".format(str(dimworld))]
        gridModule = module(includes, typeName) # , extraCMake=extraCMake)
        return gridModule.reader(constructor).leafView
except KeyError:
    def albertaGrid(constructor, dimgrid=None, dimworld=None):
        print("""
Alberta was not found during the configuration of dune-grid.
To use 'albertaGrid' install the alberta package first and then reinstall the dune-grid package.
""")
grid_registry["Alberta"] = albertaGrid

from dune.packagemetadata import getCMakeFlags
try:
    if not getCMakeFlags()["HAVE_DUNE_UGGRID"]:
        raise KeyError
    def ugGrid(constructor, dimgrid=None, **parameters):
        from .grid_generator import module, getDimgrid

        if not dimgrid:
            dimgrid = getDimgrid(constructor)
        typeName = "Dune::UGGrid< " + str(dimgrid) + " >"
        includes = ["dune/grid/uggrid.hh", "dune/grid/io/file/dgfparser/dgfug.hh"]
        gridModule = module(includes, typeName)

        return gridModule.reader(constructor).leafView
except KeyError:
    def ugGrid(constructor, dimgrid=None, **parameters):
        print("""
UGGrid was not found during the configuration of dune-grid.
To use 'ugGrid' install the dune-uggrid package first and then reinstall the dune-grid package.
""")
grid_registry["UG"] = ugGrid

if __name__ == "__main__":
    import doctest
