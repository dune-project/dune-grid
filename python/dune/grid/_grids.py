from dune.common.checkconfiguration import assertCMakeHave, ConfigurationError
from dune.typeregistry import generateTypeName

class CartesianDomain(tuple):
    def __new__ (cls, lower,upper,division,**parameters):
        from ._grid import reader
        dgf = "DGF\n"
        dgf += "INTERVAL\n"
        dgf += " ".join([str(x) for x in lower]) + "\n"
        dgf += " ".join([str(x) for x in upper]) + "\n"
        dgf += " ".join([str(x) for x in division]) + "\n"
        dgf += "#\n"
        dgf += "GRIDPARAMETER\n"
        for key in parameters:
            dgf += key + " " + str(parameters[key]) + "\n"
        dgf += "#\n"
        try:
            periodic = parameters["periodic"]
            if any(periodic):
                dim = len(lower)
                # create identity matrix in comma separated rows
                eye = ""
                for i in range(dim):
                    for j in range(dim):
                        eye += " 1" if i == j else " 0"
                    eye += " + " if i == dim-1 else ", "

                dgf += "PERIODICFACETRANSFORMATION\n"
                for i,p in enumerate(periodic):
                    if p:
                        dgf += eye
                        dgf += " ".join([str(upper[i]-lower[i]) if i==j
                            else "0" for j in range(dim)])
                        dgf += "\n"
                dgf += "#\n"
        except KeyError:
            pass
        return super(CartesianDomain, cls).__new__(cls,
                       tuple( (reader.dgfString, dgf) ) )
    def __init__(self,lower,upper,division,**parameters):
        self.dimgrid = len(lower)
        self.lower = lower
        self.upper = upper
        self.division = division
        self.param = parameters
    def dimgrid(self):
        return self.dimgrid

def onedGrid(constructor):
    from .grid_generator import module, getDimgrid

    typeName = "Dune::OneDGrid"
    includes = ["dune/grid/onedgrid.hh", "dune/grid/io/file/dgfparser/dgfoned.hh"]
    gridModule = module(includes, typeName)

    return gridModule.reader(constructor).leafView

def moduleYaspCoordinates(dim, ctype="double"):
    moduleName = "yaspcoordinates_dim{dim}".format(ct = ctype, dim = dim)
    source = """
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/stl.h>

#include <dune/common/typelist.hh>
#include <dune/grid/yaspgrid/coordinates.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/common/fvector.hh>

using namespace Dune;
using namespace Dune::Python;
namespace py = pybind11;

template<typename T>
auto registerCoords(py::module module, std::string name)
{{
  auto includes = IncludeFiles{{"dune/grid/yaspgrid/coordinates.hh"}};
  std::string nspace("Dune::");
  auto typeName = GenerateTypeName(nspace+name, MetaType<{ct}>(), "{dim}");
  auto cls = insertClass<T>(module, name, typeName, includes);
  if (cls.second)
  {{
    // cls.first.def("name", [name](const T & self)
    cls.first.def_static("name", [name]() {{ return name; }});
    cls.first.def_property_readonly_static("typeName", [typeName](py::object) {{ return typeName.name(); }});
    cls.first.def_property_readonly_static("dimgrid", [](py::object) {{ return {dim}; }});
    cls.first.def_property_readonly_static("numpy_ctype", [](py::object) {{ return py::dtype::of<{ct}>(); }});
    cls.first.def_property_readonly_static("ctype", [](py::object) {{ return "{ct}"; }});
  }}
  return cls.first;
}}

PYBIND11_MODULE({moduleName}, module)
{{
  // make sure FieldVector is known to pybind11
  registerFieldVector<{ct},{dim}>(module);

  // EquidistantCoordinates(const Dune::FieldVector<ct,dim>& upperRight, const std::array<int,dim>& s)
  //py::class_<EquidistantCoordinates<{ct},{dim}>>(module, "EquidistantCoordinates")
  registerCoords<EquidistantCoordinates<{ct},{dim}>>(module, "EquidistantCoordinates")
    .def(py::init<Dune::FieldVector<{ct},{dim}>, std::array<int,{dim}>>());

  // EquidistantOffsetCoordinates(const Dune::FieldVector<ct,dim>& lowerLeft, const Dune::FieldVector<ct,dim>& upperRight, const std::array<int,dim>& s)
  registerCoords<EquidistantOffsetCoordinates<{ct},{dim}>>(module, "EquidistantOffsetCoordinates")
    .def(py::init<Dune::FieldVector<{ct},{dim}>, Dune::FieldVector<{ct},{dim}>, std::array<int,{dim}>>());
    //.def(py::init( [] (std::array<double,{dim}>, std::array<double,{dim}>, std::array<int,{dim}>) {{ return static_cast<EquidistantOffsetCoordinates<{ct},{dim}>*>(0); }} ));

  // TensorProductCoordinates(const std::array<std::vector<ct>,dim>& c, const std::array<int,dim>& offset)
  registerCoords<TensorProductCoordinates<{ct},{dim}>>(module, "TensorProductCoordinates")
    .def(py::init<std::array<std::vector<{ct}>,{dim}>, std::array<int,{dim}>>());
}}

""".format(moduleName = moduleName, ct = ctype, dim = dim)
    from dune.generator import builder
    module = builder.load(moduleName, source, "yasp coordinates dim={dim}".format(dim = dim)) # , self.typeName[0], extraCMake)
    return module

def equidistantOffsetCoordinates(lowerleft, upperright, elements, ctype='double'):
    import numpy as np
    dim = len(elements)
    mod = moduleYaspCoordinates(dim)
    coords_ = mod.EquidistantOffsetCoordinates
    # make sure we have float values
    dtype = coords_.numpy_ctype
    lowerleft = np.array(lowerleft, dtype=np.float64)
    upperright = np.array(upperright, dtype=np.float64)
    assert(len(lowerleft) == dim)
    assert(len(upperright) == dim)
    return coords_(lowerleft, upperright, elements)

def equidistantCoordinates(upperright, elements):
    import numpy as np
    dim = len(elements)
    # make sure we have float values
    upperright = np.array(upperright, dtype=np.float64)
    assert(len(upperright) == dim)
    mod = moduleYaspCoordinates(dim)
    coords_ = mod.EquidistantCoordinates
    return coords_(upperright, elements)

def tensorProductCoordinates(coords, offset=None, ctype='double'):
    dim = len(coords)
    mod = moduleYaspCoordinates(dim)
    coords_ = mod.TensorProductCoordinates
    if offset is None:
        offset = [0,]*dim
    if len(offset) != dim:
        raise ValueError("tensorProductCoordinates: offset parameter has wrong size")
    return coords_(coords,offset)

def yaspGrid(constructor, dimgrid=None, coordinates="equidistant", ctype="double", **param):
    """create a Dune::YaspGrid"""
    from dune.generator import Constructor
    from .grid_generator import module, getDimgrid

    if isinstance(constructor, CartesianDomain):
        constructor = equidistantOffsetCoordinates(
            lowerleft  = constructor.lower,
            upperright = constructor.upper,
            elements   = constructor.division,
            ctype      = ctype
        )

    if not dimgrid:
        dimgrid = getDimgrid(constructor)

    ctype = constructor.ctype
    coordinates_type = constructor.typeName
    typeName, _ = generateTypeName("Dune::YaspGrid", str(dimgrid), constructor)
    includes = ["dune/grid/yaspgrid.hh", "dune/grid/io/file/dgfparser/dgfyasp.hh"]

    ctor = Constructor(
            [ "const " + coordinates_type + "& coordinates",
              'std::array<bool, '+str(dimgrid)+'> periodic',
              'int overlap' ],
            [ 'std::bitset<'+str(dimgrid)+'> periodic_;',
              'for (int i=0;i<'+str(dimgrid)+';++i) periodic_.set(i,periodic[i]);',
              'return new DuneType(coordinates,periodic_,overlap);' ],
            [ '"coordinates"_a', '"periodic"_a', '"overlap"_a' ])

    gridModule = module(includes, typeName, ctor)
    # try creating a grid from the constructor
    try:
        print ("WARNING: add periodic and overlap parameters")
        periodic   = [False,]*dimgrid
        overlap    = 0
        return gridModule.HierarchicalGrid(constructor,periodic,overlap).leafView
    # if it fails, we assume that the constructor needs to be fed to a reader
    except AttributeError:
        return gridModule.reader(constructor).leafView

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
