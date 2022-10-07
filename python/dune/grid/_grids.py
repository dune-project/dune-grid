# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

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
        foundRefEdge = False
        for key in parameters:
            if key.lower() == 'refinementedge':
                foundRefEdge = True
            dgf += key + " " + str(parameters[key]) + "\n"
        # set default value for refinementedge if not provided to avoid warning
        if not foundRefEdge:
            dgf += "REFINEMENTEDGE ARBITRARY\n"
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
    moduleName = "yaspcoordinates_dim{dim}_ct{ct}".format(ct = ctype, dim = dim)
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
  addToTypeRegistry<{ct}>(GenerateTypeName("{ct}"));
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
    module = builder.load(moduleName, source, "yasp coordinates dim={dim} ctype={ct}".format(ct = ctype, dim = dim))
    return module

def equidistantOffsetCoordinates(lowerleft, upperright, elements, ctype='double'):
    import numpy as np
    dim = len(elements)
    mod = moduleYaspCoordinates(dim,ctype)
    coords_ = mod.EquidistantOffsetCoordinates
    # make sure we have float values
    dtype = coords_.numpy_ctype
    lowerleft = np.array(lowerleft, dtype=dtype)
    upperright = np.array(upperright, dtype=dtype)
    assert(len(lowerleft) == dim)
    assert(len(upperright) == dim)
    return coords_(lowerleft, upperright, elements)

def equidistantCoordinates(upperright, elements):
    import numpy as np
    dim = len(elements)
    mod = moduleYaspCoordinates(dim,ctype)
    coords_ = mod.EquidistantCoordinates
    # make sure we have float values
    dtype = coords_.numpy_ctype
    upperright = np.array(upperright, dtype=dtype)
    assert(len(upperright) == dim)
    return coords_(upperright, elements)

def tensorProductCoordinates(coords, offset=None, ctype='double'):
    import numpy as np
    dim = len(coords)
    mod = moduleYaspCoordinates(dim,ctype)
    coords_ = mod.TensorProductCoordinates
    if offset is None:
        offset = [0,]*dim
    if len(offset) != dim:
        raise ValueError("tensorProductCoordinates: offset parameter has wrong size")
    dtype = coords_.numpy_ctype
    coords = np.array(coords, dtype=dtype)
    return coords_(coords,offset)

def yaspGrid(constructor, dimgrid=None, coordinates="equidistant", ctype=None,
             periodic=None, overlap=None, **param):
    """create a Dune::YaspGrid

       constructor: a Yaspgrod coordinates object
                    (created via equidistantCoordinates,
                                 equidistantOffsetCoordinates or
                                 tensorProductCoordinates)
                    or CartesianDomain
                    or a reader object that can be called via gridModule.reader(constructor)
       dimgrid:     explicitly specify the dimension of the grid (should only be set if using a reader)
       coordinates: explicitly specify coordinates template parameter (deprecated, and ignored)
       ctype:       explicitly specify grids C++ ctype template parameter (should only be set if using a reader, default: double)
       periodic:    boolean list specifying periodicity per dimension (default: no periodicity)
       overlap:     size of overlap for periodic or parallel grids (default: 1)


    parameters coordinates and ctype are deprecated and will in future
    be deduced from the constructor parameter.

    """
    from dune.generator import Constructor
    from .grid_generator import module, getDimgrid

    # we have to keep track whether we create the grid via a constructor call or via a reader
    useReader = False

    #### we try to guess what kind of "constructor" we are using now
    # CartesianDomain
    if isinstance(constructor, CartesianDomain):
        # retrieve parameters from constructor.param, explicit parameters take precedence
        ctype_     = constructor.param.get("ctype", ctype)
        if not ctype:
            ctype = ctype_
        else:
            if ctype != ctype_:
                print("WARNING: yaspGrid: ctype Parameter of CartesianDomain overwritten by explicit parameter")
        if not ctype:
            ctype = "double" # default is doubel
        # ---
        periodic_  = constructor.param.get("periodic", periodic)
        if not periodic:
            periodic = periodic_ # if periodic is not gieven on the command line, we use the parameter from CartesianDomain
        else:
            if periodic != periodic_:
                print("WARNING: yaspGrid: periodic Parameter of CartesianDomain overwritten by explicit parameter")
        # ---
        overlap_   = constructor.param.get("overlap", overlap)
        if not overlap:
            overlap = overlap_ # if overlap is not gieven on the command line, we use the parameter from CartesianDomain
        else:
            if overlap != overlap_:
                print("WARNING: yaspGrid: overlap Parameter of CartesianDomain overwritten by explicit parameter")
        # ---
        constructor = equidistantOffsetCoordinates(
            lowerleft  = constructor.lower,
            upperright = constructor.upper,
            elements   = constructor.division,
            ctype      = ctype
        )

    from dune.grid import reader
    # Coordinate object
    if (not isinstance(constructor, tuple)):
        dimgrid_    = getDimgrid(constructor)
        if dimgrid:
            if dimgrid != dimgrid_:
                raise ValueError("yaspGrid: parameter dimgrid can only be specified, when using a reader")
        else:
            dimgrid = dimgrid_
        if ctype:
            if ctype != constructor.ctype:
                raise ValueError("yaspGrid: ctype is specified via coordinate object and can not be overwritten")
        ctype = constructor.ctype
        coordinates_type = constructor.typeName
        dimgrid    = getDimgrid(constructor)
        if periodic is None:
            periodic   = [False,]*dimgrid
        if overlap is None:
            overlap   = 1
    # Reader
    elif isinstance(constructor, tuple) and isinstance(constructor[0], reader):
        if not ctype:
            ctype = "double"
        if not dimgrid:
            raise ValueError("yaspGrid: parameter dimgrid must be specified, when using a reader")
        mod = moduleYaspCoordinates(dimgrid,ctype)
        # reader only works with EquidistantOffsetCoordinates or EquidistantCoordinates
        coords = mod.EquidistantOffsetCoordinates
        coordinates_type = coords.typeName
        useReader = True
        # periodic and overlap are (if required) defined by the DGF reader
        periodic   = None
        overlap   = None
    else:
        raise ValueError("yaspGrid: unsupported constructor parameter " + str(constructor))

    # compile YaspGrid for a given dimension & coordinate type
    includes = ["dune/grid/yaspgrid.hh", "dune/grid/io/file/dgfparser/dgfyasp.hh"]
    gridTypeName, _ = generateTypeName("Dune::YaspGrid", str(dimgrid), coordinates_type)
    ctor = Constructor(
            [ "const " + coordinates_type + "& coordinates",
              'std::array<bool, '+str(dimgrid)+'> periodic',
              'int overlap' ],
            [ 'std::bitset<'+str(dimgrid)+'> periodic_;',
              'for (int i=0;i<'+str(dimgrid)+';++i) periodic_.set(i,periodic[i]);',
              'return new DuneType(coordinates,periodic_,overlap);' ],
            [ '"coordinates"_a', '"periodic"_a', '"overlap"_a' ])
    gridModule = module(includes, gridTypeName, ctor)

    # read the grid either via reader or create it directly...
    if useReader:
        return gridModule.reader(constructor).leafView
    else:
        return gridModule.HierarchicalGrid(constructor,periodic,overlap).leafView

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
