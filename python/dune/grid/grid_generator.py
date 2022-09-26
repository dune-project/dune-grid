# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import os, inspect
from ..generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dune.common import FieldVector
from dune.common.utility import isString
from dune.deprecate import deprecated
from dune.grid import gridFunction, DataType
from dune.grid import OutputType
from dune.generator.algorithm import cppType
from dune.generator import builder
from dune.deprecate import deprecated

def getDimgrid(constructor):
    dimgrid = None
    if not dimgrid:
        try:
            dimgrid = constructor.dimgrid
        except AttributeError:
            pass
    if not dimgrid:
        try:
            dimgrid = len(constructor["vertices"][0])
        except KeyError:
            pass
    if not dimgrid:
        raise ValueError("Couldn't extract dimension of grid from constructor arguments, added dimgrid parameter")
    return dimgrid

def triangulation(grid, level=0):
    if grid.dimGrid != 2:
        raise Exception("Grid must be 2-dimensional for use as matplotlib triangulation.")
    from matplotlib.tri import Triangulation
    x, triangles = grid.tessellate(level)
    return Triangulation(x[:,0], x[:,1], triangles)

_writeVTKDispatcher = []
def _writeVTK(vtk,grid,f,name,dataTag):
    done = False
    try:
        f.addToVTKWriter(name, vtk, dataTag)
        done = True
    except AttributeError:
        pass
    if not done:
        for dispatch in _writeVTKDispatcher:
            try:
                func = dispatch(grid,f)
            except:
                func = None
            if func is not None:
                func.addToVTKWriter(name,vtk,dataTag)
                done = True
                break
    if not done:
        gridFunction(grid)(f).addToVTKWriter(name, vtk, dataTag)

def writeVTK(grid, name,
             celldata=None, pointdata=None,
             cellvector=None, pointvector=None,
             number=None, subsampling=None,outputType=OutputType.appendedraw,
             write=True, nonconforming=False):
    vtk = grid.vtkWriter(nonconforming) if subsampling is None else grid.vtkWriter(subsampling)

    def addDataToVTKWriter(dataFunctions, dataName, dataTag):
        if dataFunctions is None: return
        if isinstance(dataFunctions, dict):
            for n, f in dataFunctions.items():
                if f is None: continue
                _writeVTK(vtk,grid,f,n,dataTag)
        elif isinstance(dataFunctions, list):
            for f in dataFunctions:
                if f is None: continue
                try:
                    _writeVTK(vtk,grid,f,f.name,dataTag)
                except AttributeError:
                    try:
                        _writeVTK(vtk,grid,f[0],f[1],dataTag)
                    except IndexError:
                        raise TypeError("""
Did you try to pass in a function without a name attribute?
Try using a dictionary with name:function instead.""")

        elif dataFunctions is not None:
            raise TypeError("Argument '" + dataName + "' must be a dict or list instance.")

    addDataToVTKWriter(celldata, 'celldata', DataType.CellData)
    addDataToVTKWriter(pointdata, 'pointdata', DataType.PointData)
    addDataToVTKWriter(cellvector, 'cellvector', DataType.CellVector)
    addDataToVTKWriter(pointvector, 'pointvector', DataType.PointVector)

    assert isinstance(outputType,OutputType)
    if write:
        if number is None:
            vtk.write(name, outputType)
        else:
            vtk.write(name, number, outputType)
    else:
        return vtk

class SequencedVTK:
    def __init__(self, grid, name, number,
                 celldata, pointdata, cellvector, pointvector,
                 subsampling, outputType=OutputType.appendedraw):
        self.number = number
        self.name = name
        self.vtk = grid.writeVTK(name,celldata=celldata,pointdata=pointdata,cellvector=cellvector,pointvector=pointvector,subsampling=subsampling,write=False)
        self.outputType = outputType
    def __call__(self):
        self.vtk.write(self.name, self.number, self.outputType)
        self.number += 1

def sequencedVTK(grid, name, celldata=None, pointdata=None, cellvector=None, pointvector=None,
                 number=0, subsampling=None, outputType=OutputType.appendedraw):
    return SequencedVTK(grid,name,number,
                        celldata=celldata,pointdata=pointdata,
                        cellvector=cellvector,pointvector=pointvector,
                        subsampling=subsampling,outputType=outputType)

def plot(self, function=None, *args, **kwargs):
    import dune.plotting
    if not function:
        dune.plotting.plotGrid(self, *args, **kwargs)
    else:
        if not hasattr(function,"grid"):
            function = self.function(function)
        dune.plotting.plot(solution=function,*args,**kwargs)

isGenerator = SimpleGenerator("GridViewIndexSet", "Dune::Python")
def indexSet(gv):
    try:
        return gv._indexSet
    except TypeError:
        includes = gv.cppIncludes + ["dune/python/grid/indexset.hh"]
        typeName = gv.cppTypeName+"::IndexSet"
        moduleName = "indexset_" + hashIt(typeName)
        module = isGenerator.load(includes, typeName, moduleName)
        return gv._indexSet
mcmgGenerator = SimpleGenerator("MultipleCodimMultipleGeomTypeMapper", "Dune::Python")
def mapper(gv,layout):
    includes = gv.cppIncludes + ["dune/python/grid/mapper.hh"]
    typeName = "Dune::MultipleCodimMultipleGeomTypeMapper< "+gv.cppTypeName+" >"
    moduleName = "mcmgmapper_" + hashIt(typeName)
    module = mcmgGenerator.load(includes, typeName, moduleName)
    return gv._mapper(layout)
def function(gv,callback,includeFiles=None,*args,name=None,order=None,dimRange=None):
    if name is None:
        name = "tmp"+str(gv._gfCounter)
        gv.__class__._gfCounter += 1
    if isString(callback):
        if includeFiles is None:
            raise ValueError("""if `callback` is the name of a C++ function
            then at least one include file containing that function must be
            provided""")

        # unique header guard is added further down
        source  = '#include <config.h>\n\n'
        source += '#define USING_DUNE_PYTHON 1\n\n'
        includes = []
        if isString(includeFiles):
            if not os.path.dirname(includeFiles):
                with open(includeFiles, "r") as include:
                    source += include.read()
                source += "\n"
            else:
                source += "#include <"+includeFiles+">\n"
                includes += [includeFiles]
        elif hasattr(includeFiles,"readable"): # for IOString
            with includeFiles as include:
                source += include.read()
            source += "\n"
        elif isinstance(includeFiles, list):
            for includefile in includeFiles:
                if not os.path.dirname(includefile):
                    with open(includefile, "r") as include:
                        source += include.read()
                    source += "\n"
            else:
                source += "#include <"+includefile+">\n"
                includes += [includefile]
        includes += gv.cppIncludes
        argTypes = []
        for arg in args:
            t,i = cppType(arg)
            argTypes.append(t)
            includes += i

        signature = callback + "( " + ", ".join(argTypes) + " )"
        moduleName = "gf_" + hashIt(signature) + "_" + hashIt(source)

        # add unique header guard with moduleName
        source = '#ifndef Guard_'+moduleName+'\n' + \
                 '#define Guard_'+moduleName+'\n\n' + \
                 source

        includes = sorted(set(includes))
        source += "".join(["#include <" + i + ">\n" for i in includes])
        source += "\n"
        source += '#include <dune/python/grid/function.hh>\n'
        source += '#include <dune/python/pybind11/pybind11.h>\n'
        source += '\n'

        source += "PYBIND11_MODULE( " + moduleName + ", module )\n"
        source += "{\n"
        source += "  module.def( \"gf\", [module] ( "+gv.cppTypeName + " &gv"+"".join([", "+argTypes[i] + " arg" + str(i) for i in range(len(argTypes))]) + " ) {\n"
        source += "      auto callback="+callback+"<"+gv.cppTypeName+">( "+",".join(["arg"+str(i) for i in range(len(argTypes))]) +"); \n"
        source += "      return Dune::Python::registerGridFunction<"+gv.cppTypeName+",decltype(callback)>(module,pybind11::cast(gv),\"tmp\",callback);\n"
        source += "    },"
        source += "    "+",".join(["pybind11::keep_alive<0,"+str(i+1)+">()" for i in range(len(argTypes)+1)])
        source += ");\n"
        source += "}\n"
        source += "#endif\n"
        gf = builder.load(moduleName, source, signature).gf(gv,*args)
    else:
        if len(inspect.signature(callback).parameters) == 1: # global function, turn into a local function
            callback_ = callback
            callback = lambda e,x: callback_(e.geometry.toGlobal(x))
        else:
            callback_ = None
        if dimRange is None:
            # if no `dimRange` attribute is set on the callback,
            # try to evaluate the function to determin the dimension of
            # the return value. This can fail if the function is singular in
            # the computational domain in which case an exception is raised
            e = gv.elements.__iter__().__next__()
            try:
                y = callback(e,e.referenceElement.position(0,0))
            except ArithmeticError:
                try:
                    y = callback(e,e.referenceElement.position(0,2))
                except ArithmeticError:
                    raise TypeError("can not determin dimension of range of "+
                      "given grid function due to arithmetic exceptions being "+
                      "raised. Add a `dimRange` parameter to the grid function to "+
                      "solve this issue - set `dimRange`=0 for a scalar function.")
            try:
                dimRange = len(y)
            except TypeError:
                dimRange = 0
        if dimRange > 0:
            scalar = "false"
        else:
            scalar = "true"
        FieldVector(dimRange*[0]) # register FieldVector for the return value
        if not dimRange in gv.__class__._functions.keys():
            # unique header key is added further down
            source  = '#include <config.h>\n\n'
            source += '#define USING_DUNE_PYTHON 1\n\n'
            includes = gv.cppIncludes

            signature = gv.cppTypeName+"::gf<"+str(dimRange)+">"
            moduleName = "gf_" + hashIt(signature) + "_" + hashIt(source)

            # add unique header guard with moduleName
            source = '#ifndef Guard_'+moduleName+'\n' + \
                     '#define Guard_'+moduleName+'\n\n' + \
                     source

            includes = sorted(set(includes))
            source += "".join(["#include <" + i + ">\n" for i in includes])
            source += "\n"
            source += '#include <dune/python/grid/function.hh>\n'
            source += '#include <dune/python/pybind11/pybind11.h>\n'
            source += '\n'

            source += "PYBIND11_MODULE( " + moduleName + ", module )\n"
            source += "{\n"
            source += "  typedef pybind11::function Evaluate;\n";
            source += "  Dune::Python::registerGridFunction< "+gv.cppTypeName+", Evaluate, "+str(dimRange)+" >( module, \"gf\", "+scalar+" );\n"
            source += "}\n"
            source += "#endif\n"
            gfModule = builder.load(moduleName, source, signature)
            gfFunc = getattr(gfModule,"gf"+str(dimRange))
            if callback_ is not None:
                gfFunc.localCall = gfFunc.__call__
                feval = lambda self,e,x=None: callback_(e) if x is None else self.localCall(e,x)
                subclass = type(gfFunc.__name__, (gfFunc,), {"__call__": feval})
                gv.__class__._functions[dimRange] = subclass
            else:
                gv.__class__._functions[dimRange] = gfFunc
        gf = gv.__class__._functions[dimRange](gv,callback)
    def gfPlot(gf, *args, **kwargs):
        gf.grid.plot(gf,*args,**kwargs)
    gf.plot = gfPlot.__get__(gf)
    gf.name = name
    gf.order = order
    return gf

def addAttr(module, cls):
    setattr(cls, "_module", module)
    setattr(cls, "writeVTK", writeVTK)
    setattr(cls, "sequencedVTK", sequencedVTK)
    setattr(cls, "_functions", {})
    @deprecated(name="dune.grid.GridView.tesselate", msg="Use 'tessellate' (note spelling)")
    def tesselate(gv, *args,**kwargs):
        return gv.tessellate(*args,**kwargs)
    setattr(cls,"tesselate", tesselate)
    [[deprecated("use 'tessellate' (note spelling)")]]

    if cls.dimension == 2:
        setattr(cls, "plot", plot)
        setattr(cls, "triangulation", triangulation)
    else:
        def throwFunc(msg):
            def throw(*args, **kwargs):
                raise AttributeError(msg)
            return throw
        setattr(cls, "plot", throwFunc("plot(...) only implemented on 2D grids"))
        setattr(cls, "triangulation", throwFunc("triangulation(...) only implemented on 2d grid"))

    cls.indexSet = property(indexSet)
    setattr(cls,"mapper",mapper)
    setattr(cls,"function",function)
    setattr(cls,"_gfCounter",0)
def addHAttr(module):
    # register reference element for this grid
    import dune.geometry
    for d in range(module.LeafGrid.dimension+1):
        dune.geometry.module(d)
    setattr(module.HierarchicalGrid,"levelView",levelView)
    setattr(module.HierarchicalGrid,"persistentContainer",persistentContainer)
    # setattr(module.HierarchicalGrid,"backup",backup)
    addAttr(module, module.LeafGrid)

gvGenerator = SimpleGenerator("GridView", "Dune::Python")
def viewModule(includes, typeName, *args, **kwargs):
    includes = includes + ["dune/python/grid/gridview.hh"]
    moduleName = "view_" + hashIt(typeName)
    module = gvGenerator.load(includes, typeName, moduleName, *args, **kwargs)
    addAttr(module, module.GridView)
    return module

def levelView(hgrid,level):
    includes = hgrid.cppIncludes
    typeName = "typename "+hgrid.cppTypeName+"::LevelGridView"
    viewModule(includes, typeName)
    return hgrid._levelView(level)

pcGenerator = SimpleGenerator("PersistentContainer", "Dune::Python")
def persistentContainer(hgrid,codim,dimension):
    includes = hgrid.cppIncludes + ["dune/python/grid/persistentcontainer.hh"]
    typeName = "Dune::PersistentContainer<"+hgrid.cppTypeName+", Dune::FieldVector<double,"+str(dimension)+">>"
    moduleName = "persistentcontainer_" + hashIt(typeName)
    module = pcGenerator.load(includes, typeName, moduleName)
    return module.PersistentContainer(hgrid,codim)

def module(includes, typeName, *args, **kwargs):
    try:
        generator = kwargs.pop("generator")
    except KeyError:
        generator = SimpleGenerator("HierarchicalGrid", "Dune::Python")
    includes = includes + ["dune/python/grid/hierarchical.hh"]
    typeHash = "hierarchicalgrid_" + hashIt(typeName)
    kwargs["dynamicAttr"] = True
    kwargs["holder"] = "std::shared_ptr"
    module = generator.load(includes, typeName, typeHash, *args, **kwargs)
    return module

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
