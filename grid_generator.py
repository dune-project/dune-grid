from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

from ..generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
from dune.common import _raise
from dune.deprecate import deprecated
from dune.grid import gridFunction, DataType

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
            raise ValueError("Couldn't extract dimension of grid from constructor arguments, added dimgrid parameter")
    return dimgrid

def triangulation(grid, level=0):
    if grid.dimGrid != 2:
        raise Exception("Grid must be 2-dimensional for use as matplotlib triangulation.")
    from matplotlib.tri import Triangulation
    x, triangles = grid.tesselate(level)
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
            func = dispatch(grid,f)
            if func is not None:
                func.addToVTKWriter(name,vtk,dataTag)
                done = True
                break
    if not done:
        @gridFunction(grid)
        def f_(*args,**kwargs):
            return f(*args,**kwargs)
        f_.addToVTKWriter(name, vtk, dataTag)

def writeVTK(grid, name, celldata=None, pointdata=None, cellvector=None, pointvector=None, number=None, subsampling=None, write=True):
    vtk = grid.vtkWriter() if subsampling is None else grid.vtkWriter(subsampling)

    def addDataToVTKWriter(dataFunctions, dataName, dataTag):
        if isinstance(dataFunctions, dict):
            for n, f in dataFunctions.items():
                _writeVTK(vtk,grid,f,n,dataTag)
        elif isinstance(dataFunctions, list):
            for f in dataFunctions:
                try:
                    _writeVTK(vtk,grid,f,f.name,dataTag)
                except AttributeError:
                    _writeVTK(vtk,grid,f[0],f[1],dataTag)
        elif dataFunctions is not None:
            raise TypeError("Argument '" + dataName + "' must be a dict or list instance.")

    addDataToVTKWriter(celldata, 'celldata', DataType.CellData)
    addDataToVTKWriter(pointdata, 'pointdata', DataType.PointData)
    addDataToVTKWriter(cellvector, 'cellvector', DataType.CellVector)
    addDataToVTKWriter(pointvector, 'pointvector', DataType.PointVector)

    if write:
        if number is None:
            vtk.write(name)
        else:
            vtk.write(name, number)
    else:
        return vtk

class SequencedVTK:
    def __init__(self, grid, name, number, celldata, pointdata, cellvector, pointvector, subsampling):
        self.number = number
        self.name = name
        self.vtk = grid.writeVTK(name,celldata=celldata,pointdata=pointdata,cellvector=cellvector,pointvector=pointvector,subsampling=subsampling,write=False)
    def __call__(self):
        self.vtk.write(self.name, self.number)
        self.number += 1

def sequencedVTK(grid, name, celldata=None, pointdata=None, cellvector=None, pointvector=None, number=0, subsampling=None):
    return SequencedVTK(grid,name,number,celldata=celldata,pointdata=pointdata,cellvector=cellvector,pointvector=pointvector,subsampling=subsampling)

def plot(self, function=None, *args, **kwargs):
    import dune.plotting
    if not function:
        dune.plotting.plotGrid(self, *args, **kwargs)
    else:
        try:
            grid = function.grid
            dune.plotting.plot(solution=function,*args,**kwargs)
        except AttributeError:
            dune.plotting.plot(solution=self.function(function),*args,**kwargs)

@deprecated("use the `gridFunction` decorator")
def globalGridFunction(gv, evaluator):
    return gv.function(evaluator)
@deprecated("use the `gridFunction` decorator")
def localGridFunction(gv, evaluator):
    return gv.function( lambda x: evaluator(x.entity,x.local) )
@deprecated("use the `referenceElement` attribute instead")
def domain(self):
    return self.referenceElement
@deprecated("use the `toGlobal` attribute instead")
def position(self,*arg,**kwarg):
    return self.toGlobal(*arg,**kwarg)
@deprecated("use the `toLocal` attribute instead")
def localPosition(self,*arg,**kwarg):
    return self.toLocal(*arg,**kwarg)

def addAttr(module, cls):
    setattr(cls, "_module", module)
    setattr(cls, "writeVTK", writeVTK)
    setattr(cls, "sequencedVTK", sequencedVTK)

    if cls.dimension == 2:
        setattr(cls, "plot", plot)
        setattr(cls, "triangulation", triangulation)
    else:
        setattr(cls, "plot", lambda *arg,**kwarg: _raise(AttributeError("plot only implemented on 2D grids")))
        setattr(cls, "triangulation", lambda *arg,**kwarg: _raise(AttributeError("triangulation only implemented on 2d grid")))

    setattr(cls, "globalGridFunction", globalGridFunction)
    setattr(cls, "localGridFunction", localGridFunction)

    def gfPlot(gf, *args, **kwargs):
        gf.grid.plot(gf,*args,**kwargs)
    for gf in dir(cls):
        if gf.startswith("GridFunction"):
            setattr( getattr(cls, gf), "plot", gfPlot)
    for ent in dir(cls):
        if ent.startswith("Entity"):
            Ent = getattr(cls, ent)
            Ent.domain = property(domain)
            Geo = getattr(Ent, "Geometry")
            Geo.domain = property(domain)
            setattr( Geo, "position", position)
            setattr( Geo, "localPosition", localPosition)


generator = SimpleGenerator("HierarchicalGrid", "Dune::Python")


def module(includes, typeName, *args):
    includes = includes + ["dune/python/grid/hierarchical.hh"]
    typeHash = "hierarchicalgrid_" + hashIt(typeName)
    module = generator.load(includes, typeName, typeHash, *args)
    addAttr(module, module.LeafGrid)
    addAttr(module, module.LevelGrid)

    # register reference element for this grid
    import dune.geometry
    for d in range(module.LeafGrid.dimension+1):
        dune.geometry.module(d)

    return module


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
