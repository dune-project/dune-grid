from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import dune.common as common
from ..generator.generator import SimpleGenerator
from dune.common.hashit import hashIt
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
            raise ValueError("Couldn't extract dimension of grid from constructor arguments, added dimgrid parameter")
    return dimgrid

def triangulation(grid, level=0):
    if grid.dimGrid != 2:
        raise Exception("Grid must be 2-dimensional for use as matplotlib triangulation.")
    from matplotlib.tri import Triangulation
    x, triangles = grid.tesselate(level)
    return Triangulation(x[:,0], x[:,1], triangles)


def writeVTK(grid, name, celldata=None, pointdata=None, cellvector=None, pointvector=None, number=None, subsampling=None):
    vtk = grid.vtkWriter() if subsampling is None else grid.vtkWriter(subsampling)

    def addDataToVTKWriter(dataFunctions, dataName, dataTag):
        if isinstance(dataFunctions, dict):
            for n, f in dataFunctions.items():
                try:
                    f.addToVTKWriter(n, vtk, dataTag)
                except AttributeError:
                    gf = grid.function(f)
                    gf.addToVTKWriter(n, vtk, dataTag)
        elif isinstance(dataFunctions, list):
            for f in dataFunctions:
                f.addToVTKWriter(f.name, vtk, dataTag)
        elif dataFunctions is not None:
            raise TypeError("Argument '" + dataName + "' must be a dict or list instance.")

    addDataToVTKWriter(celldata, 'celldata', common.DataType.CellData)
    addDataToVTKWriter(pointdata, 'pointdata', common.DataType.PointData)
    addDataToVTKWriter(cellvector, 'cellvector', common.DataType.CellVector)
    addDataToVTKWriter(pointvector, 'pointvector', common.DataType.PointVector)

    if number is None:
        vtk.write(name)
    else:
        vtk.write(name, number)
    return vtk

def plot(self, function=None, *args, **kwargs):
    import dune.plotting
    if not function:
        dune.plotting.plotGrid(self, *args, **kwargs)
    else:
        try:
            grid = function.grid
            dune.plotting.plotPointData(solution=function,*args,**kwargs)
        except AttributeError:
            dune.plotting.plotPointData(solution=self.function(function),*args,**kwargs)

def mapper(self,layout):
    from dune.grid.map import MultipleCodimMultipleGeomTypeMapper as Mapper
    return Mapper(self,layout)

@deprecated
def globalGridFunction(gv, evaluator):
    return gv.function(evaluator)
@deprecated
def localGridFunction(gv, evaluator):
    return gv.function( lambda x: evaluator(x.entity,x.local) )

def addAttr(module, cls):
    setattr(cls, "_module", module)
    setattr(cls, "writeVTK", writeVTK)
    setattr(cls, "mapper", mapper)

    if cls.dimension == 2:
        setattr(cls, "plot", plot)
        setattr(cls, "triangulation", triangulation)
    else:
        setattr(cls, "plot", lambda *arg,**kwarg: common._raise(AttributeError("plot only implemented on 2D grids")))
        setattr(cls, "triangulation", lambda *arg,**kwarg: common._raise(AttributeError("triangulation only implemented on 2d grid")))

    setattr(cls, "globalGridFunction", globalGridFunction)
    setattr(cls, "localGridFunction", localGridFunction)


generator = SimpleGenerator("HierarchicalGrid", "Dune::Python")


def module(includes, typeName, *args):
    includes = includes + ["dune/python/grid/hierarchical.hh"]
    typeHash = "hierarchicalgrid_" + hashIt(typeName)
    module = generator.load(includes, typeName, typeHash, *args)
    addAttr(module, module.LeafGrid)
    addAttr(module, module.LevelGrid)
    return module


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
