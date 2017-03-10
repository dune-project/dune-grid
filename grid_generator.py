from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import dune.common as common
from ..generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

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
                f.addToVTKWriter(n, vtk, dataTag)
        elif isinstance(dataFunctions, list):
            for f in dataFunctions:
                f.addToVTKWriter(f.name, vtk, dataTag)
        elif dataFunctions is not None:
            raise TypeError("Argument '" + dataName + "' must be a dict instance.")

    addDataToVTKWriter(celldata, 'celldata', common.DataType.CellData)
    addDataToVTKWriter(pointdata, 'pointdata', common.DataType.PointData)
    addDataToVTKWriter(cellvector, 'cellvector', common.DataType.CellVector)
    addDataToVTKWriter(pointvector, 'pointvector', common.DataType.PointVector)

    if number is None:
        vtk.write(name)
    else:
        vtk.write(name, number)


def addAttr(module, cls):
    setattr(cls, "_module", module)
    setattr(cls, "triangulation", triangulation)
    setattr(cls, "writeVTK", writeVTK)


generator = SimpleGenerator("HierarchicalGrid", "Dune::CorePy")


def module(includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/corepy/grid/hierarchical.hh"]
    typeHash = "hierarchicalgrid_" + hashIt(typeName)
    module = generator.load(includes, typeName, typeHash, constructors, methods)
    addAttr(module, module.HierarchicalGrid.LeafView)
    addAttr(module, module.HierarchicalGrid.LevelView)
    return module


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
