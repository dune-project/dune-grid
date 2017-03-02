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


def writeVTK(grid, name, celldata=[], pointdata=[], cellvector=[], pointvector=[], number=None, subsampling=None):
    vtk = grid.vtkWriter() if subsampling is None else grid.vtkWriter(subsampling)

    for f in celldata:
        f.addToVTKWriter(vtk, common.DataType.CellData)
    for f in pointdata:
        f.addToVTKWriter(vtk, common.DataType.PointData)
    for f in cellvector:
        f.addToVTKWriter(vtk, vtk.CellVectorData)
    for f in pointvector:
        f.addToVTKWriter(vtk, vtk.PointVectorData)

    if number is None:
        vtk.write(name)
    else:
        vtk.write(name, number)
    return vtk


def addAttr(module, cls):
    setattr(cls, "_module", module)
    setattr(cls, "triangulation", triangulation)
    setattr(cls, "writeVTK", writeVTK)


generator = SimpleGenerator("Grid", "Dune::CorePy", "LeafGrid")
fileBase = "grid"


def module(includes, typeName, constructors=None, methods=None):
    typeName = typeName + "::LeafGridView"
    includes = includes + ["dune/corepy/grid.hh"]
    typeHash = fileBase + "_" + hashIt(typeName)
    module = generator.load(includes, typeName, typeHash, constructors, methods)
    addAttr(module, module.LeafGrid)
    return module


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
