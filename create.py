from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import sys
import importlib
from types import ModuleType

import dune.common as common
from ..generator.generator import SimpleGenerator

from types import ModuleType

def triangulation(self):
    if self.dimGrid != 2 or self.dimWorld != 2:
        raise Exception("Grid must be 2-dimensional for use as matplotlib triangulation.")
    from matplotlib.tri import Triangulation
    x = self.coordinates()
    triangles = self.tesselate()
    return Triangulation(x[:,0], x[:,1], triangles)

def addAttr(module, cls):
    setattr(cls, "_module", module)
    setattr(cls, "triangulation", triangulation)

generator = SimpleGenerator("Grid", "Dune::CorePy", "LeafGrid")

def module(includes, typeName, constructors=None, methods=None):
    typeName = typeName + "::LeafGridView"
    includes = includes + ["dune/corepy/grid.hh"]
    typeHash = "grid_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, typeHash, constructors, methods)
    addAttr(module, module.LeafGrid)
    return module

gridNames = { "Alberta"        : "dune.grid.alberta",
              "OneDGrid"       : "dune.grid.aoned",
              "SPGrid"         : "dune.grid.asp",
              "UGGrid"         : "dune.grid.aug",
              "YaspGrid"       : "dune.grid.ayasp"
            }
def create(grid, *args, **kwargs):
    gridModule = importlib.import_module(gridNames[grid])
    return gridModule.create(*args,**kwargs)

#############################################
if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
