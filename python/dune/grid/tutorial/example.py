# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import math
import sys, os

# find grid files relative to example.py script
griddir = os.path.dirname(sys.argv[0])

# example of how to perform operations on a given grid
def runOnGrid(grid):

    # define a grid function and visualize
    from dune.grid import gridFunction
    @gridFunction(grid)
    def f(x):
        return math.cos(2.*math.pi/(1+abs(x[0]*x[1])))
    f.plot()
    grid.writeVTK("example", pointdata={"gf":f})

    # integrate the grid function using a quadrature rule from dune.geometry
    from dune.geometry import quadratureRules
    rules = quadratureRules(5)
    l2norm2 = 0
    for e in grid.elements:
        geo = e.geometry
        for qp in rules(e.type):
            x,w = qp.position, qp.weight
            l2norm2 += f(e,x)**2*w*geo.integrationElement(x)
    print("integral of grid function=",math.sqrt(l2norm2))

# construct ugGrid and yaspGrid via file reader
from dune.grid import ugGrid, reader
print ("construct an unstructured Grid (ugGrid) via file reader")
mshfile = os.path.join(griddir, "circle1storder.msh")
unstructuredGrid = ugGrid( (reader.gmsh, mshfile), dimgrid=2 )
if not unstructuredGrid:
    print ("WARNING: skipped ugGrid example, as dune-uggrid is not installed")
else:
    unstructuredGrid.plot()
    runOnGrid(unstructuredGrid)

print ("construct a Grid via file reader")
from dune.grid import yaspGrid, reader
mshfile = os.path.join(griddir, "test2d_offset.dgf")
print(mshfile)
dgfgrid = yaspGrid( (reader.dgf, mshfile), dimgrid=2 )
dgfgrid.plot()
runOnGrid(dgfgrid)

# construct a Cartesian grid
print ("construct a Cartesian grid")
from dune.grid import structuredGrid
grid = structuredGrid([-1,-1],[1,1],[10,10])
print("number of elements of Cartesian grid:",grid.size(0))
grid.plot()
runOnGrid(grid)

# construct YaspGrids with different coordinate types
# yaspGrid allows to specify the coordinate type
print ("construct a YaspGrid with tensor product coordinate type")
from dune.grid import yaspGrid, tensorProductCoordinates
import numpy as np
coords = tensorProductCoordinates([np.array([1,2,3,4]), np.array([10,11,33,44])], ctype='float')
ygrid = yaspGrid(coords)
print("number of elements of tensor YaspGrid grid:",ygrid.size(0))
ygrid.plot()
runOnGrid(ygrid)

# create a YaspGrid for a cartesian domain with non-standard overlap and periodicity
print ("construct a YaspGrid as a CartesianDomain with non-standard overlap and periodicity")
from dune.grid import yaspGrid, cartesianDomain
dim = 2
cartDomain = cartesianDomain([-10]*dim,[10]*dim, [20]*dim, periodic=[True]*dim, overlap=2)
p_grid = yaspGrid( cartDomain, dimgrid=dim) #, ctype='float' )
print("number of elements of periodic YaspGrid grid:",p_grid.size(0))
p_grid.plot()
runOnGrid(p_grid)
