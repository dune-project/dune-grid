# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import math
import dune.common
import dune.grid
import sys

# initialize MPI *if present*
try:
    import mpi4py
except ImportError as exception:
    pass

base = "../../../doc/grids/gmsh/"
reader = (dune.grid.reader.gmsh, base+"circle2ndorder.msh")
grid = dune.grid.ugGrid(reader, dimgrid=2)
# we either have UG, ot this test should not be run at all 
assert(grid)

exact = 6.27593206157460
# import dune.plotting

def integrateBoundary(grid):
    print(grid)
    sum = 0.0
    for intersection in grid.boundaryIntersections:
        sum += intersection.geometry.volume
    return sum


len = integrateBoundary(grid)
eoc = 0.0
for l in range(5):
    grid.hierarchicalGrid.globalRefine(1)
    len_refined = integrateBoundary(grid)
    err0 = abs(len-exact)
    err  = abs(len_refined-exact)
    eoc  = math.log2(err0/err) # h0/h = 2
    len = len_refined;
    print("boundary length: ", len, "\t", exact, "\terr=", err, "\teoc=", eoc)

# we expect 2nd order convergence
if math.isclose(eoc,2.0,abs_tol=0.02):
    print("Geometry approximation: 2d order convergence")
else:
    print("Geometry approximation does not show expected 2d order convergence")
    sys.exit(-1)

sys.exit(0)
