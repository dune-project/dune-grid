# constructed a Cartesian grid
from dune.grid import structuredGrid
grid = structuredGrid([-1,-1],[1,1],[10,10])
grid.plot()
# define a grid function and visualize
from dune.grid import gridFunction
@gridFunction(grid)
def f(x):
    return math.cos(2.*math.pi/(1+abs(x[0]*x[1])))
f.plot()
