import math

# constructed a Cartesian grid
from dune.grid import structuredGrid
grid = structuredGrid([-1,-1],[1,1],[10,10])
print("number of elements of Cartesian grid:",grid.size(0))
grid.plot()

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
