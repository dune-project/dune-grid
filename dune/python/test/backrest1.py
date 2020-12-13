import pickle, numpy
import dune.generator

def backup():
    from dune.grid import structuredGrid
    grid = structuredGrid([0,0],[1,1],[2,2])
    grid.hierarchicalGrid.globalRefine(2)
    a = numpy.array([1,2,3])

    pickle.dump([a,"hallo",grid.hierarchicalGrid,10], open("dumpA",'wb'))
    return grid

grid = backup()

def restore():
    return pickle.load(open("dumpA","rb"))

[b,string,otherHGrid,value] = restore()
otherGrid = otherHGrid.leafView

print("leaf after refine", grid.size(0),otherGrid.size(0))
print("level 1 after refine",
      grid.hierarchicalGrid.levelView(1).size(0),
      otherGrid.hierarchicalGrid.levelView(1).size(0))
otherGrid.hierarchicalGrid.globalRefine(-2)
print("coarsen other", grid.size(0),otherGrid.size(0))
grid.hierarchicalGrid.globalRefine(-2)
print("coarsen original", grid.size(0),otherGrid.size(0))
