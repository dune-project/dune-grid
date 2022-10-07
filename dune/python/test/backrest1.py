# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import pickle, numpy
import dune.generator

def backup():
    from dune.grid import structuredGrid
    grid = structuredGrid([0,0],[1,1],[2,2])
    grid.hierarchicalGrid.globalRefine(2)
    a = numpy.array([1,2,3])

    pickle.dump([a,"hallo",grid.hierarchicalGrid,10], open("dumpA",'wb'))
    return grid

def restore():
    return pickle.load(open("dumpA","rb"))

class Test:
  def __init__(self,g,og):
      # note: only classes containing HGrids can be pickeled
      self.hg = g.hierarchicalGrid
      self.hog = og.hierarchicalGrid
  def run(self):
      return self.hg.leafView.size(0) == self.hog.leafView.size(0)


if __name__ == "__main__":
    grid = backup()

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

    test = Test(grid,otherGrid)
    print("test:",test.run())
    pickle.dump(test,open("dumpB","wb"))

    del test, grid, otherGrid

    test = pickle.load(open("dumpB","rb"))
    print("test:",test.run())
