# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import numpy
from dune.grid import gridFunction, structuredGrid

fcalls = 0
gcalls = 0
def testGF_first(gridView):
    global fcalls
    global gcalls
    @gridFunction(gridView)
    def f(x):
        global fcalls
        fcalls += 1
        return x[0]*x[1]
    @gridFunction(gridView)
    def g(e,x):
        global gcalls
        gcalls += 1
        return e.geometry.toGlobal(x)

    fcalls = 0
    gcalls = 0
    e = gridView.elements.__next__()
    xLoc = numpy.array([[0,0.1,0.2,0.3],[0,0.4,0.6,0.8]])
    xGlb = e.geometry.toGlobal(xLoc)

    fcalls = 0
    gcalls = 0
    y=f(xGlb)
    # print( y, fcalls)
    assert fcalls == 1
    y = g(e, xLoc)
    # print( y, gcalls)
    assert gcalls == 1

    fcalls = 0
    y = f(e,xLoc)
    # print( y, fcalls)
    assert fcalls == 1

    fcalls = 0
    gcalls = 0
    lf = f.localFunction()
    lg = g.localFunction()
    lf.bind(e)
    lg.bind(e)
    y = lf(xLoc)
    # print( y, fcalls)
    assert fcalls == 1
    y = lg(xLoc)
    # print( y, gcalls)
    assert gcalls == 1
    lg.unbind()
    lf.unbind()

if __name__ == "__main__":
    gridView = structuredGrid([0,0],[1,1],[10,10])
    testGF_first(gridView)
