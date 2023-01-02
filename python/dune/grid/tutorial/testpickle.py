# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import sys, numpy
import dune.common.pickle
try:
    from dune.alugrid import aluConformGrid as view
    from dune.grid import Marker
except ImportError:
    from dune.grid import yaspGrid as view
    Marker = None
from dune.grid import cartesianDomain, gridFunction

def test(fileName):
    # at the moment only locally defined gridfunctions can be pickled also
    # function needs to be defined outside of 'main' for unpickling
    from picklefunc import globalF, localF, TimeDependent
    grid = view( cartesianDomain([-2,-2,-2],[2,2,2],[3,3,3]) )
    grid.hierarchicalGrid.globalRefine(2)
    if Marker is not None:
        def indicator(e,y):
            x = e.geometry.toGlobal(y)
            return numpy.exp( -(x*x-1)**2 )
        for i in range(5):
            grid.hierarchicalGrid.mark(lambda e:
                 Marker.refine if indicator(e,[1./3.,1./3.]) > 0.9
                               else Marker.coarsen)
            grid.hierarchicalGrid.adapt()
    print("size of adapted grid:", grid.size(0))

    # make localF into a gridFunction
    gf = gridFunction(grid, name="gf", order=3)(localF)
    # glabal version does not work: gf = gridFunction(grid, name="gf", order=3)(globalF)
    with open(fileName+".dbf","wb") as f:
        dune.common.pickle.dump([gf],f)

    func = TimeDependent()
    timeGf = gridFunction(grid, name="gf", order=3)(func.__call__)

    series = dune.common.pickle.SeriesPickler(fileName+"_ts", [timeGf])

    tsps = [0,0.1,0.4,0.8,2,4]
    for i,tsp in enumerate(tsps):
        func.t = tsp
        series.dump({"time":tsp})


fbase = "dump.test"
if len(sys.argv)<2 or (not sys.argv[1] == 'load'):
    test(fbase)
else:
    with open(fbase+".dbf","rb") as f:
        dump = dune.common.pickle.load(f)
    dump[0].plot()
