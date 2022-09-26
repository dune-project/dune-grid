# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

from io import StringIO
import math
from dune.common import FieldVector
from dune.grid import gridFunction, structuredGrid

codeFunc = """
#include <cmath>
#include <dune/common/fvector.hh>

template <class GridView>
auto myFunction(double a)
{
  return [a](const auto& en,const auto& x) -> auto
  {
    auto y = en.geometry().global(x);
    return std::sin(a*M_PI*(y[0]+y[1]));
  };
}
"""
codeVecFunc = """
#include <dune/common/fvector.hh>

template <class GridView, class GF>
auto myVecFunction(Dune::FieldVector<double,1> &a, const GF &gf)
{
  return [&a,lgf=localFunction(gf)](const auto& en,const auto& x) mutable -> auto
  {
    lgf.bind(en);
    auto v = lgf(x);
    auto y = en.geometry().global(x);
    return Dune::FieldVector<double,2>{v[0],(y[0]-0.5)*a[0]};
  };
}
"""

def testGF_second(gridView):
    gf1 = gridView.function(lambda e,x:\
              math.sin(math.pi*(e.geometry.toGlobal(x)[0]+e.geometry.toGlobal(x)[1])))

    if True:
        a = 2.
        gf1 = gridView.function(lambda e,x:\
                math.sin(a*math.pi*(e.geometry.toGlobal(x)[0]+e.geometry.toGlobal(x)[1])),
                name="gf1")
        lgf1 = gf1.localFunction()
        average1 = 0
        for e in gridView.elements:
            lgf1.bind(e)
            average1 += lgf1([0.5,0.5])*e.geometry.volume
        # print(average1)
        # gf1.plot()
        gf2 = gridView.function("myFunction",StringIO(codeFunc),a,name="gf2")
        lgf2 = gf2.localFunction()
        average2 = 0
        for e in gridView.elements:
           lgf2.bind(e)
           average2 += lgf2([0.5,0.5])*e.geometry.volume
        # print(average2)
        # gf2.plot()
        # assert abs(average1-average2)<1e-12
        diff = 0
        for e in gridView.elements:
           lgf1.bind(e)
           lgf2.bind(e)
           diff += abs(lgf1([0.5,0.5])-lgf2([0.5,0.5]))
        assert diff<1e-12

    if True:
        gf1 = gridView.function(lambda e,x:\
                [math.sin(2*math.pi*(e.geometry.toGlobal(x)[0]+e.geometry.toGlobal(x)[1])),\
                 (e.geometry.toGlobal(x)[0]-0.5)*2],
                name="gf1")
        lgf1 = gf1.localFunction()
        average1 = 0
        for e in gridView.elements:
            lgf1.bind(e)
            average1 += sum(lgf1([0.5,0.5]))*e.geometry.volume
        # print(average1)
        # gf1.plot()
        a = FieldVector([2])
        gf2 = gridView.function("myVecFunction",StringIO(codeVecFunc),a,gf1,name="gf2")
        lgf2 = gf2.localFunction()
        average2 = 0
        for e in gridView.elements:
            lgf2.bind(e)
            average2 += sum(lgf2([0.5,0.5]))*e.geometry.volume
        # print(average2)
        # gf2.plot()
        assert abs(average1-average2)<1e-12
        diff = 0
        for e in gridView.elements:
            lgf1.bind(e)
            lgf2.bind(e)
            diff += abs(lgf1([0.5,0.5]).two_norm-lgf2([0.5,0.5]).two_norm)
        assert diff<1e-12
        a[0] = 3
        diff = 0
        for e in gridView.elements:
            lgf1.bind(e)
            lgf2.bind(e)
            v = lgf1([0.5,0.5])
            v[1] *= 3./2.
            diff += abs(v.two_norm-lgf2([0.5,0.5]).two_norm)
        assert diff<1e-12

        gridView.writeVTK("test_gf",pointdata=[gf1,gf2])

    if False:
        a = 2.
        @gridFunction(gridView)
        def gf1(e,x):
            return math.sin(a*math.pi*(e.geometry.toGlobal(x)[0]+e.geometry.toGlobal(x)[1]))
        lgf1 = gf1.localFunction()
        average1 = 0
        for e in gridView.elements:
            lgf1.bind(e)
            average1 += lgf1([0.5,0.5])*e.geometry.volume
        # print(average1)
        # gf1.plot()
        @gridFunction(gridView)
        def gf2(x):
            return math.sin(a*math.pi*(x[0]+x[1]))
        gf2 = gridView.function("myFunction",StringIO(codeFunc),a,name="gf2")
        lgf2 = gf2.localFunction()
        average2 = 0
        for e in gridView.elements:
            lgf2.bind(e)
            average2 += lgf2([0.5,0.5])*e.geometry.volume
        # print(average2)
        # gf2.plot()
        assert abs(average1-average2)<1e-12
        diff = 0
        for e in gridView.elements:
            lgf1.bind(e)
            lgf2.bind(e)
            diff += abs(lgf1([0.5,0.5])-lgf2([0.5,0.5]))
        assert diff<1e-12

if __name__ == "__main__":
    gridView = structuredGrid([0,0],[1,1],[10,10])
    testGF_second(gridView)
