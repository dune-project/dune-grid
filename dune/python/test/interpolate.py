# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

import time, math, numpy, io
import dune
from dune.generator import algorithm
import dune.geometry
from dune.grid import cartesianDomain
from dune.grid import yaspGrid as gridView

def testPyQuad(view,rules,error):
    start = time.time()
    l2norm2 = 0
    for e in view.elements:
        hatxs, hatws = rules(e.type).get()
        weights = hatws * e.geometry.integrationElement(hatxs)
        l2norm2 += numpy.sum(error(e, hatxs)**2 * weights, axis=-1)
    # print("Python:",math.sqrt(l2norm2),flush=True)
    # print("time used:", round(time.time()-start,2),flush=True)
    return l2norm2

code="""
#include <cstddef>
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/numpy.h>

template< class GridView, class Rules, class GF >
double l2norm2 ( const GridView &gridView, const Rules &rules, const GF& gf )
{
  auto lf = localFunction( gf );
  double l2norm2 = 0;
  for( const auto &entity : elements( gridView ) )
  {
    const auto geo = entity.geometry();
    typedef typename decltype(geo)::LocalCoordinate LocalCoordinate;
    lf.bind( entity );
    pybind11::object pyrule = rules( geo.type() );
    pybind11::object pyPW = pyrule.attr("get")();
    auto pointsWeights = pyPW.template cast<
       std::pair<pybind11::array_t<double>,
                 pybind11::array_t<double>> >();

    const auto &valuesArray = lf( pointsWeights.first ).template cast< pybind11::array_t< double > >();
    // check shape here...
    auto values = valuesArray.template unchecked< 1 >();
    for( std::size_t i = 0, sz = pointsWeights.second.size(); i < sz; ++i )
    {
      LocalCoordinate hatx(0);
      for (std::size_t c=0;c<LocalCoordinate::size();++c)
        hatx[c] = pointsWeights.first.at(c,i);
      double weight = pointsWeights.second.at( i ) * geo.integrationElement( hatx );
      l2norm2 += (values[ i ] * values[ i ]) * weight;
    }
    lf.unbind();
  }
  return l2norm2;
}
"""
def testCppQuad(view,rules,error):
    algo = algorithm.load('l2norm2', io.StringIO(code), view, rules, error)
    start = time.time()
    l2norm2 = algo(view,rules,error)
    # print("C++:",math.sqrt(l2norm2),flush=True)
    # print("time used:", round(time.time()-start,2),flush=True)
    return l2norm2

if __name__ == "__main__":
    domain = cartesianDomain([0, 0], [1, 1], [50, 50])
    view = gridView(domain)

    @dune.grid.gridFunction(view)
    def function(x):
        return numpy.cos(numpy.pi*x[0])*numpy.cos(numpy.pi*x[1])
    def interpolate(grid):
        mapper = grid.mapper({dune.geometry.vertex: 1})
        data = numpy.zeros(mapper.size)
        for v in grid.vertices:
            data[mapper.index(v)] = function(v.geometry.center)
        return mapper, data
    mapper, data = interpolate(view)
    @dune.grid.gridFunction(view)
    def p12dEvaluate(element,x):
        indices = mapper(element)
        if element.type == dune.geometry.cube(2):
            bary = (1-x[0])*(1-x[1]), x[0]*(1-x[1]), (1-x[0])*x[1], x[0]*x[1]
        elif element.type == dune.geometry.simplex(2):
            bary = 1-x[0]-x[1], x[0], x[1]
        else:
            raise RuntimeError("basis functions not implemented for type "+str(element.type))
        assert len(indices) == len(bary)
        return sum( b*data[i] for b,i in zip(bary,indices) )

    @dune.grid.gridFunction(view)
    def error(element,x):
        return p12dEvaluate(element,x)-function(element,x)

    rules = dune.geometry.quadratureRules(5)
    value1 = testPyQuad(view,rules,error)
    value2 = testCppQuad(view,rules,error)
    assert abs(value1-value2) < 1e-14
    assert value1 < 1e-6
