// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_TEST_EVOLVE_HH
#define DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_TEST_EVOLVE_HH

#include <algorithm>
#include <cstddef>
#include <future>
#include <limits>
#include <thread>

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/utility/gridviewpartitioning/strided.hh>

struct EvolveOnInteriorIntersectionPlain
{
  template<class Intersection, class GV, class Mapper, class V>
  void operator()(const Intersection &is, const GV &gv, const Mapper &mapper,
                  int indexi, double factor, double normalFlux, const V &c,
                  V &update) const
  {
    // access neighbor
    int indexj = mapper.map(*is.outside());

    if (normalFlux<0)                // inflow
      update[indexi] -= c[indexj]*factor;
    else                         // outflow
      update[indexi] -= c[indexi]*factor;
  }
};

struct EvolveOnInteriorIntersectionOptimized
{
  template<class Intersection, class GV, class Mapper, class V>
  void operator()(const Intersection &is, const GV &gv, const Mapper &mapper,
                  int indexi, double factor, double normalFlux, const V &c,
                  V &update) const
  {
    // access neighbor
    typename GV::template Codim<0>::EntityPointer outside = is.outside();
    int indexj = mapper.map(*outside);
    double upwindc;
    if(normalFlux < 0)
      upwindc = c[indexj];
    else
      upwindc = c[indexi];

    if(gv.contains(*outside)) {
      // optimized
      if(indexi < indexj)
      {
        // compute factor in neighbor
        double nbfactor = normalFlux/outside->geometry().volume();

        update[indexi] -= upwindc*factor;
        update[indexj] += upwindc*nbfactor;
      }
    }
    else
    {
      // one-sided assembly on thread-partition boundaries
      update[indexi] -= upwindc*factor;
    }
  }
};

template<class Intersection, class GV, class Entity, class Mapper, class V,
         class EvolveOnInteriorIntersection>
void evolveOnIntersection(const Intersection &is, const GV &gv,
                          const Mapper &mapper, const Entity &e, int indexi,
                          double volume, double &sumfactor, double t,
                          const V &c, V &update,
                          const EvolveOnInteriorIntersection &evolveOnII)
{
  // intersection geometry
  typedef typename Intersection::Geometry IntersectionGeometry;

  typedef typename IntersectionGeometry::ctype ct;
  static const std::size_t dimworld = IntersectionGeometry::dimensionworld;

  // get geometry type of face
  const IntersectionGeometry &igeo = is.geometry();

  // center of face in global coordinates
  Dune::FieldVector<ct,dimworld> faceglobal = igeo.center();

  // evaluate velocity at face center
  double normalFlux =
    ( u(faceglobal,t) * is.centerUnitOuterNormal() ) * igeo.volume();

  // compute factor occuring in flux formula
  double factor = normalFlux/volume;

  // for time step calculation
  if (factor>=0) sumfactor += factor;

  // handle interior face
  if (is.neighbor())             // "correct" version /*@\label{evh:neighbor}@*/
    evolveOnII(is, gv, mapper, indexi, factor, normalFlux, c, update);

  // handle boundary face
  if (is.boundary())                               /*@\label{evh:bndry}@*/
  {
    if (normalFlux<0)                 // inflow, apply boundary condition
      update[indexi] -= b(faceglobal,t)*factor;
    else                 // outflow
      update[indexi] -= c[indexi]*factor;
  }
}             // end all intersections             /*@\label{evh:flux1}@*/

template<class Entity, class GridView, class Mapper, class V,
         class EvolveOnInteriorIntersection>
void evolveOnEntity(const Entity &e, const GridView &gridView,
                    const Mapper &mapper, double t, const V &c, V &update,
                    double &dt, const EvolveOnInteriorIntersection &evolveOnII)
{
  // cell geometry
  const typename Entity::Geometry geo = e.geometry();


  // cell volume, assume linear map here
  double volume = geo.volume();

  // cell index
  int indexi = mapper.map(e);

  // variable to compute sum of positive factors
  double sumfactor = 0.0;

  // run through all intersections with neighbors and boundary
  typedef typename GridView::IntersectionIterator IntersectionIterator;
  IntersectionIterator isend = gridView.iend(e);       /*@\label{evh:flux0}@*/
  for (IntersectionIterator is = gridView.ibegin(e); is!=isend; ++is)
    evolveOnIntersection(*is, gridView, mapper, e, indexi, volume, sumfactor,
                         t, c, update, evolveOnII);

  // compute dt restriction
  dt = std::min(dt,1.0/sumfactor);                   /*@\label{evh:dt}@*/

}       // end grid traversal                        /*@\label{evh:loop1}@*/

struct SeqEvolve
{
  template<class G, class M, class V, class EvolveOnInteriorIntersection>
  void operator()(const G& grid, const M& mapper, V& c, double t,
                  double& dt,
                  const EvolveOnInteriorIntersection &evolveOnII) const
  {
    // type of grid view on leaf part
    typedef typename G::LeafGridView GridView;

    // element iterator type
    typedef typename GridView::template Codim<0>::Iterator Iterator;

    // get grid view on leaf part
    GridView gridView = grid.leafView();

    // allocate a temporary vector for the update
    V update(c.size());                                  /*@\label{evh:update}@*/
    for (typename V::size_type i=0; i<c.size(); i++) update[i] = 0;

    // initialize dt very large
    dt = 1E100;

    // compute update vector and optimum dt in one grid traversal
    Iterator endit = gridView.template end<0>();     /*@\label{evh:loop0}@*/
    for (Iterator it = gridView.template begin<0>(); it!=endit; ++it)
      evolveOnEntity(*it, gridView, mapper, t, c, update, dt, evolveOnII);

    // scale dt with safety factor
    dt *= 0.99;                                          /*@\label{evh:.99}@*/

    // update the concentration vector
    for (unsigned int i=0; i<c.size(); ++i)
      c[i] += dt*update[i];                              /*@\label{evh:updc}@*/

    return;
  }
};

#if HAVE_TBB
// evolve with the help of TBB.  Use up to 2 threads.
struct TBBEvolve
{
  template<class G, class M, class V, class EvolveOnInteriorIntersection>
  void operator()(const G& grid, const M& mapper, V& c, double t,
                  double& dt,
                  const EvolveOnInteriorIntersection &evolveOnII) const
  {
    // type of grid view on leaf part
    typedef typename G::LeafGridView HostView;
    typedef Dune::StridedGridViewPartitioning<HostView> Partitioning;

    // get grid view on leaf part
    HostView hostView = grid.leafView();
    Partitioning partitioning(hostView, 2);

    // allocate a temporary vector for the update
    V update(c.size(), 0);                                  /*@\label{evh:update}@*/

    dt = tbb::parallel_reduce
      ( tbb::blocked_range<int>(0,2),
        std::numeric_limits<double>::infinity(),
        [&](const tbb::blocked_range<int> &range, double mydt)
        {
          for(int workerId = range.begin(); workerId != range.end();
              ++workerId)
          {
            // compute update vector and optimum dt in one grid traversal
            auto gridView = partitioning.gridView(workerId);
            auto endit = gridView.template end<0>();
            for (auto it = gridView.template begin<0>(); it!=endit; ++it)
              evolveOnEntity(*it, gridView, mapper, t, c, update, mydt,
                             evolveOnII);
          }
          return mydt;
        },
        [](double a, double b) { return std::min(a,b); });

    // scale dt with safety factor
    dt *= 0.99;                                          /*@\label{evh:.99}@*/

    // update the concentration vector
    for (unsigned int i=0; i<c.size(); ++i)
      c[i] += dt*update[i];                              /*@\label{evh:updc}@*/
  }
};
#endif // HAVE_TBB

#endif // DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_TEST_EVOLVE_HH
