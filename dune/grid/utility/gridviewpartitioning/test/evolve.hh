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

#include <dune/grid/utility/filteringentityset.hh>

struct EvolveOnInteriorIntersectionPlain
{
  template<class EntitySet, class Mapper, class V>
  void operator()(const typename EntitySet::GridView::Intersection &is,
                  const EntitySet &entiySet, const Mapper &mapper, int indexi,
                  double factor, double normalFlux, const V &c,
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
  template<class EntitySet, class Mapper, class V>
  void operator()(const typename EntitySet::GridView::Intersection &is,
                  const EntitySet &entitySet, const Mapper &mapper, int indexi,
                  double factor, double normalFlux, const V &c,
                  V &update) const
  {
    // access neighbor
    auto outside = is.outside();
    int indexj = mapper.map(*outside);
    double upwindc;
    if(normalFlux < 0)
      upwindc = c[indexj];
    else
      upwindc = c[indexi];

    if(entitySet.contains(*outside)) {
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

template<class EntitySet, class Mapper, class V,
         class EvolveOnInteriorIntersection>
void evolveOnIntersection(const typename EntitySet::GridView::Intersection &is,
                          const EntitySet &entitySet, const Mapper &mapper,
                          const typename EntitySet::Element &e, int indexi,
                          double volume, double &sumfactor, double t,
                          const V &c, V &update,
                          const EvolveOnInteriorIntersection &evolveOnII)
{
  // get geometry of face
  const auto &igeo = is.geometry();

  // center of face in global coordinates
  auto faceglobal = igeo.center();

  // evaluate velocity at face center
  double normalFlux =
    ( u(faceglobal,t) * is.centerUnitOuterNormal() ) * igeo.volume();

  // compute factor occuring in flux formula
  double factor = normalFlux/volume;

  // for time step calculation
  if (factor>=0) sumfactor += factor;

  // handle interior face
  if (is.neighbor())             // "correct" version /*@\label{evh:neighbor}@*/
    evolveOnII(is, entitySet, mapper, indexi, factor, normalFlux, c, update);

  // handle boundary face
  if (is.boundary())                               /*@\label{evh:bndry}@*/
  {
    if (normalFlux<0)                 // inflow, apply boundary condition
      update[indexi] -= b(faceglobal,t)*factor;
    else                 // outflow
      update[indexi] -= c[indexi]*factor;
  }
}             // end all intersections             /*@\label{evh:flux1}@*/

template<class EntitySet, class Mapper, class V,
         class EvolveOnInteriorIntersection>
void evolveOnEntity(const typename EntitySet::Element &e,
                    const EntitySet &entitySet, const Mapper &mapper, double t,
                    const V &c, V &update, double &dt,
                    const EvolveOnInteriorIntersection &evolveOnII)
{
  typedef typename EntitySet::Element Entity;
  typedef typename EntitySet::GridView GridView;
  const GridView &gridView = entitySet.gridView();

  // cell geometry
  const typename Entity::Geometry& geo = e.geometry();

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
    evolveOnIntersection(*is, entitySet, mapper, e, indexi, volume, sumfactor,
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
    typedef Dune::StridedEntitySet<GridView, 0> EntitySet;

    // get grid view on leaf part
    GridView gridView = grid.leafView();
    EntitySet entitySet(gridView);

    // allocate a temporary vector for the update
    V update(c.size());                                  /*@\label{evh:update}@*/
    for (typename V::size_type i=0; i<c.size(); i++) update[i] = 0;

    // initialize dt very large
    dt = 1E100;

    // compute update vector and optimum dt in one grid traversal
    for (const auto &e : entitySet)
      evolveOnEntity(e, entitySet, mapper, t, c, update, dt, evolveOnII);

    // scale dt with safety factor
    dt *= 0.99;                                          /*@\label{evh:.99}@*/

    // update the concentration vector
    for (unsigned int i=0; i<c.size(); ++i)
      c[i] += dt*update[i];                              /*@\label{evh:updc}@*/

    return;
  }
};

#if HAVE_TBB
// evolve with the help of TBB.
class TBBEvolve
{
  std::size_t maxStride_;

public:
  TBBEvolve(std::size_t maxStride = 0) :
    maxStride_(maxStride)
  { }

  template<class G, class M, class V, class EvolveOnInteriorIntersection>
  void operator()(const G& grid, const M& mapper, V& c, double t,
                  double& dt,
                  const EvolveOnInteriorIntersection &evolveOnII) const
  {
    // type of grid view on leaf part
    typedef typename G::LeafGridView GridView;
    typedef Dune::StridedEntitySet<GridView, 0> EntitySet;

    // get grid view on leaf part
    GridView gridView = grid.leafView();

    // allocate a temporary vector for the update
    V update(c.size(), 0);                                  /*@\label{evh:update}@*/

    dt = tbb::parallel_reduce
      ( EntitySet(gridView, maxStride_),
        std::numeric_limits<double>::infinity(),
        [&](const EntitySet &entitySet, double mydt)
        {
          // compute update vector and optimum dt in one grid traversal
          for (const auto &e : entitySet)
            evolveOnEntity(e, entitySet, mapper, t, c, update, mydt,
                           evolveOnII);
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
