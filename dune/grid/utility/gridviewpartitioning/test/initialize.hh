// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_TEST_INITIALIZE_HH
#define DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_TEST_INITIALIZE_HH

#include <dune/common/fvector.hh>

//! initialize the vector of unknowns with initial value
template<class G, class M, class V>
void initialize (const G& grid, const M& mapper, V& c)
{
  // first we extract the dimensions of the grid
  //const int dim = G::dimension;
  const int dimworld = G::dimensionworld;

  // type used for coordinates in the grid
  typedef typename G::ctype ct;

  // type of grid view on leaf part
  typedef typename G::LeafGridView GridView;

  // leaf iterator type
  typedef typename GridView::template Codim<0>::Iterator LeafIterator;

  // geometry type
  typedef typename LeafIterator::Entity::Geometry Geometry;

  // get grid view on leaf part
  GridView gridView = grid.leafView();

  // iterate through leaf grid an evaluate c0 at cell center
  LeafIterator endit = gridView.template end<0>();
  for (LeafIterator it = gridView.template begin<0>(); it!=endit; ++it)
  {
    // get geometry
    const Geometry geo = it->geometry();

    // get global coordinate of cell center
    Dune::FieldVector<ct,dimworld> global = geo.center();

    // initialize cell concentration
    c[mapper.map(*it)] = c0(global);
  }
}

#endif // DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_TEST_INITIALIZE_HH
