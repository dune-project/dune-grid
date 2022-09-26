// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef COORDFUNCTION

#if defined __GNUC__ and not defined __clang__ and not defined __ICC
  #define GCCPOOL
#endif

#include <dune/common/timer.hh>

#include <dune/common/poolallocator.hh>
#include <dune/common/debugallocator.hh>
#ifdef GCCPOOL
#include <ext/pool_allocator.h>
#endif

#include <dune/grid/geometrygrid.hh>
#include <dune/grid/geometrygrid/cachedcoordfunction.hh>
#include <dune/grid/io/file/dgfparser.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/dgfparser/dgfgeogrid.hh>

#include "functions.hh"

#include "gridcheck.hh"
#include "checkcommunicate.hh"
#include "checkgeometryinfather.hh"
#include "checkintersectionit.hh"
#include "checkiterators.hh"
#include "checkpartition.hh"
#include "checkgeometry.hh"

namespace Dune
{

  template< int dim, int dimworld >
  class AlbertaGrid;

}


template< int dim, int dimworld >
struct EnableLevelIntersectionIteratorCheck< Dune::AlbertaGrid< dim, dimworld > >
{
  static const bool v = false;
};

template< class HostGrid, class CoordFunction >
struct EnableLevelIntersectionIteratorCheck< Dune::GeometryGrid< HostGrid, CoordFunction > >
{
  static const bool v = EnableLevelIntersectionIteratorCheck< HostGrid >::v;
};


typedef Dune::COORDFUNCTION AnalyticalCoordFunction;

typedef GRIDTYPE Grid;

#if CACHECOORDFUNCTION
typedef Dune::CachedCoordFunction< Grid, AnalyticalCoordFunction > CoordFunction;
#else
typedef AnalyticalCoordFunction CoordFunction;
#endif

typedef Dune::GeometryGrid< Grid, CoordFunction > GeometryGrid;
typedef Dune::GeometryGrid< Grid, CoordFunction, Dune::PoolAllocator< char, 16384 > > GeometryGridWithPoolAllocator;
#ifdef GCCPOOL
typedef Dune::GeometryGrid< Grid, CoordFunction, __gnu_cxx::__pool_alloc<char> > GeometryGridWithGCCPoolAllocator;
#endif
typedef Dune::GeometryGrid< Grid, CoordFunction, Dune::DebugAllocator<char> > GeometryGridWithDebugAllocator;

template <class HostGridView>
struct DeformationFunction
  : public Dune::DiscreteCoordFunction< double, HostGridView::dimension, DeformationFunction<HostGridView> >
{
  static const int dim = HostGridView::dimension;

  DeformationFunction(const HostGridView& gridView,
                      const std::vector<Dune::FieldVector<double,dim> >& deformedPosition)
    : indexSet_(gridView.indexSet()),
      deformedPosition_(deformedPosition)
  {}

  /** \brief Evaluate the position at the corner of an entity
   *
   * \note This method needs to work for entities of all codimensions,
   *    not just for elements!
   */
  template <class HostEntity>
  void evaluate (const HostEntity& hostEntity, unsigned int corner,
                 Dune::FieldVector<double,dim> &y ) const
  {
    auto idx = indexSet_.subIndex(hostEntity, corner,dim);
    y = deformedPosition_[idx];
  }

private:
  const typename HostGridView::IndexSet& indexSet_;
  const std::vector< Dune::FieldVector<double,dim> > deformedPosition_;
};

template <class GeometryGridType>
void test(const std::string& gridfile)
{
  Dune::GridPtr< GeometryGridType > pgeogrid(gridfile);
  GeometryGridType &geogrid = *pgeogrid;

  geogrid.globalRefine( 1 );
  geogrid.loadBalance();

  std::cerr << "Checking grid..." << std::endl;
  gridcheck( geogrid );

  std::cerr << "Checking geometry... " << std::endl;
  Dune::GeometryChecker< GeometryGridType > checker;

  checker.checkGeometry( geogrid.leafGridView() );
  for( int i = 0; i <= geogrid.maxLevel(); ++i )
    checker.checkGeometry( geogrid.levelGridView( i ) );

  std::cerr << "Checking geometry in father..." << std::endl;
  checkGeometryInFather( geogrid );
  std::cerr << "Checking intersections..." << std::endl;
  checkIntersectionIterator( geogrid, !EnableLevelIntersectionIteratorCheck< Grid >::v );

  checkIterators( geogrid.leafGridView() );
  for( int i = 0; i <= geogrid.maxLevel(); ++i )
    checkIterators( geogrid.levelGridView( i ) );

  checkPartitionType( geogrid.leafGridView() );
  for( int i = 0; i <= geogrid.maxLevel(); ++i )
    checkPartitionType( geogrid.levelGridView( i ) );

  std::cerr << "Checking geometry lifetime..." << std::endl;
  checkGeometryLifetime( geogrid.leafGridView() );

  std::cerr << "Checking communication..." << std::endl;
  checkCommunication( geogrid, -1, std::cout );
  if( EnableLevelIntersectionIteratorCheck< Grid >::v )
  {
    for( int i = 0; i <= geogrid.maxLevel(); ++i )
      checkCommunication( geogrid, i, std::cout );
  }

}

void testNestedGeometryGrid(const std::string& gridfile) {

  using NestedGeometryGrid = Dune::GeometryGrid< GeometryGrid, Dune::IdenticalCoordFunction< Grid::ctype, Grid::dimensionworld > >;
  using NestedGgLeafView = NestedGeometryGrid::LeafGridView;

  Dune::GridPtr< NestedGeometryGrid > pgeogrid(gridfile);
  NestedGeometryGrid &geogrid = *pgeogrid;
  geogrid.globalRefine( 1 );

  // creating different variables for storing grid views
  NestedGgLeafView gv1 = geogrid.leafGridView();
  std::shared_ptr<NestedGgLeafView> pgv2;
  NestedGgLeafView gv3 = geogrid.leafGridView();

  // use geometry grid views copy constructor and copy assignment operator
  {
    NestedGgLeafView tmpGv = geogrid.leafGridView();
    // copying/assigning temporary grid view would be dangerous with default implementation of copy constructor or copy assignment operator from geomrtry grid view
    pgv2 = std::make_shared<NestedGgLeafView> (tmpGv);
    gv3 = tmpGv;
  }
  // The following test would now be guaranteed to fail with the default implementation of geometry grid views copy constructor
  // assert(&(pgv2->indexSet().hostIndexSet()) == &(pgv2->impl().hostGridView().indexSet()))
  // but this would not compile since the member functions impl() and hostIndexSet() are protected/private.
  // So we do the following

  // get index sets and do some arbitrary test to check whether they work correctly
  auto const& is1 = gv1.indexSet();
  auto const& is2 = pgv2->indexSet();
  auto const& is3 = gv3.indexSet();
  for(auto const& e : elements(gv1)) {
    bool idxCorrect = (is1.index(e)==is2.index(e));
    idxCorrect &= (is1.index(e)==is3.index(e));
    assert(idxCorrect);
  }
}

int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance( argc, argv );

  std::string gridfile = DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/cube-2.dgf";
  if(argc >= 2)
  {
    gridfile = argv[1];
  }

  Dune::Timer watch;

  watch.reset();
  test<GeometryGrid>(gridfile);
  std::cout << "=== GeometryGrid took " << watch.elapsed() << " seconds\n";

  watch.reset();
  testNestedGeometryGrid(gridfile);
  std::cout << "=== NestedGeometryGrid took " << watch.elapsed() << " seconds\n";

  // compile, but do not actually call, because it is not working yet
  if (false)
  {
    watch.reset();
    test<GeometryGridWithPoolAllocator>(gridfile);
    std::cout << "=== GeometryGridWithPoolAllocator took " << watch.elapsed() << " seconds\n";
  }

#ifdef GCCPOOL
  watch.reset();
  test<GeometryGridWithGCCPoolAllocator>(gridfile);
  std::cout << "=== GeometryGridWithGCCPoolAllocator took " << watch.elapsed() << " seconds\n";
#endif

  watch.reset();
  test<GeometryGridWithDebugAllocator>(gridfile);
  std::cout << "=== GeometryGridWithDebugAllocator took " << watch.elapsed() << " seconds\n";

  // Check with a discrete coordinate function
  Dune::GridPtr< Grid > hostGridPtr(gridfile);
  auto& hostGrid = *hostGridPtr;
  const int dim = Grid::dimension;

  typedef Grid::LeafGridView HostGridView;
  const HostGridView& hostGridView = hostGrid.leafGridView();

  std::vector<Dune::FieldVector<double,dim> > vertexPositions(hostGridView.size(dim));

  for (const auto& vertex : vertices(hostGridView))
    vertexPositions[hostGridView.indexSet().index(vertex)] = vertex.geometry().corner(0);

  typedef Dune::GeometryGrid<Grid, DeformationFunction<HostGridView> > DiscretelyTransformedGrid;
  DeformationFunction<HostGridView> discreteTransformation(hostGridView, vertexPositions);
  DiscretelyTransformedGrid discretelyTransformedGrid(hostGrid, discreteTransformation);

  gridcheck(discretelyTransformedGrid);

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch( const std::exception &e )
{
  std::cerr << e.what() << std::endl;
  return 1;
}
catch( ... )
{
  std::cerr << "Unknown exception raised." << std::endl;
  return 1;
}

#else
#error "COORDFUNCTION not defined (e.g., Helix, Circle; see functions.hh)"
#endif
