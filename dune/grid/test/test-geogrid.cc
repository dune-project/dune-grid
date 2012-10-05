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
#include <dune/grid/io/file/dgfparser/dgfgeogrid.hh>

#include "functions.hh"

#include "gridcheck.cc"
#include "checkcommunicate.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkiterators.cc"
#include "checkpartition.cc"
#include "checkgeometry.cc"

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

typedef Dune::GridSelector::GridType Grid;

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
  checkGeometry( geogrid.leafView() );
  for( int i = 0; i <= geogrid.maxLevel(); ++i )
    checkGeometry( geogrid.levelView( i ) );

  std::cerr << "Checking geometry in father..." << std::endl;
  checkGeometryInFather( geogrid );
  std::cerr << "Checking intersections..." << std::endl;
  checkIntersectionIterator( geogrid, !EnableLevelIntersectionIteratorCheck< Grid >::v );

  checkIterators( geogrid.leafView() );
  for( int i = 0; i <= geogrid.maxLevel(); ++i )
    checkIterators( geogrid.levelView( i ) );

  checkPartitionType( geogrid.leafView() );
  for( int i = 0; i <= geogrid.maxLevel(); ++i )
    checkPartitionType( geogrid.levelView( i ) );

  std::cerr << "Checking geometry lifetime..." << std::endl;
  checkGeometryLifetime( geogrid.leafView() );

  std::cerr << "Checking communication..." << std::endl;
  checkCommunication( geogrid, -1, std::cout );
  if( EnableLevelIntersectionIteratorCheck< Grid >::v )
  {
    for( int i = 0; i <= geogrid.maxLevel(); ++i )
      checkCommunication( geogrid, i, std::cout );
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
