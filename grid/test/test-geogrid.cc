// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#ifdef COORDFUNCTION

#include <dune/grid/geogrid.hh>
#include <dune/grid/geogrid/cachedcoordfunction.hh>

#if HAVE_DUNE_PSG
#include <dune/grid/io/file/dgfparser/dgfpsggridtype.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include "functions.hh"

#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkcommunicate.cc>
#include <dune/grid/test/checkgeometryinfather.cc>
#include <dune/grid/test/checkintersectionit.cc>

using namespace Dune;

typedef COORDFUNCTION AnalyticalCoordFunctionType;

#if CACHECOORDFUNCTION
typedef CachedCoordFunction< GridType, AnalyticalCoordFunctionType > CoordFunctionType;
#else
typedef AnalyticalCoordFunctionType CoordFunctionType;
#endif

typedef GeometryGrid< GridType, CoordFunctionType > GeometryGridType;


int main ( int argc, char **argv )
try
{
  //MPIHelper &mpi = MPIHelper :: instance( argc, argv );
  MPIHelper :: instance( argc, argv );

  if( argc < 2 )
  {
    std :: cerr << "Usage: " << argv[ 0 ] << " <dgffile>" << std::endl;
    return 1;
  }

  // create Grid from DGF parser
  GridPtr< GridType > grid( argv[ 1 ] );

#if CACHECOORDFUNCTION
  AnalyticalCoordFunctionType analyticalFunction;
  CoordFunctionType coordFunction( *grid, analyticalFunction );
#else
  CoordFunctionType coordFunction;
#endif
  GeometryGridType geogrid( *grid, coordFunction );
  //GridType &geogrid = *grid;

  geogrid.globalRefine( 1 );
  geogrid.loadBalance();

  std :: cerr << "Checking grid..." << std :: endl;
  gridcheck( geogrid );
  std :: cerr << "Checking geometry in father..." << std :: endl;
  checkGeometryInFather( geogrid );
  std :: cerr << "Checking intersections..." << std :: endl;
  checkIntersectionIterator( geogrid );

  std :: cerr << "Checking communication..." << std :: endl;
  checkCommunication( geogrid, -1, std :: cout );
  for( int i = 0; i <= geogrid.maxLevel(); ++i )
    checkCommunication( geogrid, i, std :: cout );

  return 0;
}
catch( const Exception &e )
{
  std :: cerr << e << std :: endl;
  return 1;
}
catch( ... )
{
  std :: cerr << "Unknown exception raised." << std :: endl;
  return 1;
}

#else
#error "COORDFUNCTION not defined (e.g., Helix, Circle; see functions.hh)"
#endif
