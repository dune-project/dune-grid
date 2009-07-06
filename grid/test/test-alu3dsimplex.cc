// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <sstream>
#include <string>

#include <dune/common/mpihelper.hh>

#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include "gridcheck.cc"

#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkcommunicate.cc"
#include "checktwists.cc"
#include "check-albertareader.cc"

#include "basicunitcube.hh"

using namespace Dune;

template< ALU3dGridElementType type >
struct MapTwistAlu2Dune;

template<>
struct MapTwistAlu2Dune< tetra >
{
  int operator() ( const int twist ) const
  {
    const int map[ 6 ] = { -2, -3, -1, 0, 2, 1 };

    assert( (twist >= -3) && (twist <= 2) );
    return map[ twist+3 ];
  }
};


template <class GridType>
void makeNonConfGrid(GridType &grid,int level,int adapt) {
  int myrank = grid.comm().rank();
  grid.loadBalance();
  grid.globalRefine(level);
  grid.loadBalance();
  for (int i=0; i<adapt; i++)
  {
    if (myrank==0)
    {
      typedef typename GridType :: template Codim<0> ::
      template Partition<Interior_Partition> :: LeafIterator LeafIterator;

      LeafIterator endit = grid.template leafend<0,Interior_Partition>   ();
      int nr = 0;
      int size = grid.size(0);
      for(LeafIterator it    = grid.template leafbegin<0,Interior_Partition> ();
          it != endit ; ++it,nr++ )
      {
        grid.mark( 1, *it );
        if (nr>size*0.8) break;
      }
    }
    typename GridType::template Codim<0>::LeafIterator it = grid.template leafbegin<0>();
    grid.adapt();
    grid.postAdapt();
    grid.loadBalance();
  }
}

#if 0
template <class GridType>
void checkIteratorAssignment(GridType & grid)
{
  // check Iterator assignment
  {
    enum { dim = GridType :: dimension };
    typedef typename GridType :: template Codim<dim> :: LevelIterator
    IteratorType;

    IteratorType it = grid.template lbegin<dim>(0);
    if( grid.maxLevel() > 0 ) it = grid.template lbegin<dim>(1);
  }

  {
    enum { dim = GridType :: dimension };
    typedef typename GridType :: template Codim<dim> :: LevelIterator
    IteratorType;
    typedef typename GridType::Traits::template Codim<dim>::EntityPointer EntityPointerType;

    IteratorType it = grid.template lbegin<dim>(0);

    if( it != grid.template lend<dim>(0) )
    {
      assert( it->level() == 0 );
      EntityPointerType p = it;

      assert( p.level()  == 0 );
      assert( p->level() == 0 );

      if( grid.maxLevel() > 0 )
      {
        it = grid.template lbegin<dim>(1);
        p = it;

        assert( it->level() == 1 );
        assert( p.level()   == 1 );
        assert( p->level()  == 1 );
      }
    }
  }
}
#endif

#if 0
template <class GridType>
void checkLevelIndexNonConform(GridType & grid)
{
  typedef typename GridType :: template Codim<0> :: LeafIterator
  IteratorType;
  {
    IteratorType end = grid.template leafend<0>();
    for(IteratorType it = grid.template leafbegin<0>(); it!=end; ++it)
    {
      // call index of level index set
      grid.levelIndexSet(it->level()).index(*it);
    }
  }

  {
    IteratorType it = grid.template leafbegin<0>();
    if( it != grid.template leafend<0>() )
    {
      // mark first entity
      grid.mark( 1, *it );
    }
  }

  grid.preAdapt();
  grid.adapt();
  grid.postAdapt();

  {
    IteratorType end = grid.template leafend<0>();
    for(IteratorType it = grid.template leafbegin<0>(); it!=end; ++it)
    {
      // call index of level index set
      grid.levelIndexSet(it->level()).index(*it);
    }
  }
}
#endif

template< class GridType >
void checkALUSerial ( GridType &grid, int maxLevel = 1 )
{
  // be careful, each global refine create 8 x maxlevel elements
  if( grid.comm().rank() == 0 )
    std::cout << ">>> Checking macro grid..." << std::endl;
  gridcheck( grid );

  for( int i = 0; i < maxLevel; ++i )
  {
    grid.globalRefine( DGFGridInfo<GridType> :: refineStepsForHalf() );
    if( grid.comm().rank() == 0 )
      std::cout << ">>> Checking grid refined " << (i+1) << " times..." << std::endl;
    gridcheck( grid );
  }

  // check also non-conform grids
  makeNonConfGrid(grid,0,1);
  gridcheck(grid);

  // check the method geometryInFather()
  checkGeometryInFather(grid);

#if 1
  // check the intersection iterator and the geometries it returns
  checkIntersectionIterator(grid);
#endif
  checkTwists( grid.leafView(), MapTwistAlu2Dune< GridType::elementType >() );

#if 0
  // some checks for assignment of iterators
  checkIteratorAssignment(grid);

  checkLevelIndexNonConform(grid);
#endif
}

#if 0
#if HAVE_MPI
template <class GridType>
void checkALUParallel(GridType & grid, int gref, int mxl = 3)
{
  makeNonConfGrid(grid,gref,mxl);

  // -1 stands for leaf check
  checkCommunication(grid, -1, std::cout);

  for(int l=0; l<= mxl; ++l)
    checkCommunication(grid, l , Dune::dvverb);
}
#else
template <class GridType>
void checkALUParallel(GridType & grid, int gref, int mxl = 3)
{}
#endif
#endif

int main ( int argc, char **argv )
try {
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper &mpihelper = MPIHelper::instance( argc, argv );
  int myrank = mpihelper.rank();
  int mysize = mpihelper.size();

  // extra-environment to check destruction
  {
    typedef ALUSimplexGrid< 3, 3 > GridType;
    #if 0
    // check empty grid
    {
      if( myrank == 0 )
        std::cout << "Check empty grid" << std::endl;
      GridType grid;
      checkALUSerial( grid );
    }

    {
      if( myrank == 0 )
        std::cerr << ">>> Checking twisted unit cube..." << std::endl;
      GridFactory< GridType > factory;
      BasicUnitCube< GridType::dimension >::insertVertices( factory );
      BasicUnitCube< GridType::dimension >::insertSimplices( factory );
      GridType *grid = factory.createGrid();
      checkALUSerial( *grid );
      delete grid;
    }

    if( myrank == 0 )
      checkAlbertaReader< GridType >();
    #endif
    {
      std::string filename;
      if( argc > 1 )
        filename = argv[ 1 ];
      else if( mysize <= 2 )
        filename = "simplex-testgrid-3-3.dgf";
      else
        filename = "simplex-testgrid-3-3-large.dgf";

      GridPtr< GridType > gridPtr( filename );
      GridType &grid = *gridPtr;
      grid.loadBalance();

      if (myrank == 0)
        std::cout << "Check serial grid" << std::endl;
      checkALUSerial( grid, 1 );

#if 0
      // perform parallel check only when more then one proc
      if( mysize > 1 )
      {
        if (myrank == 0) std::cout << "Check conform grid" << std::endl;
        checkALUParallel(grid,0,0);  //1,3
        if (myrank == 0) std::cout << "Check non-conform grid" << std::endl;
        checkALUParallel(grid,0,2);  //1,3
      }
#endif
    }
  }

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
