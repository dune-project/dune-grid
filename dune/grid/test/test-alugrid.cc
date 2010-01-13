// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#define NEW_SUBENTITY_NUMBERING 1
// #undef DUNE_ENABLE_OLD_NUMBERING

#include <iostream>
#include <sstream>
#include <string>

#include <dune/common/mpihelper.hh>

#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfwriter.hh>

#include "gridcheck.cc"

#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkcommunicate.cc"
#include "checktwists.cc"

#include <dune/grid/io/visual/grapegriddisplay.hh>

using namespace Dune;

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
        grid.mark(1, *it );
        if (nr>size*0.8) break;
      }
    }
    grid.adapt();
    grid.postAdapt();
    grid.loadBalance();
  }
}

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
      grid.mark(1, *it);
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

template <class GridView>
void writeFile( const GridView& gridView )
{
  DGFWriter< GridView > writer( gridView );
  writer.write( "dump.dgf" );
}

template <class GridType>
void checkALUSerial(GridType & grid, int mxl = 2, const bool display = false)
{
  //writeFile( grid.leafView() );

  if( display )
  {
    GrapeGridDisplay< GridType > grape( grid );
    grape.display();
  }

  std::cout << "  CHECKING: Refined" << std::endl;
  // be careful, each global refine create 8 x maxlevel elements
  std::cout << "  CHECKING: Macro" << std::endl;
  gridcheck(grid);

  // only check twists for simplex grids
  // const bool checkTwist = grid.geomTypes(0)[0].isSimplex();

  // if( checkTwist )
  //  checkTwists( grid.leafView(), NoMapTwist() );

  for(int i=0; i<mxl; i++)
  {
    grid.globalRefine( DGFGridInfo<GridType> :: refineStepsForHalf() );
    std::cout << "  CHECKING: Refined" << std::endl;
    gridcheck(grid);
    // if( checkTwist )
    //  checkTwists( grid.leafView(), NoMapTwist() );

    if( display )
    {
      GrapeGridDisplay< GridType > grape( grid );
      grape.display();
    }
  }

  // check also non-conform grids
  makeNonConfGrid(grid,0,1);

  if( display )
  {
    GrapeGridDisplay< GridType > grape( grid );
    grape.display();
  }

  std::cout << "  CHECKING: non-conform" << std::endl;
  gridcheck(grid);
  std::cout << "  CHECKING: twists " << std::endl;
  // if( checkTwist )
  //  checkTwists( grid.leafView(), NoMapTwist() );

  // check the method geometryInFather()
  std::cout << "  CHECKING: geometry in father" << std::endl;
  checkGeometryInFather(grid);
  // check the intersection iterator and the geometries it returns
  std::cout << "  CHECKING: intersections" << std::endl;
  checkIntersectionIterator(grid);

  // some checks for assignment of iterators
  checkIteratorAssignment(grid);

  checkLevelIndexNonConform(grid);

  std::cout << std::endl << std::endl;
}

template <class GridType>
void checkALUParallel(GridType & grid, int gref, int mxl = 3)
{
#if HAVE_MPI
  makeNonConfGrid(grid,gref,mxl);

  // -1 stands for leaf check
  checkCommunication(grid, -1, std::cout);

  for(int l=0; l<= mxl; ++l)
    checkCommunication(grid, l , Dune::dvverb);
#endif
}


int main (int argc , char **argv) {

  // this method calls MPI_Init, if MPI is enabled
  MPIHelper & mpihelper = MPIHelper::instance(argc,argv);
  int myrank = mpihelper.rank();
  int mysize = mpihelper.size();

  try {
    /* use grid-file appropriate for dimensions */

    std::string key;
    bool initialize = true ;
    if( argc >= 2 )
    {
      key = argv[1];
      initialize = false;
    }
    else
    {
      std::cout << "usage:" << argv[0] << " <2d|2dsimp|2dconf|3d|3dsimp|3dcube> <display>" << std::endl;
    }

    const bool display = (argc > 2);

    bool testALU2dSimplex = initialize ;
    bool testALU2dConform = initialize ;
    bool testALU3dSimplex = initialize ;
    bool testALU3dCube    = initialize ;

    if( key == "2d" )
    {
      testALU2dSimplex = true ;
      testALU2dConform = true ;
    }

    if( key == "2dsimp" ) testALU2dSimplex = true ;
    if( key == "2dconf" ) testALU2dConform = true ;

    if( key == "3d" )
    {
      testALU3dSimplex = true ;
      testALU3dCube    = true ;
    }
    if( key == "3dsimp" ) testALU3dSimplex = true ;
    if( key == "3dcube" ) testALU3dCube    = true ;

    // extra-environment to check destruction
    {
      factorEpsilon = 5.e+5;
      // check empty grid

#ifndef NO_3D
      if (myrank == 0 && (testALU3dCube || testALU3dSimplex) )
        std::cout << "Check empty grids" << std::endl;

      if( testALU3dCube )
      {
        ALUCubeGrid<3,3> grid;
        checkALUSerial(grid);
      }

      if( testALU3dSimplex )
      {
        ALUSimplexGrid<3,3> grid;
        checkALUSerial(grid);
      }
#endif

#ifndef NO2D
      // check non-confrom ALUGrid for 2d
      if( testALU2dSimplex )
      {
        typedef ALUSimplexGrid<2,2> GridType;
        std::string filename("simplex-testgrid-2-2.dgf");
        std::cout << "READING from " << filename << std::endl;
        GridPtr<GridType> gridPtr(filename);
        checkALUSerial(*gridPtr, 2, display);

        //CircleBoundaryProjection<2> bndPrj;
        //GridType grid("alu2d.triangle", &bndPrj );
        //checkALUSerial(grid,2);
      }

      // check confrom ALUGrid for 2d
      if( testALU2dConform )
      {
        typedef ALUConformGrid<2,2> GridType;
        std::string filename("simplex-testgrid-2-2.dgf");
        GridPtr<GridType> gridPtr(filename);
        checkALUSerial(*gridPtr, 2, display);

        //CircleBoundaryProjection<2> bndPrj;
        //GridType grid("alu2d.triangle", &bndPrj );
        //checkALUSerial(grid,2);
      }
#endif

#ifndef NO_3D
      if( testALU3dCube )
      {
        std::string filename( "simplex-testgrid-3-3.dgf" );

        typedef ALUCubeGrid<3,3> GridType;
        GridPtr<GridType> gridPtr(filename);
        GridType & grid = *gridPtr;

        {
          std::cout << "Check serial grid" << std::endl;
          checkALUSerial(grid,
                         (mysize == 1) ? 1 : 0,
                         (mysize == 1) ? display : false);
        }

        // perform parallel check only when more then one proc
        if(mysize > 1)
        {
          if (myrank == 0) std::cout << "Check conform grid" << std::endl;
          checkALUParallel(grid,1,0);
          if (myrank == 0) std::cout << "Check non-conform grid" << std::endl;
          checkALUParallel(grid,0,2);
        }
      }

      if( testALU3dSimplex )
      {
        std::string filename( "simplex-testgrid-3-3.dgf" );

        typedef ALUSimplexGrid<3,3> GridType;
        GridPtr<GridType> gridPtr(filename);
        GridType & grid = *gridPtr;

        {
          std::cout << "Check serial grid" << std::endl;
          checkALUSerial(grid,
                         (mysize == 1) ? 1 : 0,
                         (mysize == 1) ? display : false);
        }

        // perform parallel check only when more then one proc
        if(mysize > 1)
        {
          if (myrank == 0) std::cout << "Check conform grid" << std::endl;
          checkALUParallel(grid,0,0);  //1,3
          if (myrank == 0) std::cout << "Check non-conform grid" << std::endl;
          checkALUParallel(grid,0,2);  //1,3
        }
      }
#endif
    };

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
