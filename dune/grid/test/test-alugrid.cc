// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

// #define NO_2D
// #define NO_3D

#define DISABLE_DEPRECATED_METHOD_CHECK 1

#include <iostream>
#include <sstream>
#include <string>

#include <dune/common/tupleutility.hh>
#include <dune/common/tuples.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfwriter.hh>

#include "gridcheck.cc"

#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkcommunicate.cc"
//#include "checktwists.cc"

#include <dune/grid/io/visual/grapegriddisplay.hh>

#if ALU3DGRID_PARALLEL && HAVE_MPI
#define USE_PARALLEL_TEST 1
#endif

using namespace Dune;

template<>
struct EnableLevelIntersectionIteratorCheck< Dune::ALUGrid<3, 3, simplex, conforming> >
{
  static const bool v = false;
};

template <bool leafconform, class Grid>
void checkCapabilities(const Grid& grid)
{
  dune_static_assert( Dune::Capabilities::hasSingleGeometryType< Grid > :: v == true,
                      "hasSingleGeometryType is not set correctly");
  dune_static_assert( Dune::Capabilities::isLevelwiseConforming< Grid > :: v == ! leafconform,
                      "isLevelwiseConforming is not set correctly");
  dune_static_assert( Dune::Capabilities::isLeafwiseConforming< Grid > :: v == leafconform,
                      "isLevelwiseConforming is not set correctly");
  static const bool hasEntity = Dune::Capabilities::hasEntity<Grid, 1> :: v == true;
  dune_static_assert( hasEntity,
                      "hasEntity is not set correctly");
  dune_static_assert( Dune::Capabilities::hasBackupRestoreFacilities< Grid > :: v == true,
                      "hasBackupRestoreFacilities is not set correctly");

  static const bool reallyParallel =
#if ALU3DGRID_PARALLEL
    Grid :: dimension == 3;
#else
    false ;
#endif
  dune_static_assert( Dune::Capabilities::isParallel< Grid > :: v == reallyParallel,
                      "isParallel is not set correctly");

  static const bool reallyCanCommunicate =
#if ALU3DGRID_PARALLEL
    Grid :: dimension == 3;
#else
    false ;
#endif
  static const bool canCommunicate = Dune::Capabilities::canCommunicate< Grid, 1 > :: v
                                     == reallyCanCommunicate;
  dune_static_assert( canCommunicate,
                      "canCommunicate is not set correctly");

}

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
      EntityPointerType p( it );

      assert( p.level()  == 0 );
      assert( p->level() == 0 );

      if( grid.maxLevel() > 0 )
      {
        it = grid.template lbegin<dim>(1);
        if (grid.size(1,0)>0)
        {
          p = it;
          assert( it->level() == 1 );
          assert( p.level()   == 1 );
          assert( p->level()  == 1 );
        }
      }
    }
  }
}


template <class EntityType, class LocalGeometryType>
int aluTwistCheck(const EntityType& en, const LocalGeometryType& localGeom,
                  const int face, const bool neighbor, const bool output )
{
  enum { dim = EntityType :: dimension };
  typedef typename EntityType :: ctype ctype;

  typedef FaceTopologyMapping<tetra> SimplexFaceMapping;
  typedef FaceTopologyMapping<hexa>  CubeFaceMapping;

  // get reference element
  const ReferenceElement< ctype, dim > &refElem =
    ReferenceElements< ctype, dim >::general( en.type() );

  const int vxSize = refElem.size( face, 1, dim );
  typedef  FieldVector<ctype,dim> CoordinateVectorType;

  // now calculate twist by trial and error for all possible twists
  // the calculated twist is with respect to the ALUGrid
  // reference face, see twistprovider.cc
  int twistFound = -66;
  for(int twist = -vxSize; twist<vxSize; ++twist)
  {
    bool twistOk = true;
    // now check mapping with twist
    for(int i=0; i<vxSize; ++i)
    {
      int twistedDuneIndex = -1;
      if( localGeom.type().isCube() )
      {
        twistedDuneIndex = CubeFaceMapping::twistedDuneIndex( i, twist );
      }
      else
      {
        twistedDuneIndex = SimplexFaceMapping::twistedDuneIndex( i, twist );
      }

      // get face vertices of number in self face
      int vxIdx = refElem.subEntity( face, 1 , twistedDuneIndex , dim);

      // get position in reference element of vertex i
      CoordinateVectorType refPos = refElem.position( vxIdx, dim );

      // check coordinates again
      CoordinateVectorType localPos = localGeom.corner( i );
      if( (refPos - localPos).infinity_norm() > 1e-8 )
      {
        twistOk = false;
        break;
      }
    }

    if( twistOk )
    {
      twistFound = twist;
      break ;
    }
  }

  // if no twist found, then something is wrong
  if( twistFound == -66 )
  {
    assert(false);
    DUNE_THROW(GridError,"Not matching twist found");
  }

  if( output )
  {
    std::string twistIn( (neighbor) ? "twistInNeighbor()" : "twistInSelf" );
    std::string numberIn( (neighbor) ? "indexInOutside()" : "indexInInside" );
    std::cout << "ERROR: Face "<< face << " : twist = "<< twistFound << std::endl;
    std::cout << "\nPut twist = "<< twistFound << " In TwistUtility::"<< twistIn << " for " << numberIn << " = " << face << " ! \n";
    std::cout << "******************************************\n";
  }

  return twistFound;
}

template <class GridView>
void checkALUTwists( const GridView& gridView, const bool verbose = false )
{

  typedef typename GridView :: template Codim< 0 > :: Iterator Iterator ;
  typedef typename Iterator :: Entity Entity ;
  typedef typename GridView :: IntersectionIterator IntersectionIterator ;

  const Iterator endit = gridView.template end< 0 >();
  for( Iterator it = gridView.template begin< 0 >(); it != endit ; ++it )
  {
    const Entity& entity = *it ;
    const IntersectionIterator endnit = gridView.iend( entity );
    for( IntersectionIterator nit = gridView.ibegin( entity ); nit != endnit; ++nit )
    {
      typedef typename IntersectionIterator :: Intersection Intersection;
      const Intersection& intersection = * nit ;

      // check twist of inside geometry
      const int twistInside = aluTwistCheck( entity, intersection.geometryInInside(),
                                             intersection.indexInInside(), false, verbose );
      const int twistIn = gridView.grid().getRealIntersection( intersection ).twistInInside();

      if( twistInside != twistIn )
        std::cerr << "Error: inside twists " << twistInside << " (found)  and  " << twistIn << " (given) differ" << std::endl;

      if( intersection.neighbor() )
      {
        // check twist of inside geometry
        const int twistOutside = aluTwistCheck( *intersection.outside(), intersection.geometryInOutside(),
                                                intersection.indexInOutside(), true, verbose );
        const int twistOut = gridView.grid().getRealIntersection( intersection ).twistInOutside();
        if( twistOutside != twistOut )
          std::cerr << "Error: outside twists " << twistOutside << " (found)  and  " << twistOut << " (given) differ" << std::endl;
      }
    }
  }
}

template <int codim, class GridType>
void checkIteratorCodim(GridType & grid)
{
#ifdef NO_2D
  typedef typename GridType::template Codim<codim>::
  template Partition<Dune::InteriorBorder_Partition>::LeafIterator
  IteratorInteriorBorder;

  typedef typename GridType::template Codim<codim>:: Geometry Geometry ;
  typedef typename GridType:: ctype ctype;

  /** Loop only over the interior elements, not over ghost elements. */
  IteratorInteriorBorder endIterator = grid.template leafend<codim,InteriorBorder_Partition>();
  for (IteratorInteriorBorder iter =
         grid.template leafbegin<codim,InteriorBorder_Partition>();
       iter!=endIterator; ++iter)
  {
    /** Provide geometry type of element. */
    const Geometry& geo = iter->geometry();
    if( geo.corners() > 1 )
    {
      Dune::FieldVector<ctype, GridType::dimension>
      diff( geo.corner(0) - geo.corner(1) );
      if( diff.two_norm() < 1e-8 )
      {
        std::cout << diff << " twonorm = " << diff.two_norm() << " point 0 and 1 do not differ! " << std::endl;
        assert( diff.two_norm() > 1e-8 );
      }
    }
  }
#endif
}

template <class GridType>
void checkIterators( GridType& grid )
{
  checkIteratorCodim< 0 > ( grid );
  checkIteratorCodim< 1 > ( grid );
  checkIteratorCodim< 2 > ( grid );
  checkIteratorCodim< GridType :: dimension > ( grid );
}

template <int codim, class GridType>
void checkPersistentContainerCodim(GridType & grid)
{
  typedef PersistentContainer< GridType, int > ContainerType ;

  ContainerType persistentContainer ( grid, codim ) ;
  typedef typename ContainerType :: Iterator iterator;

  // clear container
  const iterator end = persistentContainer.end();
  for( iterator it = persistentContainer.begin(); it != end; ++ it )
    *it = 0;

  typedef typename GridType::template Codim<codim>::LeafIterator Iterator;

  typedef typename GridType::template Codim<codim>:: Entity Entity;
  typedef typename GridType::template Codim<codim>:: Geometry Geometry ;
  typedef typename GridType:: ctype ctype;

  /** Loop only over the interior elements, not over ghost elements. */
  const Iterator endIterator = grid.template leafend<codim>();
  for (Iterator iter =
         grid.template leafbegin<codim>();
       iter!=endIterator; ++iter)
  {
    const Entity& entity = *iter ;
    persistentContainer[ entity ] = 1 ;
  }

  int sum = 0;
  for( iterator it = persistentContainer.begin(); it != end; ++ it )
    sum += *it;

  // the number of leaf entities should equal to what we just stored.
  if( grid.size( codim ) != sum )
    DUNE_THROW(InvalidStateException,"PersistentContainer for codim " << codim<< " gives wrong results!");
}

template <class GridType>
void checkPersistentContainer( GridType& grid )
{
  checkPersistentContainerCodim< 0 > ( grid );
  checkPersistentContainerCodim< 1 > ( grid );
  checkPersistentContainerCodim< 2 > ( grid );
  checkPersistentContainerCodim< GridType :: dimension > ( grid );
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
  const bool skipLevelIntersections = ! EnableLevelIntersectionIteratorCheck< GridType > :: v ;
  {
    GridType* gr = new GridType();
    assert( gr );
    delete gr;
  }

  //writeFile( grid.leafView() );

  if( display )
  {
    GrapeGridDisplay< GridType > grape( grid );
    grape.display();
  }

  std::cout << "  CHECKING: grid size = " << grid.size( 0 ) << std::endl;

  // be careful, each global refine create 8 x maxlevel elements
  std::cout << "  CHECKING: Macro" << std::endl;
  gridcheck(grid);
  std::cout << "  CHECKING: Macro-intersections" << std::endl;
  checkIntersectionIterator(grid, skipLevelIntersections);

  if( GridType :: dimension == 3 )
  {
    // this only works for ALUGrid 3d
    std::cout << "  CHECKING: 3d Twists " << std::endl;
    checkALUTwists( grid.leafView() );
  }

  // only check twists for simplex grids
  // const bool checkTwist = grid.geomTypes(0)[0].isSimplex();

  //if( checkTwist )
  //  checkTwists( grid.leafView(), NoMapTwist() );

  for(int i=0; i<mxl; i++)
  {
    grid.globalRefine( DGFGridInfo<GridType> :: refineStepsForHalf() );
    std::cout << "  CHECKING: Refined" << std::endl;
    gridcheck(grid);
    std::cout << "  CHECKING: intersections" << std::endl;
    checkIntersectionIterator(grid, skipLevelIntersections);
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

  // check iterators
  checkIterators( grid );

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
  checkIntersectionIterator(grid, skipLevelIntersections);

  // some checks for assignment of iterators
  checkIteratorAssignment(grid);

  // check level index sets on nonconforming grids
  checkLevelIndexNonConform(grid);

  // check life time of geometry implementation
  std::cout << "  CHECKING: geometry lifetime" << std::endl;
  checkGeometryLifetime( grid.leafView() );

  // check persistent container
  checkPersistentContainer( grid );

  std::cout << std::endl << std::endl;
}

template <class GridType>
void checkALUParallel(GridType & grid, int gref, int mxl = 3)
{
#if USE_PARALLEL_TEST
  makeNonConfGrid(grid,gref,mxl);

  // check iterators
  checkIterators( grid );

  // -1 stands for leaf check
  checkCommunication(grid, -1, std::cout);

  if( Capabilities :: isLevelwiseConforming< GridType > :: v )
  {
    for(int l=0; l<= mxl; ++l)
      checkCommunication(grid, l , Dune::dvverb);
  }
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
#ifdef ALUGRID_SURFACE_2D
      std::cout << "usage:" << argv[0] << " <2d|2dsimp|2dcube|2dconf|3d|3dsimp|3dcube> <display>" << std::endl;
#else
      std::cout << "usage:" << argv[0] << " <2d|2dsimp|2dconf|3d|3dsimp|3dcube> <display>" << std::endl;
#endif
    }

    const char *newfilename = 0;
    bool display = false;
    if( argc > 2 )
    {
      display = (std::string( argv[ 2 ] ) == "display");
      newfilename = (display ? 0 : argv[ 2 ]);
      if( newfilename && argc > 3 )
        display = true ;
    }

#ifndef NO_2D
    bool testALU2dSimplex = initialize ;
    bool testALU2dConform = initialize ;
    bool testALU2dCube    = initialize ;
    if( key == "2d" )
    {
      testALU2dSimplex = true ;
      testALU2dConform = true ;
      testALU2dCube   = true ;
    }
    if( key == "2dsimp" ) testALU2dSimplex = true ;
    if( key == "2dconf" ) testALU2dConform = true ;
    if( key == "2dcube" ) testALU2dCube    = true ;
#endif // #ifndef NO_2D

#ifndef NO_3D
    bool testALU3dSimplex = initialize ;
    bool testALU3dConform = initialize ;
    bool testALU3dCube    = initialize ;
    if( key == "3d" )
    {
      testALU3dSimplex = true ;
      testALU3dConform = true ;
      testALU3dCube    = true ;
    }
    if( key == "3dnonc" )
    {
      testALU3dSimplex = true ;
      testALU3dCube    = true ;
    }
    if( key == "3dsimp" ) testALU3dSimplex = true ;
    if( key == "3dconf" ) testALU3dConform = true ;
    if( key == "3dcube" ) testALU3dCube    = true ;
#endif // #ifndef NO_3D

    // extra-environment to check destruction
    {
      factorEpsilon = 5.e+5;
      // check empty grid

#ifndef NO_3D
      if (myrank == 0 && (testALU3dCube || testALU3dSimplex) )
        std::cout << "Check empty grids" << std::endl;

      if( testALU3dCube )
      {
        ALUGrid<3,3, cube, nonconforming > grid;
        checkALUSerial(grid);
      }

      if( testALU3dSimplex )
      {
        ALUGrid<3,3, simplex, nonconforming > grid;
        checkALUSerial(grid);
      }

      if( testALU3dConform )
      {
        ALUGrid<3,3, simplex, conforming > grid;
        checkALUSerial(grid);
      }
#endif

#ifndef NO_2D
#ifdef ALUGRID_SURFACE_2D
      // check non-conform ALUGrid for 2d
      if( testALU2dCube )
      {
        typedef ALUGrid<2,2, cube, nonconforming > GridType;
        //typedef ALUCubeGrid<2,2> GridType;
        std::string filename( DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/cube-testgrid-2-2.dgf" );
        std::cout << "READING from " << filename << std::endl;
        GridPtr<GridType> gridPtr(filename);
        checkCapabilities< false >( *gridPtr );
        checkALUSerial(*gridPtr, 2, display);

        //CircleBoundaryProjection<2> bndPrj;
        //GridType grid("alu2d.triangle", &bndPrj );
        //checkALUSerial(grid,2);

        typedef ALUGrid< 2, 3, cube, nonconforming > SurfaceGridType;
        //typedef ALUCubeGrid< 2, 3 > SurfaceGridType;
        std::string surfaceFilename( DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/cube-testgrid-2-3.dgf" );
        std::cout << "READING from '" << surfaceFilename << "'..." << std::endl;
        GridPtr< SurfaceGridType > surfaceGridPtr( surfaceFilename );
        checkCapabilities< false >( *surfaceGridPtr );
        checkALUSerial( *surfaceGridPtr, 1, display );
      }
#endif // #ifdef ALUGRID_SURFACE_2D
       // check non-conform ALUGrid for 2d
      if( testALU2dSimplex )
      {
        typedef ALUGrid< 2, 2, simplex, nonconforming > GridType;
        //typedef ALUSimplexGrid< 2, 2 > GridType;
        std::string filename(DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/simplex-testgrid-2-2.dgf");
        std::cout << "READING from " << filename << std::endl;
        GridPtr<GridType> gridPtr(filename);
        checkCapabilities< false >( *gridPtr );
        checkALUSerial(*gridPtr, 2, display);

        //CircleBoundaryProjection<2> bndPrj;
        //GridType grid("alu2d.triangle", &bndPrj );
        //checkALUSerial(grid,2);

#ifdef ALUGRID_SURFACE_2D
        typedef ALUGrid< 2, 3, simplex, nonconforming > SurfaceGridType;
        //typedef ALUSimplexGrid< 2, 3 > SurfaceGridType;
        std::string surfaceFilename( DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/simplex-testgrid-2-3.dgf" );
        std::cout << "READING from '" << surfaceFilename << "'..." << std::endl;
        GridPtr< SurfaceGridType > surfaceGridPtr( surfaceFilename );
        checkCapabilities< false >( *surfaceGridPtr );
        checkALUSerial( *surfaceGridPtr, 1, display );
#endif // #ifdef ALUGRID_SURFACE_2D
      }

      // check conform ALUGrid for 2d
      if( testALU2dConform )
      {
        typedef ALUGrid< 2, 2, simplex, conforming > GridType;
        //typedef ALUConformGrid< 2, 2 > GridType;
        std::string filename(DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/simplex-testgrid-2-2.dgf");
        GridPtr<GridType> gridPtr(filename);
        checkCapabilities< true >( *gridPtr );
        checkALUSerial(*gridPtr, 2, display);

        //CircleBoundaryProjection<2> bndPrj;
        //GridType grid("alu2d.triangle", &bndPrj );
        //checkALUSerial(grid,2);

#ifdef ALUGRID_SURFACE_2D
        typedef ALUGrid< 2, 3, simplex, conforming > SurfaceGridType;
        //typedef ALUConformGrid< 2, 3 > SurfaceGridType;
        std::string surfaceFilename( DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/simplex-testgrid-2-3.dgf" );
        std::cout << "READING from '" << surfaceFilename << "'..." << std::endl;
        GridPtr< SurfaceGridType > surfaceGridPtr( surfaceFilename );
        checkCapabilities< true >( *surfaceGridPtr );
        checkALUSerial( *surfaceGridPtr, 1, display );
#endif // #ifdef ALUGRID_SURFACE_2D
      }
#endif // #ifndef NO_2D

#ifndef NO_3D
      if( testALU3dCube )
      {
        std::string filename;
        if( newfilename )
          filename = newfilename;
        else
          filename = DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/simplex-testgrid-3-3.dgf";

        typedef ALUGrid<3,3, cube, nonconforming > GridType;
        //typedef ALUCubeGrid<3,3> GridType;
        GridPtr<GridType> gridPtr(filename);
        GridType & grid = *gridPtr;
        grid.loadBalance();

        checkCapabilities< false >( grid );

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
        std::string filename;
        if( newfilename )
          filename = newfilename;
        else
          filename = DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/simplex-testgrid-3-3.dgf";

        typedef ALUGrid<3,3, simplex, nonconforming > GridType;
        //typedef ALUSimplexGrid<3,3> GridType;
        GridPtr<GridType> gridPtr(filename);
        GridType & grid = *gridPtr;
        grid.loadBalance();
        checkCapabilities< false >( grid );

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

#ifdef ALUGRID_3D_CONFORMING_REFINEMENT
      if( testALU3dConform )
      {
        std::string filename;
        if( newfilename )
          filename = newfilename;
        else
          filename = DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/simplex-testgrid-3-3.dgf";

        typedef ALUGrid<3,3, simplex, conforming > GridType;
        GridPtr<GridType> gridPtr(filename);
        GridType & grid = *gridPtr;
        grid.loadBalance();
        checkCapabilities< true >( grid );

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
