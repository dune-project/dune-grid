// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$
#ifndef DUNE_GRIDCHECK_CC
#define DUNE_GRIDCHECK_CC

/** \file
    \brief Implements a generic grid check


   \todo check return types

 */

#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/capabilities.hh>

#include "staticcheck.hh"
#include "checkindexset.cc"
#include "checkgeometry.cc"
#include "checkentityseed.cc"


#include <limits>

// machine epsilon is multiplied by this factor
static double factorEpsilon = 1.e8;

class CheckError : public Dune::Exception {};


// check
// Entity::geometry()[c] == Entity::entity<dim>.geometry()[0]
// for codim=cd
template <int cd, class Grid, class Entity, bool doCheck>
struct subIndexCheck
{
  subIndexCheck ( const Grid &g, const Entity &e )
  {
    {
      typedef typename Grid::template Codim< Entity::codimension >::EntityPointer EntityPointer;
      typedef typename Grid::template Codim< Entity::codimension >::EntitySeed EntitySeed;
      #ifndef NDEBUG
      EntitySeed seed = e.seed();

      EntityPointer ep1 ( e );
      // regain entity pointer and check equality
      EntityPointer ep2 = g.entityPointer( seed );
      #endif
      assert( ep1 == ep2 );
    }

    typedef typename Grid::template Codim< cd >::EntityPointer EntityPointer;
    const int imax = e.template count<cd>();
    for( int i = 0; i < imax; ++i )
    {
      // check construction of entity pointers
      EntityPointer ep( *(e.template subEntity< cd >( i ) ) );
      assert( ep == e.template subEntity< cd >( i ) );

      #ifndef NDEBUG
      typedef typename Grid::template Codim< cd >::EntitySeed EntitySeed;
      EntitySeed seed = ep->seed();

      // regain entity pointer and check equality
      EntityPointer ep2 = g.entityPointer( seed );
      #endif
      assert( ep == ep2 );

      const typename Grid::LevelGridView &levelGridView = g.levelView(e.level());

      if( !levelGridView.contains( *ep ) )
      {
        std::cerr << "Error: Level grid view does not contain all subentities." << std::endl;
        assert( false );
      }

      const typename Grid::LevelIndexSet &levelIndexSet = g.levelIndexSet( e.level() );

      if( !levelIndexSet.contains( *ep ) )
      {
        std::cerr << "Error: Level index set does not contain all subentities." << std::endl;
        assert( false );
      }
      if( levelIndexSet.index( *ep ) != levelIndexSet.subIndex( e, i, cd ) )
      {
        int id_e = levelIndexSet.index( e );
        int id_e_i = levelIndexSet.index( *ep );
        int subid_e_i = levelIndexSet.subIndex( e, i, cd );
        std::cerr << "Error: levelIndexSet.index( *(e.template subEntity< cd >( i ) ) ) "
                  << "!= levelIndexSet.subIndex( e, i, cd )  "
                  << "[with cd=" << cd << ", i=" << i << "]" << std::endl;
        std::cerr << "       ... index( e ) = " << id_e << std::endl;
        std::cerr << "       ... index( e.subEntity< cd >( i ) ) = " << id_e_i << std::endl;
        std::cerr << "       ... subIndex( e, i, cd ) = " << subid_e_i << std::endl;
        assert( false );
      }

    }

    subIndexCheck< cd-1, Grid, Entity, Dune::Capabilities::hasEntity< Grid, cd-1 >::v > sick( g, e );
  }
};


// end recursion of subIndexCheck
template <class Grid, class Entity, bool doCheck>
struct subIndexCheck<-1, Grid, Entity, doCheck>
{
  subIndexCheck (const Grid & g, const Entity & e)
  {
    return;
  }
};
// do nothing if doCheck==false
template <int cd, class Grid, class Entity>
struct subIndexCheck<cd, Grid, Entity, false>
{
  subIndexCheck (const Grid & g, const Entity & e)
  {
    subIndexCheck<cd-1,Grid,Entity,
        Dune::Capabilities::hasEntity<Grid,cd-1>::v> sick(g,e);
  }
};
template <class Grid, class Entity>
struct subIndexCheck<-1, Grid, Entity, false>
{
  subIndexCheck (const Grid & g, const Entity & e)
  {
    return;
  }
};

// name says all
template <class Grid>
void zeroEntityConsistency (Grid &g)
{
  typedef typename Grid::ctype ctype;
  const int dimGrid = Grid::dimension;
  // const int dimWorld = Grid::dimensionworld;

  typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
  typedef typename Grid::template Codim<0>::Geometry Geometry;
  typedef typename Grid::template Codim<0>::Entity Entity;
  typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;

  const typename Grid::LevelIndexSet &levelIndexSet = g.levelIndexSet( g.maxLevel() );
  const typename Grid::LeafIndexSet &leafIndexSet = g.leafIndexSet();
  const typename Grid::LocalIdSet &localIdSet = g.localIdSet();
  const typename Grid::GlobalIdSet &globalIdSet = g.globalIdSet();

  LevelIterator it = g.template lbegin<0>(g.maxLevel());
  const LevelIterator endit = g.template lend<0>(g.maxLevel());

  for (; it!=endit; ++it)
  {
    // check construction from entity
    EntityPointer ep( *it ) ;
    assert( ep == it );

    // Entity::subEntity<0>(0) == Entity
    EntityPointer subEntity = it->template subEntity< 0 >( 0 );
    if( levelIndexSet.index( *subEntity ) != levelIndexSet.index( *it ) )
    {
      std::cerr << "Error: Level index for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }
    if( leafIndexSet.index( *subEntity ) != leafIndexSet.index( *it ) )
    {
      std::cerr << "Error: Leaf index for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }
    if( globalIdSet.id( *subEntity ) != globalIdSet.id( *it ) )
    {
      std::cerr << "Error: Global id for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }
    if( localIdSet.id( *subEntity ) != localIdSet.id( *it ) )
    {
      std::cerr << "Error: Local id for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }

    if( subEntity->level() != it->level() )
    {
      std::cerr << "Error: Level for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }

    // Entity::count<dim>() == Entity::geometry().corners();
    // Entity::geometry()[c] == Entity::entity<dim>.geometry()[0];
    const int numCorners = it->template count< dimGrid >();
    if( numCorners != it->geometry().corners() )
    {
      std::cerr << "Error: Entity::count< dimGrid >() != Entity::geometry().corners()." << std::endl;
      assert( false );
    }
    for( int c = 0; c < numCorners; ++c )
    {
      typename Geometry::GlobalCoordinate c1( it->geometry().corner( c ) );
      typename Geometry::GlobalCoordinate c2( it->template subEntity< dimGrid >( c )->geometry().corner( 0 ) );

      if( (c2-c1).two_norm() > 10*std::numeric_limits< ctype >::epsilon() )
      {
        std::cerr << "Error: | geometry.corner( c ) - subEntity< dimGrid >( c ) | = | "
                  << c1 << " - " << c2 << " | = " << (c2-c1).two_norm()
                  << " != 0  [with c = " << c << " ]" << std::endl;
        assert( false );
      }
    }

    // check that corner coordinates are distinct
    const int corners= it->geometry().corners();
    for (int c=0; c<corners; ++c)
    {
      for(int d=c+1; d<corners; ++d)
      {
        const typename Geometry::GlobalCoordinate c1 = it->geometry().corner( c );
        const typename Geometry::GlobalCoordinate c2 = it->geometry().corner( d );
        if( (c2-c1).two_norm() < 10 * std::numeric_limits<typename Grid::ctype>::epsilon() )
          DUNE_THROW( CheckError, "geometry.corner( " << c << " ) == geometry.corner( " << d << " )" );
      }
    }

    // check hasSingleGeometryType capability
    const bool hasSingleGeomType = Dune :: Capabilities :: hasSingleGeometryType< Grid > :: v ;
    if( hasSingleGeomType )
    {
      // check that geometry type generated from constant topo Id is the same as of the codim 0 entity
      const Dune::GeometryType constantType( Dune :: Capabilities :: hasSingleGeometryType< Grid > :: topologyId , dimGrid );
      if( it->type() != constantType )
      {
        DUNE_THROW(Dune::InvalidStateException,"it->type() " << it->type() << " differs from singleGeometryType " << constantType);
      }
    }

    const bool isCartesian = Dune :: Capabilities :: isCartesian< Grid > :: v ;
    // check geometry for Cartesian grids
    if( isCartesian )
    {
      if( ! it->geometry().affine() )
        DUNE_THROW( Dune::GridError, "Geometry is not affine, although isCartesian is true");
    }

    subIndexCheck<Grid::dimension, Grid, Entity, true> sick(g,*it);
  }
}

/*
 * search the LevelIterator for each IntersectionIterator
 */
template <class Grid>
void assertNeighbor (Grid &g)
{
  typedef typename Grid::LevelGridView GridView;
  typedef typename GridView::IntersectionIterator IntersectionIterator;

  typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
  enum { dim = Grid::dimension };
  //typedef typename Grid::ctype ct DUNE_UNUSED;

  typedef typename Grid::GlobalIdSet GlobalIdSet;
  const GlobalIdSet & globalid = g.globalIdSet();

  typedef typename Grid::template Codim< 0 >::Entity Entity;
  typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;

  GridView gridView = g.levelView( 0 );

  LevelIterator e = g.template lbegin<0>(0);
  const LevelIterator eend = g.template lend<0>(0);

  LevelIterator next = e;
  if (next != eend)
  {
    ++next;

    if( !EnableLevelIntersectionIteratorCheck< Grid >::v )
    {
      static int called = 0;
      if( called++ == 0 )
        std::cerr << "Warning: Skipping level neighbor check." << std::endl;
      return;
    }

    // small check if LevenIntersectionIterators
    // from work reassigned EntityPointers
    // after creation of LevelIterator on different level
    if (g.maxLevel()>0)
    {
      EntityPointer p( g.template lbegin<0>(0) );
      p = g.template lbegin<0>(1);
      LevelIterator it = g.template lbegin<0>(0);
      ++it;
#if !DISABLE_DEPRECATED_METHOD_CHECK
      p->ilevelbegin();
#endif // #if !DISABLE_DEPRECATED_METHOD_CHECK
    }

    // iterate over level
    for (; e != eend; ++e)
    {
      const Entity &entity = *e;

      const bool isGhost = (entity.partitionType() == Dune::GhostEntity);

      // call global id
      globalid.id( entity );

      if( !entity.isLeaf() )
      {
#if !DISABLE_DEPRECATED_METHOD_CHECK
        if( entity.ileafbegin() != entity.ileafend() )
        {
          DUNE_THROW(CheckError, "On non-leaf entities ileafbegin should be equal to ileafend!");
        }
#endif // #if !DISABLE_DEPRECATED_METHOD_CHECK
      }

      const int numFaces = entity.template count< 1 >();
      // flag vector for elements faces
      std::vector< bool > visited( numFaces, false );

      // loop over intersections
      const IntersectionIterator endit = gridView.iend( entity );
      IntersectionIterator it = gridView.ibegin( entity );

      if( it == endit )
      {
        // ALU has no intersection iterators on ghosts, so make this case non-fatal
        if( !isGhost )
        {
          std::cerr << "Error: Entity without intersections encountered." << std::endl;
          assert( false );
        }
        else
          std::cerr << "Warning: Ghost Entity without intersections encountered." << std :: endl;
      }

      // for all intersections
      for(; it != endit; ++it)
      {
        // numbering
        const int numberInSelf = it->indexInInside();
        if( (numberInSelf < 0) || (numberInSelf >= numFaces) )
        {
          std :: cout << "Error: Invalid numberInSelf: " << numberInSelf
                      << " (should be between 0 and " << (numFaces-1) << ")"
                      << std :: endl;
          assert( false );
        }
        else
        {
          // mark face visited
          visited[ numberInSelf ] = true;
        }

        // index / id of boundary segment
        if( it->boundary() )
        {
          it->boundarySegmentIndex();
        }

        // check id
        assert( globalid.id(*(it->inside())) == globalid.id( entity ) );


        // geometry
        it->geometryInInside();
        it->geometry();

        // normal vectors
        typename IntersectionIterator::Intersection::LocalCoordinate v(0);

        it->outerNormal(v);
        it->integrationOuterNormal(v);
        it->unitOuterNormal(v);

        if( isGhost && !it->neighbor() )
        {
          std::cerr << "Error: On ghosts, all intersections must have neighbors." << std::endl;
          assert( false );
        }

        if( it->neighbor() )
        {
          const EntityPointer outsidePtr = it->outside();
          const Entity &outside = *outsidePtr;

          if( isGhost && (outside.partitionType() == Dune::GhostEntity) )
          {
            std::cerr << "Error: Intersections between ghosts shall not be returned." << std::endl;
            assert( false );
          }

          // geometry
          it->geometryInOutside();

          // numbering
          const int indexInOutside = it->indexInOutside();
          const int numFaces = outside.template count< 1 >();
          if( (indexInOutside < 0) || (indexInOutside >= numFaces) )
          {
            std :: cout << "Error: Invalid indexInOutside: " << indexInOutside
                        << " (should be between 0 and " << (numFaces-1) << ")"
                        << std :: endl;
            assert( false );
          }

          // search neighbouring cell
          if( outsidePtr->partitionType() == Dune::InteriorEntity )
          {
            assert( globalid.id( outside ) != globalid.id( entity ) );
            const Dune::PartitionIteratorType pitype = Dune::InteriorBorder_Partition;

            typedef typename Grid::template Codim< 0 >
            ::template Partition< pitype >::LevelIterator
            LevelIterator;

            const int level = entity.level();
            bool foundNeighbor = false;
            LevelIterator nit = g.template lbegin< 0, pitype >( level );
            const LevelIterator nend = g.template lend< 0, pitype > ( level );
            for( ; nit != nend; ++nit )
            {
              if( nit->partitionType() != Dune::InteriorEntity )
              {
                std::cerr << "Error: LevelIterator for InteriorBorder_Partition "
                          << "stops on non-interior entity." << std :: endl;
                assert( false );
              }

              if( nit != outsidePtr )
                assert( globalid.id( outside ) != globalid.id( *nit ) );
              else
                foundNeighbor = true;
            }
            if( !foundNeighbor )
            {
              std :: cerr << "Error: Interior neighbor returned by "
                          << "LevelIntersectionIterator not found on that level."
                          << std :: endl;
              assert( false );
            }
          }
        }
      }

      // check that all faces were visited
      // note: This check is wrong on ghosts, where only intersections with
      //       the domain are allowed.
      if( !isGhost )
      {
        for(size_t i=0; i<visited.size(); i++) assert(visited[i] == true);
      }
    }
  }
}

template <class GridType, bool c>
struct CheckMark
{
  template <class IteratorType>
  static void check(GridType & grid, IteratorType & it)
  {
    // last marker is 0, so the grid is not changed after this check
    const int refCount[4] = {1,0,-1,0};
    for(int k=0; k<4; ++k)
    {
      // mark entity
      bool marked = grid.mark( refCount[k] , *it );
      // if element was marked, check that the marker was set correctly
      if(marked)
      {
        // now getMark should return the mark we just set, otherwise error
        if( grid.getMark( *it ) < refCount[k] )
          DUNE_THROW(CheckError,"mark/getMark method not working correctly!");
      }
    }
  }
};

template <class GridType>
struct CheckMark<GridType,false>
{
  template <class IteratorType>
  static void check(const GridType & grid, IteratorType & it )
  {}
};

/*
 * Iterate over the grid und do some runtime checks
 */

template <bool checkMark , class Grid>
void iterate(Grid &g)
{
  typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
  //typedef typename Grid::template Codim<0>::EntityPointer EntityPointer DUNE_UNUSED;
  //typedef typename Grid::HierarchicIterator HierarchicIterator DUNE_UNUSED;
  typedef typename Grid::template Codim<0>::Geometry Geometry;
  int l = g.maxLevel();
  LevelIterator it = g.template lbegin<0>(l);
  const LevelIterator endit = g.template lend<0>(l);

  typename Geometry::LocalCoordinate origin(1);
  typename Geometry::LocalCoordinate result;

  for (; it != endit; ++it)
  {
    LevelIterator l1 = it;
    LevelIterator l2 = l1; ++l1;
    assert( (l2 == it) && (it == l2) );
    assert( (l1 != it) && (it != l1) );
    ++l2;
    assert( (l1 == l2) && (l2 == l1) );

    const Geometry &geo = it->geometry();

    result = geo.local( geo.global( origin ) );
    typename Grid::ctype error = (result-origin).two_norm();
    if(error >= factorEpsilon * std::numeric_limits<typename Grid::ctype>::epsilon())
    {
      DUNE_THROW(CheckError, "|| geom.local(geom.global(" << origin
                                                          << ")) - origin || != 0 ( || " << result << " - origin || ) = " << error);
    };
    geo.integrationElement( origin );

    if((int)Geometry::coorddimension == (int)Geometry::mydimension)
      geo.jacobianInverseTransposed( origin );

    if( geo.type() != it->type() )
    {
      std::cerr << "inconsistent geometry type: entity.type() = " << it->type()
                << ", geometry.type() = " << geo.type() << "." << std::endl;
      assert( false );
    }
    geo.corners();
    geo.corner( 0 );
  }

  typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
  LeafIterator lit = g.template leafbegin<0>();
  const LeafIterator lend = g.template leafend<0>();

  // if empty grid, do nothing
  if(lit == lend) return;

  for (; lit != lend; ++lit)
  {
    LeafIterator l1 = lit;
    LeafIterator l2 = l1; ++l1;
    assert(l2 == lit);
    assert(l1 != lit);
    ++l2;
    assert(l1 == l2);

    // leaf check
    if( !lit->isLeaf() )
      DUNE_THROW(CheckError,"LeafIterator gives non-leaf entity!");

    // check adaptation mark for leaf entity mark
    CheckMark<Grid,checkMark>::check(g,lit);

    result = lit->geometry().local(lit->geometry().global(origin));
    typename Grid::ctype error = (result-origin).two_norm();
    if(error >= factorEpsilon * std::numeric_limits<typename Grid::ctype>::epsilon())
    {
      DUNE_THROW(CheckError, "|| geom.local(geom.global(" << origin
                                                          << ")) - origin || != 0 ( || " << result << " - origin || ) = " << error);
    };
    lit->geometry().integrationElement(origin);
    if((int)Geometry::coorddimension == (int)Geometry::mydimension)
      lit->geometry().jacobianInverseTransposed(origin);

    lit->geometry().type();
    lit->type();
    lit->geometry().corners();
    lit->geometry().corner( 0 );
  }

}

template <class Grid>
void iteratorEquals (Grid &g)
{
  typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
  typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
  typedef typename Grid::HierarchicIterator HierarchicIterator;
  typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;

  typedef typename Grid::LeafGridView LeafGridView;
  typedef typename LeafGridView::IntersectionIterator LeafIntersectionIterator;

  LeafGridView leafView = g.leafView();

  // assignment tests
  LevelIterator l1 = g.template lbegin<0>(0);
  LevelIterator l2 = g.template lbegin<0>(0);
  LeafIterator L1 = g.template leafbegin<0>();
  LeafIterator L2 = g.template leafbegin<0>();

  // if grid empty, leave
  if (l1 == g.template lend<0>(0))
    return;

  // check '==' consistency
  EntityPointer a( g.template levelView<Dune::All_Partition>(0).template begin<0>() );
  EntityPointer i( g.template levelView<Dune::Interior_Partition>(0).template begin<0>() );

  assert(
    (g.levelIndexSet(0).index(*a) != g.levelIndexSet(0).index(*i)) // index equal
    == (a != i) // entitypointer equal
    );
  assert(
    (g.levelIndexSet(0).index(*a) == g.levelIndexSet(0).index(*i)) // index equal
    == (a == i) // entitypointer equal
    );

  HierarchicIterator h1 = l1->hbegin(99);
  HierarchicIterator h2 = l2->hbegin(99);
  LeafIntersectionIterator i1 = leafView.ibegin( *l1 );
  LeafIntersectionIterator i2 = leafView.ibegin( *l2 );
  EntityPointer e1( l1 );
  EntityPointer e2( h2 );

  // assign
  l1 = l2;
  L1 = L2;
  h1 = h2;
  i1 = i2;
  e1 = e2;

  // equals
  #define TestEquals(i) { \
    { bool DUNE_UNUSED tmp = (i == e2); }                         \
    { bool DUNE_UNUSED tmp = (i == l2); }                         \
    { bool DUNE_UNUSED tmp = (i == h2); }                         \
    { bool DUNE_UNUSED tmp = (i == L2); }                         \
    if (i2 != leafView.iend( *l2 )) bool DUNE_UNUSED tmp = (i == i2->inside()); \
    if (i2 != leafView.iend( *l2 ) && i2->neighbor()) bool DUNE_UNUSED tmp = (i == i2->outside()); \
}
  TestEquals(e1);
  TestEquals(l1);
  TestEquals(h1);
  TestEquals(L1);
  if (i1 != leafView.iend( *l2 )) TestEquals(i1->inside());
  if (i1 != leafView.iend( *l2 ) && i1->neighbor()) TestEquals(i1->outside());
}


template< class Grid, class Entity, class Intersection >
void checkBoundarySegmentIndexProlongation ( const Grid &grid, const Entity &entity, const Intersection &intersection )
{
  typedef typename Entity::HierarchicIterator HierarchicIterator;
  typedef typename Entity::LocalGeometry GeometryInFather;
  typedef typename Grid::LevelGridView GridView;
  typedef typename GridView::IntersectionIterator IntersectionIterator;

  const typename Intersection::LocalGeometry &geoInInside = intersection.geometryInInside();

  const GridView gridView = grid.levelView( entity.level()+1 );
  const HierarchicIterator hend = entity.hend( entity.level()+1 );
  for( HierarchicIterator hit = entity.hbegin( entity.level()+1 ); hit != hend; ++hit )
  {
    GeometryInFather geoInFather = hit->geometryInFather();
    const Dune::ReferenceElement< typename Grid::ctype, Grid::dimension > &refElement
      = Dune::ReferenceElements< typename Grid::ctype, Grid::dimension >::general( geoInFather.type() );

    const IntersectionIterator iend = gridView.iend( *hit );
    for( IntersectionIterator iit = gridView.ibegin( *hit ); iit != iend; ++iit )
    {
      if( !iit->boundary() )
        continue;
      const typename GeometryInFather::GlobalCoordinate x
        = geoInFather.global( refElement.position( iit->indexInInside(), 1 ) );
      const typename GeometryInFather::GlobalCoordinate y = geoInInside.global( geoInInside.local( x ) );
      if( (x-y).two_norm2() > 1e-10 )
        continue;
      if( iit->boundarySegmentIndex() != intersection.boundarySegmentIndex() )
      {
        std::cerr << "Error: Boundary segment index " << intersection.boundarySegmentIndex()
                  << " not prolonged down (level " << entity.level() << " to " << hit->level() << ")."
                  << std::endl;
      }
      if( !hit->isLeaf() )
        checkBoundarySegmentIndexProlongation( grid, *hit, *iit );
    }
  }
}


// check for consistency of hasFather(), father(), and father().level()
template< class Grid >
void checkFatherLevel ( Grid &grid )
{
  for(int level=0; level<=grid.maxLevel(); ++level)
  {
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::template Codim<0>::Iterator Iterator;
    typedef typename GridView::template Codim<0>::EntityPointer EntityPointer;

    GridView gv = grid.levelView(level);
    Iterator it = gv.template begin<0>();
    Iterator end = gv.template end<0>();
    for(; it!=end; ++it)
    {
      if (level==0)
      {
        if (it->hasFather())
        {
          std::cerr << "Error: hasFather() is true for element " << gv.indexSet().index(*it)
                    << " on level 0"
                    << std::endl;
          assert(false);
        }
      }
      else
      {
        if (it->hasFather())
        {
          // check level of father for newly created entitypointer
          {
            EntityPointer f = it->father();
            if (f.level() != level-1)
            {
              std::cerr << "Error: father().level()=" << f.level()
                        << " for element on level " << level << std::endl;
              assert(false);
            }
            if (f->level() != level-1)
            {
              std::cerr << "Error: father()->level()=" << f.level()
                        << " for element on level " << level << std::endl;
              assert(false);
            }
          }
          // check level of father after reassignment
          {
            EntityPointer f(*it);
            f = it->father();
            if (f.level() != level-1)
            {
              std::cerr << "Error: father().level()=" << f.level()
                        << " for element on level " << level << " with reassigned father pointer" << std::endl;
              assert(false);
            }
            if (f->level() != level-1)
            {
              std::cerr << "Error: father()->level()=" << f.level()
                        << " for element on level " << level << " with reassigned father pointer" << std::endl;
              assert(false);
            }
          }
        }
      }
    }
  }
}


template< class GridView >
void checkBoundarySegmentIndex ( const GridView &gridView )
{
  typedef typename GridView::template Codim< 0 >::Iterator Iterator;
  typedef typename GridView::IntersectionIterator IntersectionIterator;

  size_t numBoundarySegments = gridView.grid().numBoundarySegments();
  size_t countBoundarySegments = 0;
  std::vector< int > count( numBoundarySegments, 0 );
  bool error = false;

  const Iterator end = gridView.template end< 0 >();
  for( Iterator it = gridView.template begin< 0 >(); it != end; ++it )
  {
    if( !it->hasBoundaryIntersections() )
      continue;
    const IntersectionIterator iend = gridView.iend( *it );
    for( IntersectionIterator iit = gridView.ibegin( *it ); iit != iend; ++iit )
    {
      if( !iit->boundary() )
        continue;
      ++countBoundarySegments;
      const size_t index = iit->boundarySegmentIndex();
      if( index >= numBoundarySegments )
      {
        std::cerr << "Error: Boundary segment index out of bounds: (index: "
                  << index << ", size: " << numBoundarySegments << ")."
                  << std::endl;
        error = true;
      }
      else
        ++count[ index ];
      if( !it->isLeaf() )
        checkBoundarySegmentIndexProlongation( gridView.grid(), *it, *iit );
    }
  }

  if( countBoundarySegments != numBoundarySegments )
  {
    std::cerr << "Error: Wrong number of boundary segments (reported: "
              << numBoundarySegments << ", counted: "
              << countBoundarySegments << ")." << std::endl;
    error = true;
  }

  for( size_t i = 0; i < count.size(); ++i )
  {
    if( count[ i ] == 0 )
    {
      std::cerr << "Error: Boundary segment indices not consecutive (index "
                << i << " not used)." << std::endl;
      error = true;
    }
    else if( count[ i ] != 1 )
    {
      std::cerr << "Error: Boundary segment index not unique within macro grid"
                << " (index " << i << " used " << count[ i ] << " times)."
                << std::endl;
      error = true;
    }
  }

  if (error)
  {
    std::string msg("Error encountered during checkBoundarySegmentIndex.");
    if( gridView.grid().comm().size() > 1 )
      std::cerr << msg << std::endl;
    else
      DUNE_THROW(Dune::Exception, msg );
  }
}


template <class Grid>
void gridcheck (Grid &g)
{
  /*
   * first do the compile-test: this will not produce any code but
   * fails if an interface-component is missing
   */
  GridInterface<Grid>();

  enum { dim      = Grid :: dimension };
  enum { dimworld = Grid :: dimensionworld };
  typedef typename Grid  :: ctype ctype;
  typedef typename Grid  :: GridFamily GridFamily;

  /*
   * Check whether we have a nonzero number of entities in each view
   * in case the grid is non-empty
   */
  if( g.size( 0 ) > 0 )
  {
    for (int i=0; i<=dim; i++)
    {
      assert(g.size(i)>0);
      for (int j=0; j<=g.maxLevel(); j++)
        if (g.comm().size() == 0)
          assert(g.size(j,i)>0);
        else
          assert(g.size(j,i)>=0);   // in parallel this could be zero (depends on definition of maxLevel() method in parallel)
    }
  }

#if DUNE_VERBOSE_TESTS
  // print infos
  Dune::gridinfo(g, "GridInfo");
  for(int l=0; l<g.maxLevel(); l++)
    Dune::gridlevellist(g, l, "GridLevelInfo");
  Dune::gridleaflist(g, "GridLeafInfo");
#endif

  // type of GridInterface == GridDefaultImplementation
  typedef Dune::GridDefaultImplementation<dim,dimworld,ctype,GridFamily> GridIF;
  const GridIF & gridIF = g;
  // check functionality when grid is interpreted as reference to interface
  GridInterface<GridIF>::check(gridIF);
  /*
   * now the runtime-tests
   */
  const Grid & cg = g;
  iteratorEquals(g);
  iteratorEquals(cg);
  iterate<true>(g);
  iterate<false>(cg);
  zeroEntityConsistency(g);
  zeroEntityConsistency(cg);
  assertNeighbor(g);
  assertNeighbor(cg);
  checkFatherLevel(g);
  checkFatherLevel(cg);

  // check geometries of macro level and leaf level
  checkGeometry( g.levelView( 0 ) );
  checkGeometry( g.leafView() );

  // check entity seeds
  Dune::checkEntitySeed( g.leafView(), std::cerr );
  for( int level = 0; level <= g.maxLevel(); ++level )
    Dune::checkEntitySeed( g.levelView( level ), std::cerr );

  // note that for some grid this might fail
  // then un comment this test
  Dune :: checkIndexSet( g, g.leafView(), Dune :: dvverb );
  for( int level = 0; level <= g.maxLevel(); ++level )
    Dune :: checkIndexSet( g, g.levelView( level ), Dune :: dvverb, true );

  // check at least if the subId method is there
  {
    typename Grid::Traits::template Codim<0>::LevelIterator it = g.template lbegin<0>(g.maxLevel());
    typename Grid::Traits::template Codim<0>::LevelIterator end = g.template lend<0>(g.maxLevel());
    for (; it != end; ++it)
    {
      g.globalIdSet().subId(*it,0,dim);
      if(g.globalIdSet().subId(*it,0,dim) != g.globalIdSet().id(*(it->template subEntity<dim>(0))))
      {
        std::cerr << "Error: Inconsistent global subId for vertex 0 (id(subEntity)=" << g.globalIdSet().id(*(it->template subEntity<dim>(0)))
                  << ", subId=" << g.globalIdSet().subId(*it,0,dim) << ") of cell " << g.globalIdSet().id(*it) << std::endl;
        assert(false);
      }
    }
  }

  if( EnableLevelIntersectionIteratorCheck< Grid >::v )
    checkBoundarySegmentIndex( g.levelView( 0 ) );
  else if( g.maxLevel() == 0 )
    checkBoundarySegmentIndex( g.leafView() );
  else
    std::cout << "Warning: Skipping boundary segment index check (missing level intersection iterator)." << std::endl;
}

#endif // #ifndef GRIDCHECK_CC
