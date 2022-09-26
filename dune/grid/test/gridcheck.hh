// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_GRIDCHECK_HH
#define DUNE_GRID_TEST_GRIDCHECK_HH

/** \file
    \brief Implements a generic grid check


   \todo check return types

 */

#include <cstdlib>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/stdstreams.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/rangegenerators.hh>

#include "staticcheck.hh"
#include "checkindexset.hh"
#include "checkgeometry.hh"
#include "checkentityseed.hh"
#include "checkentitylifetime.hh"
#include <dune/grid/test/checkidset.hh>
#include "checkintersectionlifetime.hh"

#include <limits>

// machine epsilon is multiplied by this factor
static double factorEpsilon = 1.e8;

class CheckError : public Dune::Exception {};



namespace Dune {

  // This flag can be used to disable the BoundarySegmentIndexCheck
  // Disabling this flag is important for meta grids like SubGrid and
  // MultiDomainGrid that cannot guarantee the usual boundary segment
  // index semantics
  //
  // By default, this test is on.
  template< class Grid >
  struct EnableBoundarySegmentIndexCheck
    : public std::true_type
  {};

  template< class Grid >
  struct EnableBoundarySegmentIndexCheck< const Grid >
    : public EnableBoundarySegmentIndexCheck<Grid>
  {};

}


// check
// Entity::geometry()[c] == Entity::entity<dim>.geometry()[0]
// for codim=cd
template <int cd, class Grid, class Entity, bool doCheck>
struct subIndexCheck
{

  template<typename E>
  void checkEntitySeedRecovery([[maybe_unused]] const Grid& g, const E& e)
  {
    typedef typename Grid::template Codim< E::codimension >::EntitySeed EntitySeed;
    EntitySeed seed = e.seed();

    E e1( e );
    E e2 = g.entity( seed );

    if( e1 != e2 )
      DUNE_THROW( Dune::GridError, "EntitySeed recovery failed");
  }


  subIndexCheck ( const Grid &g, const Entity &e )
  {

    checkEntitySeedRecovery(g,e);

    // check subEntity range size
    if( subEntities(e, Dune::Codim<cd>{}).size() != e.subEntities(cd) )
    {
      std::cerr << "Error: Number of subentities in range "
                << "does not match the number of subentities of the element." << std::endl;
      assert( false );
    }

    // We use the both a range-based for loop and an index counter for testing purposes here.
    // In regular code you will often only need the subIndex in which case an index loop
    // is more convenient, or you only need the actual sub-entities in which
    // case the range-based for loop might be more readable.
    std::size_t subIndex = 0;
    for (const auto& se : subEntities(e, Dune::Codim<cd>{}))
    {
      // check construction of entities
      [[maybe_unused]] auto seCopy = se;
      assert( seCopy == se );

      checkEntitySeedRecovery(g,se);

      const typename Grid::LevelGridView &levelGridView = g.levelGridView(e.level());

      if( !levelGridView.contains( se ) )
      {
        std::cerr << "Error: Level grid view does not contain all subentities." << std::endl;
        assert( false );
      }

      const typename Grid::LevelIndexSet &levelIndexSet = g.levelIndexSet( e.level() );

      if( !levelIndexSet.contains( se ) )
      {
        std::cerr << "Error: Level index set does not contain all subentities." << std::endl;
        assert( false );
      }
      if( levelIndexSet.index( se ) != levelIndexSet.subIndex( e, subIndex, cd ) )
      {
        int id_e = levelIndexSet.index( e );
        int id_e_i = levelIndexSet.index( se );
        int subid_e_i = levelIndexSet.subIndex( e, subIndex, cd );
        std::cerr << "Error: levelIndexSet.index( *(e.template subEntity< cd >( subIndex ) ) ) "
                  << "!= levelIndexSet.subIndex( e, subIndex, cd )  "
                  << "[with cd=" << cd << ", subIndex=" << subIndex << "]" << std::endl;
        std::cerr << "       ... index( e ) = " << id_e << std::endl;
        std::cerr << "       ... index( e.subEntity< cd >( subIndex ) ) = " << id_e_i << std::endl;
        std::cerr << "       ... subIndex( e, subIndex, cd ) = " << subid_e_i << std::endl;
        assert( false );
      }

      ++subIndex;
    }

    subIndexCheck< cd-1, Grid, Entity, Dune::Capabilities::hasEntity< Grid, cd-1 >::v > sick( g, e );
  }
};


// end recursion of subIndexCheck
template <class Grid, class Entity, bool doCheck>
struct subIndexCheck<-1, Grid, Entity, doCheck>
{
  subIndexCheck (const Grid & /* g */, const Entity & /* e */)
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
  subIndexCheck (const Grid & , const Entity & )
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

  const typename Grid::LevelIndexSet &levelIndexSet = g.levelIndexSet( g.maxLevel() );
  const typename Grid::LeafIndexSet &leafIndexSet = g.leafIndexSet();
  const typename Grid::LocalIdSet &localIdSet = g.localIdSet();
  const typename Grid::GlobalIdSet &globalIdSet = g.globalIdSet();

  LevelIterator it = g.levelGridView(g.maxLevel()).template begin<0>();
  const LevelIterator endit = g.levelGridView(g.maxLevel()).template end<0>();

  for (; it!=endit; ++it)
  {
    // check entity copy ctor
    Entity e( *it ) ;
    assert( e == *it );

    Entity subEntity = it->template subEntity< 0 >( 0 );
    // Entity::subEntity<0>(0) == Entity
    if( levelIndexSet.index( subEntity ) != levelIndexSet.index( *it ) )
    {
      std::cerr << "Error: Level index for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }
    if( leafIndexSet.index( subEntity ) != leafIndexSet.index( *it ) )
    {
      std::cerr << "Error: Leaf index for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }
    if( globalIdSet.id( subEntity ) != globalIdSet.id( *it ) )
    {
      std::cerr << "Error: Global id for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }
    if( localIdSet.id( subEntity ) != localIdSet.id( *it ) )
    {
      std::cerr << "Error: Local id for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }

    if( subEntity.level() != it->level() )
    {
      std::cerr << "Error: Level for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }

    // Entity::subEntities(dim) == Entity::geometry().corners();
    // Entity::geometry()[c] == Entity::entity<dim>.geometry()[0];
    const int numCorners = it->subEntities(dimGrid);
    if( numCorners != it->geometry().corners() )
    {
      std::cerr << "Error: Entity::subEntities< dimGrid >() != Entity::geometry().corners()." << std::endl;
      assert( false );
    }
    for( int c = 0; c < numCorners; ++c )
    {
      typename Geometry::GlobalCoordinate c1( it->geometry().corner( c ) );
      typename Geometry::GlobalCoordinate c2( it->template subEntity< dimGrid >( c ).geometry().corner( 0 ) );

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
  [[maybe_unused]] constexpr static int dim = Grid::dimension;
  // [[maybe_unused]] typedef typename Grid::ctype ct;

  typedef typename Grid::GlobalIdSet GlobalIdSet;
  const GlobalIdSet & globalid = g.globalIdSet();

  typedef typename Grid::template Codim< 0 >::Entity Entity;

  GridView gridView = g.levelGridView( 0 );

  LevelIterator eit = g.levelGridView(0).template begin<0>();
  const LevelIterator eend = g.levelGridView(0).template end<0>();

  LevelIterator next = eit;
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

    // small check if LevelIntersectionIterators
    // from work reassigned Entity
    // after creation of LevelIterator on different level
    if (g.maxLevel()>0
        && g.levelGridView(0).template begin<0>() != g.levelGridView(0).template end<0>()
        && g.levelGridView(1).template begin<0>() != g.levelGridView(1).template end<0>()
        )
    {
      Entity e( *g.levelGridView(0).template begin<0>() );
      e = *g.levelGridView(1).template begin<0>();
      LevelIterator it = g.levelGridView(0).template begin<0>();
      ++it;
      g.levelGridView(e.level()).ibegin(e);
    }

    // iterate over level
    for (; eit != eend; ++eit)
    {
      const Entity &entity = *eit;

      const bool isGhost = (entity.partitionType() == Dune::GhostEntity);

      // call global id
      globalid.id( entity );

      const int numFaces = entity.subEntities(1);
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
        typedef typename IntersectionIterator::Intersection Intersection;
        const Intersection& is = *it;

        // check intersection copy ctor
        [[maybe_unused]] Intersection is2 = is;

        // numbering
        const int numberInSelf = is.indexInInside();
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
        assert( globalid.id(it->inside()) == globalid.id( entity ) );

        // geometry
        if( !it->type().isNone() )
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
          const Entity outside = it->outside();

          if( isGhost && (outside.partitionType() == Dune::GhostEntity) )
          {
            std::cerr << "Error: Intersections between ghosts shall not be returned." << std::endl;
            assert( false );
          }

          // geometry
          if( !it->type().isNone() )
            it->geometryInOutside();

          // numbering
          const int indexInOutside = it->indexInOutside();
          const int outNumFaces = outside.subEntities(1);
          if( (indexInOutside < 0) || (indexInOutside >= outNumFaces) )
          {
            std :: cout << "Error: Invalid indexInOutside: " << indexInOutside
                        << " (should be between 0 and " << (outNumFaces-1) << ")"
                        << std :: endl;
            assert( false );
          }

          // search neighbouring cell
          if( outside.partitionType() == Dune::InteriorEntity )
          {
            assert( globalid.id( outside ) != globalid.id( entity ) );
            const Dune::PartitionIteratorType pitype = Dune::InteriorBorder_Partition;

            typedef typename Grid::template Codim< 0 >
            ::template Partition< pitype >::LevelIterator
            LevelIteratorI;

            const int level = entity.level();
            bool foundNeighbor = false;
            LevelIteratorI nit = g.levelGridView(level).template begin< 0, pitype >();
            const LevelIteratorI nend = g.levelGridView(level).template end< 0, pitype > ();
            for( ; nit != nend; ++nit )
            {
              if( nit->partitionType() != Dune::InteriorEntity )
              {
                std::cerr << "Error: LevelIterator for InteriorBorder_Partition "
                          << "stops on non-interior entity." << std :: endl;
                assert( false );
              }

              if( *nit != outside )
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
  static void check(const GridType &, IteratorType & )
  {}
};

/*
 * Iterate over the grid und do some runtime checks
 */

template <bool checkMark , class Grid>
void iterate(Grid &g)
{
  typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
  //[[maybe_unused]] typedef typename Grid::HierarchicIterator HierarchicIterator;
  typedef typename Grid::template Codim<0>::Geometry Geometry;
  int maxLevel = g.maxLevel();
  LevelIterator it = g.levelGridView(maxLevel).template begin<0>();
  const LevelIterator endit = g.levelGridView(maxLevel).template end<0>();

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

    const Geometry geo = it->geometry();

    if( geo.type() != it->type() )
    {
      std::cerr << "inconsistent geometry type: entity.type() = " << it->type()
                << ", geometry.type() = " << geo.type() << "." << std::endl;
      assert( false );
    }

    if( !geo.type().isNone() )
    {
      origin = referenceElement(geo).position(0,0);
      typename Geometry::LocalCoordinate origin2 =
        referenceElement( *it ).position(0,0);
      if( (origin - origin2 ).two_norm() > 1e-10 )
        DUNE_THROW(CheckError, "referenceElement( entity ) != referenceElement( geo )");

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
    }

    geo.corners();
    geo.corner( 0 );
  }

  typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
  LeafIterator lit = g.leafGridView().template begin<0>();
  const LeafIterator lend = g.leafGridView().template end<0>();

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

    if( !lit->type().isNone() )
    {
      origin = referenceElement(lit->geometry()).position(0,0);
      result = lit->geometry().local(lit->geometry().global(origin));
      typename Grid::ctype error = (result-origin).two_norm();
      if(error >= factorEpsilon * std::numeric_limits<typename Grid::ctype>::epsilon())
      {
        DUNE_THROW(CheckError, "|| geom.local(geom.global(" << origin
                                                            << ")) - origin || != 0 ( || " << result << " - origin || ) = " << error);
      }
      lit->geometry().integrationElement(origin);
      if((int)Geometry::coorddimension == (int)Geometry::mydimension)
        lit->geometry().jacobianInverseTransposed(origin);
    }

    lit->geometry().type();
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

  typedef typename Grid::LeafGridView LeafGridView;
  typedef typename Grid::LevelGridView LevelGridView;
  typedef typename LeafGridView::IntersectionIterator LeafIntersectionIterator;
  typedef typename LevelGridView::IntersectionIterator LevelIntersectionIterator;

  LeafGridView leafGridView = g.leafGridView();
  LevelGridView levelGridView = g.levelGridView(0);

  // assignment tests
  LevelIterator l1 = levelGridView.template begin<0>();
  LevelIterator l2 = levelGridView.template begin<0>();
  LeafIterator L1 = leafGridView.template begin<0>();
  LeafIterator L2 = leafGridView.template begin<0>();

  // if grid empty, leave
  if (l1 == levelGridView.template end<0>())
    return;

  // check '==' consistency
  typedef typename Grid::template Codim<0>::Entity Entity;
  Entity a( *levelGridView.template begin<0,Dune::All_Partition>() );
  Entity i( *levelGridView.template begin<0,Dune::InteriorBorder_Partition>() );

  assert(
    (levelGridView.indexSet().index(a) != levelGridView.indexSet().index(i)) // index equal
    == (a != i) // entitypointer equal
    );
  assert(
    (levelGridView.indexSet().index(a) == levelGridView.indexSet().index(i)) // index equal
    == (a == i) // entitypointer equal
    );

  HierarchicIterator h1 = l1->hbegin(99);
  HierarchicIterator h2 = l2->hbegin(99);
  LeafIntersectionIterator i1 = leafGridView.ibegin( *l1 );
  LeafIntersectionIterator i2 = leafGridView.ibegin( *l2 );

  // assign
  l1 = l2;
  L1 = L2;
  h1 = h2;
  i1 = i2;

  // iterator comparisons
  assert(l1 == l2 && !(l1 != l2));
  assert(L1 == L2 && !(L1 != L2));
  assert(h1 == h2 && !(h1 != h2));
  assert(i1 == i2 && !(i1 != i2));

  // entity comparisons
  assert(*l1 == *l2 && !(*l1 != *l2));
  assert(*L1 == *L2 && !(*L1 != *L2));
  if (h1 != l1->hend(99))
    assert(*h1 == *h2 && !(*h1 != *h2));

  // intersection comparison
  if (i1 != leafGridView.iend(*l1))
    assert(*i1 == *i2 && !(*i1 != *i2));

  if( EnableLevelIntersectionIteratorCheck< Grid >::v )
  {
    LevelIntersectionIterator li1 = levelGridView.ibegin( *l1 );
    LevelIntersectionIterator li2 = levelGridView.ibegin( *l2 );
    assert(li1 == li2 && !(li1 != li2));
    if( !( li1 == li2 && !(li1 != li2) ) )
      DUNE_THROW( Dune::GridError, "LevelIntersectionIterator check failed" );

    if (li1 != levelGridView.iend(*l1))
    {
      assert(*li1 == *li2 && !(*li1 != *li2));
      if( ! ( *li1 == *li2 && !(*li1 != *li2) ) )
        DUNE_THROW( Dune::GridError, "LevelIntersectionIterator dereferencing check failed" );
    }
  }

}


template< class Grid, class Entity, class Intersection >
void checkBoundarySegmentIndexProlongation ( const Grid &grid, const Entity &entity, const Intersection &intersection )
{
  typedef typename Entity::HierarchicIterator HierarchicIterator;
  typedef typename Entity::LocalGeometry GeometryInFather;
  typedef typename Grid::LevelGridView GridView;
  typedef typename GridView::IntersectionIterator IntersectionIterator;

  const typename Intersection::LocalGeometry &geoInInside = intersection.geometryInInside();

  const GridView gridView = grid.levelGridView( entity.level()+1 );
  const HierarchicIterator hend = entity.hend( entity.level()+1 );
  for( HierarchicIterator hit = entity.hbegin( entity.level()+1 ); hit != hend; ++hit )
  {
    GeometryInFather geoInFather = hit->geometryInFather();
    auto refElement = referenceElement( geoInFather );

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
    typedef typename Iterator::Entity Entity;

    GridView gv = grid.levelGridView(level);
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
            Entity f = it->father();

            if (f.level() != level-1)
            {
              std::cerr << "Error! father().level()=" << f.level()
                        << " for element on level " << level << std::endl;
              assert(false);
            }

          }
          // check level of father after reassignment
          {
            Entity f(*it);
            f = f.father();

            if (f.level() != level-1)
            {
              std::cerr << "Error! father().level()=" << f.level()
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
  bool error = false;

  const Iterator end = gridView.template end< 0 >();
  for( Iterator it = gridView.template begin< 0 >(); it != end; ++it )
  {
    if( !it->hasBoundaryIntersections() )
      continue;
    const IntersectionIterator iend = gridView.iend( *it );
    for( IntersectionIterator iit = gridView.ibegin( *it ); iit != iend; ++iit )
    {
      if( iit->boundary() )
        ++countBoundarySegments;
    }
  }

  if( countBoundarySegments != numBoundarySegments )
  {
    std::cerr << "Error: Wrong number of boundary segments (reported: "
              << numBoundarySegments << ", counted: "
              << countBoundarySegments << ")." << std::endl;
    std::abort();
  }

  std::vector< int > count( numBoundarySegments, 0 );
  for( Iterator it = gridView.template begin< 0 >(); it != end; ++it )
  {
    if( !it->hasBoundaryIntersections() )
      continue;
    const IntersectionIterator iend = gridView.iend( *it );
    for( IntersectionIterator iit = gridView.ibegin( *it ); iit != iend; ++iit )
    {
      if( !iit->boundary() )
        continue;
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


// Test whether the local-to-global map for codim 1 subentities is a bijective
// mapping from the reference element onto the global element.
// This is hard to test in general, but for elements that lie in an axis-parallel plane
// we can just ignore one coordinate and check whether the determinant of the
// Jacobian with respect to the other two coordinates is non-zero everywhere. As we
// can't actually check "everywhere", test whether the determinant changes sign between
// different quadrature points. As the determinant should be continuous, a sign flip
// implies that the determinant is zero at some point.
template <class Grid>
typename std::enable_if<Grid::dimension == 3, void>::type checkCodim1Mapping(const Grid &g)
{
  bool error = false;
  const auto& gv = g.leafGridView();
  for (const auto& entity : Dune::elements(gv)) {
    for (unsigned int index = 0; index < entity.subEntities(1); ++index) {
      const auto& subEntity = entity.template subEntity<1>(index);
      const auto& subGeom = subEntity.geometry();
      const auto numCorners = subGeom.corners();
      std::bitset<3> equalCoords("111");
      const auto& firstCornerCoords = subGeom.corner(0);
      for (int c = 1; c < numCorners; ++c)
        for (size_t i = 0; i < 3; ++i)
          equalCoords[i] = equalCoords[i] && Dune::FloatCmp::eq(subGeom.corner(c)[i],
                                                                firstCornerCoords[i]);
      if (equalCoords.count() == 1) {
        size_t detPositive(0), detNegative(0), detZero(0), size(0);
        for (const auto& ip : Dune::QuadratureRules<double,2>::rule(subEntity.type(),2)) {
          ++size;
          const Dune::FieldMatrix< typename Grid::ctype, 2, 3 > full_jacT( subGeom.jacobianTransposed( ip.position() ) );
          Dune::FieldMatrix<double,2,2> jacT(0);
          size_t k = 0;
          for (size_t j = 0; j < 3; ++j) {
            if (!equalCoords[j]) {
              for (size_t i = 0; i < 2; ++i)
                jacT[i][k] += full_jacT[i][j];
              ++k;
            }
          }
          const auto det = jacT.determinant();
          det > 0 ? ++detPositive : (det < 0 ? ++detNegative : ++detZero);
        }
        if (!(detPositive == size || detNegative == size))
        {
          std::cerr << "Error: Local-to-global mapping for codim 1 subEntity "
                    << index << " of element with index "
                    << gv.indexSet().index(entity) << " is not invertible!"
                    << std::endl;
          error = true;
        }
      }
    }
  }
  if (error)
  {
    std::string msg("Error encountered during checkCodim1Mapping!.");
    if( g.comm().size() > 1 )
      std::cerr << msg << std::endl;
    else
      DUNE_THROW(Dune::Exception, msg );
  }
}

template <class GridView>
void checkConforming(const GridView& view)
{
  view.isConforming();
}

template<class Grid>
typename std::enable_if<Grid::dimension != 3, void>::type checkCodim1Mapping(const Grid &/*g*/)
{}

template <class Grid>
void gridcheck (Grid &g)
{
  /*
   * first do the compile-test: this will not produce any code but
   * fails if an interface-component is missing
   */
  GridInterface<Grid>();

  constexpr static int dim = Grid :: dimension;
  constexpr static int dimworld = Grid :: dimensionworld;
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
      // if entities for a specific codimension are not available the size can be 0
      if( g.size(i) > 0 )
      {
        for (int j=0; j<=g.maxLevel(); j++)
          if (g.comm().size() == 1)
            assert(g.size(j,i)>0);
          else
            assert(g.size(j,i)>=0);   // in parallel this could be zero (depends on definition of maxLevel() method in parallel)
      }
    }
  }

#if DUNE_VERBOSE_TESTS
  // print infos
  Dune::gridinfo(g, "GridInfo");
  for(int l=0; l<g.maxLevel(); l++)
    Dune::gridlevellist(g, l, "GridLevelInfo");
  Dune::gridleaflist(g, "GridLeafInfo");
#endif

  // check functionality when grid is interpreted as reference to interface
  GridInterface< Dune::Grid< dim, dimworld, ctype, GridFamily > >();

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
  checkCodim1Mapping(g);

  // check conforming
  checkConforming( g.levelGridView( 0 ) );
  checkConforming( g.leafGridView() );

  // check geometries of macro level and leaf level
  Dune::GeometryChecker<Grid> checker;
  checker.checkGeometry( g.levelGridView( 0 ) );
  checker.checkGeometry( g.leafGridView() );

  // check entity seeds
  Dune::checkEntitySeed( g.leafGridView(), std::cerr );
  for( int level = 0; level <= g.maxLevel(); ++level )
    Dune::checkEntitySeed( g.levelGridView( level ), std::cerr );

  // note that for some grid this might fail
  // then un comment this test
  Dune :: checkIndexSet( g, g.leafGridView(), Dune :: dvverb );
  for( int level = 0; level <= g.maxLevel(); ++level )
    Dune :: checkIndexSet( g, g.levelGridView( level ), Dune :: dvverb, true );

  // check id sets
  checkIdSet(g, g.localIdSet());
  checkIdSet(g, g.globalIdSet());

  // check at least if the subId method is there
  {
    typename Grid::Traits::template Codim<0>::LevelIterator it =
      g.levelGridView( g.maxLevel() ).template begin<0>();
    typename Grid::Traits::template Codim<0>::LevelIterator end =
      g.levelGridView( g.maxLevel() ).template end<0>();
    for (; it != end; ++it)
    {
      g.globalIdSet().subId(*it,0,dim);
      if(g.globalIdSet().subId(*it,0,dim) != g.globalIdSet().id(it->template subEntity<dim>(0)))
      {
        std::cerr << "Error: Inconsistent global subId for vertex 0 (id(subEntity)=" << g.globalIdSet().id(it->template subEntity<dim>(0))
                  << ", subId=" << g.globalIdSet().subId(*it,0,dim) << ") of cell " << g.globalIdSet().id(*it) << std::endl;
        assert(false);
      }
    }
  }

  if(Dune::EnableBoundarySegmentIndexCheck<Grid>::value)
    {
      if( EnableLevelIntersectionIteratorCheck< Grid >::v )
        checkBoundarySegmentIndex( g.levelGridView( 0 ) );
      else if( g.maxLevel() == 0 )
        checkBoundarySegmentIndex( g.leafGridView() );
      else
        std::cout << "Warning: Skipping boundary segment index check (missing level intersection iterator)." << std::endl;
    }
  else
    std::cout << "Warning: Skipping boundary segment index check because it has been explicitly disabled." << std::endl;

  // check for range-based for iteration and correct entity / intersection lifetime
  checkEntityLifetime(g.levelGridView(g.maxLevel()));
  if (EnableLevelIntersectionIteratorCheck< Grid >::v)
    checkIntersectionLifetime(g.levelGridView(g.maxLevel()));
  checkEntityLifetime(g.leafGridView());
  checkIntersectionLifetime(g.leafGridView());
}

#endif // #ifndef DUNE_GRID_TEST_GRIDCHECK_HH
