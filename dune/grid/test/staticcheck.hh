// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_STATICCHECK_HH
#define DUNE_GRID_TEST_STATICCHECK_HH

/** \file
    \brief Implements static grid checks


   \todo check return types

 */

#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/deprecated.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridenums.hh>

template< class Geometry, int codim, int dim >
struct GeometryInterface
{
  static void check ( const Geometry &geo )
  {
    static_assert( (Geometry::mydimension == dim-codim), "" );

    typedef typename Geometry::ctype ctype DUNE_UNUSED;

    geo.type();
    geo.affine();
    geo.corners();
    geo.corner( 0 );

    typename Geometry::LocalCoordinate v(0.0);
    geo.global(v);
    typename Geometry::GlobalCoordinate g(0.0);
    geo.local(g);
    geo.integrationElement(v);
    geo.jacobianTransposed( v );
    geo.jacobianInverseTransposed( v );
  }

  GeometryInterface ()
  {
    c = check;
  }

private:
  void (*c)( const Geometry & );
};

// --- compile-time check of entity-interface

// tests that should work on entities of all codimensions
template <class Entity>
void DoEntityInterfaceCheck (Entity &e)
{
  // exported types
  typedef typename Entity::Geometry Geometry DUNE_UNUSED;
  typedef typename Entity::EntitySeed EntitySeed DUNE_UNUSED;

  // methods on each entity
  e.level();
  e.partitionType();
  e.geometry();

  // check interface of attached element-interface
  GeometryInterface< typename Entity::Geometry, Entity::codimension, Entity::dimension >();
}

// recursive check of codim-0-entity methods subEntities(), entity()
template <class Grid, int cd, bool hasEntity>
struct ZeroEntityMethodCheck
{
  typedef typename Grid::template Codim<0>::Entity Entity;
  static void check(Entity &e)
  {
    // check types
    typedef typename Entity::HierarchicIterator HierarchicIterator DUNE_UNUSED;

    e.subEntities(cd);
    e.template subEntity<cd>(0);

    // recursively check on
    ZeroEntityMethodCheck<Grid, cd - 1,
        Dune::Capabilities::hasEntity<Grid, cd - 1>::v >();
  }
  ZeroEntityMethodCheck ()
  {
    c = check;
  }
  void (*c)(Entity &e);
};

// just the recursion if the grid does not know about this codim-entity
template<class Grid, int cd>
struct ZeroEntityMethodCheck<Grid, cd, false>
{
  typedef typename Grid::template Codim<0>::Entity Entity;
  static void check(Entity &)
  {
    // check types
    typedef typename Entity::HierarchicIterator HierarchicIterator DUNE_UNUSED;

    // recursively check on
    ZeroEntityMethodCheck<Grid, cd - 1,
        Dune::Capabilities::hasEntity<Grid, cd - 1>::v >();
  }
  ZeroEntityMethodCheck ()
  {
    c = check;
  }
  void (*c)(Entity &e);
};

// end recursive checking
template <class Grid>
struct ZeroEntityMethodCheck<Grid, 0, true>
{
  typedef typename Grid::template Codim<0>::Entity Entity;
  static void check(Entity &e)
  {
    // check types
    typedef typename Entity::HierarchicIterator HierarchicIterator DUNE_UNUSED;

    e.subEntities(0);
    e.template subEntity<0>(0);

  }
  ZeroEntityMethodCheck ()
  {
    c = check;
  }
  void (*c)(Entity &e);
};

// end recursive checking - same as true
// ... codim 0 is always needed
template <class Grid>
struct ZeroEntityMethodCheck<Grid, 0, false>
{
  typedef typename Grid::template Codim<0>::Entity Entity;
  static void check(Entity &e)
  {
    // check types
    typedef typename Entity::HierarchicIterator HierarchicIterator DUNE_UNUSED;

    e.subEntities(0);
    e.template subEntity<0>(0);
  }
  ZeroEntityMethodCheck ()
  {
    c = check;
  }
  void (*c)(Entity &e);
};

// IntersectionIterator interface check
template <class Grid,class IntersectionIterator>
struct IntersectionIteratorInterface
{
  enum { dim = Grid::dimension };
  typedef typename Grid::ctype ct;

  static void check (IntersectionIterator &i)
  {
    // increment / equality / ...
    IntersectionIterator j = i;
    ++j;
    { bool DUNE_UNUSED tmp = (i == j); }
    { bool DUNE_UNUSED tmp = (i != j); }
    j = i;

    // state
    typedef typename IntersectionIterator::Intersection Intersection;
    const Intersection &inter = *i;
    inter.boundary();
    inter.neighbor();

    inter.boundarySegmentIndex();

    // neighbouring elements
    inter.inside();
    if(inter.neighbor()) inter.outside();

    // geometry
    inter.geometryInInside();
    if(inter.neighbor()) inter.geometryInOutside();
    inter.geometry();

    inter.indexInInside();
    if(inter.neighbor()) inter.indexInOutside();

    typename Intersection::LocalCoordinate v(0);

    inter.outerNormal(v);
    inter.integrationOuterNormal(v);
    inter.unitOuterNormal(v);
  }
  IntersectionIteratorInterface ()
  {
    c = check;
  }
private:
  void (*c)(IntersectionIterator&);
};

// check codim-entity and pass on to codim + 1
template <class Grid, int codim, int dim, bool hasEntity>
struct EntityInterface
{
  typedef typename Grid::template Codim<codim>::Entity Entity;

  static void check ( const Entity &e )
  {
    // consistent?
    static_assert( (Entity::codimension == codim), "" );
    static_assert( (Entity::dimension == dim), "" );

    // do the checking
    DoEntityInterfaceCheck(e);

    // recursively check sub-entities
    EntityInterface<Grid, codim + 1, dim,
        Dune::Capabilities::hasEntity<Grid, codim + 1>::v >();
  }
  EntityInterface ()
  {
    c = check;
  }

private:
  void (*c)( const Entity & );
};

// just the recursion if the grid does not know about this codim-entity
template <class Grid, int codim, int dim>
struct EntityInterface<Grid, codim, dim, false>
{
  typedef typename Grid::template Codim<codim>::Entity Entity;

  static void check (Entity &)
  {
    // recursively check sub-entities
    EntityInterface<Grid, codim + 1, dim,
        Dune::Capabilities::hasEntity<Grid, codim + 1>::v >();
  }
  EntityInterface ()
  {
    c = check;
  }
  void (*c)(Entity&);
};

// codim-0 entities have different interface
template <class Grid, int dim>
struct EntityInterface<Grid, 0, dim, true>
{
  typedef typename Grid::template Codim<0>::Entity Entity;

  static void check ( const Entity &e, bool checkLevelIter = true )
  {
    // consistent?
    static_assert( (Entity::codimension == 0), "" );
    static_assert( (Entity::dimension == dim), "" );

    // do the common checking
    DoEntityInterfaceCheck(e);

    // special codim-0-entity methods which are parametrized by a codimension
    ZeroEntityMethodCheck
    <Grid, dim, Dune::Capabilities::hasEntity<Grid, dim>::v >();

    // grid hierarchy
    if ( e.hasFather() )
    {
      const Entity father = e.father();
      father.hbegin(0);
      e.geometryInFather();
    }

    // intersection iterator
    if (checkLevelIter) {
      // #if 0 // WARNING must be updated to new interface
      IntersectionIteratorInterface<Grid,typename Grid::LevelIntersectionIterator>();
      // #endif
    }
    // #if 0 // WARNING must be updated to new interface
    if(e.isLeaf())
      IntersectionIteratorInterface<Grid,typename Grid::LevelIntersectionIterator>();
    // #endif

    // hierarchic iterator
    e.hbegin(0);
    e.hend(0);

    // adaption
    e.isNew();
    e.mightVanish();

    // recursively check sub-entities
    EntityInterface<Grid, 1, dim,
        Dune::Capabilities::hasEntity<Grid, 1>::v >();
  }
  EntityInterface ()
  {
    c = check;
  }

private:
  void (*c)( const Entity &, bool );
};

// non existinng codim-0 entity
template <class Grid, int dim>
struct EntityInterface<Grid, 0, dim, false>
{
  typedef typename Grid::template Codim<0>::Entity Entity;

  static void check (Entity &e)
  {
    // recursively check sub-entities
    EntityInterface<Grid, 1, dim,
        Dune::Capabilities::hasEntity<Grid, 1>::v >();
  }
  EntityInterface ()
  {
    c = check;
  }
  void (*c)(Entity&);
};

// end the recursion over entity-codimensions
template <class Grid, int dim>
struct EntityInterface<Grid, dim, dim, true>
{
  typedef typename Grid::template Codim<dim>::Entity Entity;

  // end recursion
  static void check ( const Entity &e )
  {
    // consistent?
    static_assert( (Entity::codimension == dim), "" );
    static_assert( (Entity::dimension == dim), "" );

    // run common test
    DoEntityInterfaceCheck(e);
  }

  EntityInterface()
  {
    c = check;
  }
  void (*c)( const Entity & );
};

// end the recursion over entity-codimensions
// ... codim dim entity does not exist
template <class Grid, int dim>
struct EntityInterface<Grid, dim, dim, false>
{
  typedef typename Grid::template Codim<dim>::Entity Entity;

  // end recursion
  static void check (Entity &e)
  {}

  EntityInterface()
  {
    c = check;
  }
  void (*c)(Entity&);
};

template<class Grid>
struct LeafInterface
{
  static void check(Grid &g)
  {
    DUNE_UNUSED_PARAMETER(g);
  }

  LeafInterface()
  {
    c = check;
  }
  void (*c)(Grid&);
};

template< class GridView >
struct GridViewInterface
{

  static void check ( const GridView &gv )
  {
    const int dimension = GridView::dimension;

    typedef typename GridView::Grid Grid DUNE_UNUSED;
    typedef typename GridView::IndexSet IndexSet DUNE_UNUSED;

    typedef typename GridView::Intersection Intersection DUNE_UNUSED;
    typedef typename GridView::IntersectionIterator IntersectionIterator DUNE_UNUSED;

    gv.grid();

    gv.size( 0 );
    gv.size( Dune::GeometryTypes::cube( GridView::dimension ) );

    gv.template begin< 0 >();
    gv.template end< 0 >();

    // index set
    gv.indexSet();
    using namespace Dune::Hybrid;
    forEach(std::make_integer_sequence< int, dimension+1 >(), [&](auto codim) {
      typedef typename GridView::template Codim< codim >::Entity Entity;
      typedef typename GridView::template Codim< codim >::Iterator Iterator DUNE_UNUSED;

      if( gv.template begin< 0 >() == gv.template end< 0 >() )
        return;

      const Entity &entity = gv.template begin< 0 >()->template subEntity< codim >( 0 );
      gv.indexSet().index( entity );
      gv.indexSet().contains( entity );
      try
      {
        forEach(std::make_integer_sequence< int, GridView::dimension+1 - codim>(), [&](auto subCodim) {
          gv.indexSet().subIndex( entity, 0, codim+subCodim);
        });
      }
      catch( Dune::NotImplemented )
      {
        // ignore Dune::NotImplemented for higher codimension
        if( codim == 0 )
          throw;
      }

      gv.indexSet().types( codim );
    });

    // intersections
    if( gv.template begin< 0 >() != gv.template end< 0 >() )
    {
      gv.ibegin( *gv.template begin< 0 >() );
      gv.iend( *gv.template begin< 0 >() );
    }

    IntersectionIteratorInterface< Grid, IntersectionIterator >();

    // parallel interface
    typedef typename GridView::template Codim< 0 >::template Partition< Dune::Ghost_Partition >::Iterator GhostIterator DUNE_UNUSED;
    typedef typename GridView::CollectiveCommunication CollectiveCommunication DUNE_UNUSED;

    gv.template begin< 0, Dune::Ghost_Partition >();
    gv.template end< 0, Dune::Ghost_Partition >();

    gv.overlapSize( 0 );
    gv.ghostSize( 0 );

    gv.comm();
  }

  GridViewInterface () : c( check ) {}

  void (*c)( const GridView & );
};

/** \brief Instantiate this class for a full static interface check */
template <class Grid>
struct GridInterface
{
  static void check (const Grid &g)
  {
    // check for exported types
    typedef typename Grid::LevelGridView LevelGridView DUNE_UNUSED;
    typedef typename Grid::LeafGridView LeafGridView DUNE_UNUSED;

    typedef typename Grid::ctype ctype DUNE_UNUSED;

#if !DISABLE_DEPRECATED_METHOD_CHECK
    typedef typename Grid::template Codim<0>::LevelIterator LevelIterator DUNE_UNUSED;
    typedef typename Grid::template Codim<0>::LeafIterator LeafIterator DUNE_UNUSED;
#endif // #if !DISABLE_DEPRECATED_METHOD_CHECK

    // check for grid views
    g.levelGridView( 0 );
    g.leafGridView();
    levelGridView(g, 0);
    leafGridView(g);

    // check for member functions
    g.maxLevel();
    // number of grid entities of a given codim on a given level
    g.size(0,0);
    // number of leaf entities per codim in this process
    g.size(0);
    // number of entities per level and geometry type in this process
    g.size(0, Dune::GeometryTypes::cube(Grid::dimension));
    // number of leaf entities per geometry type in this process
    g.size(Dune::GeometryTypes::cube(Grid::dimension));

    // check for iterator functions
    g.levelGridView(0).template begin<0>();
    g.levelGridView(0).template end<0>();

    LeafInterface< Grid >();

    GridViewInterface< LevelGridView >();
    GridViewInterface< LeafGridView >();

    // Check for index sets
    typedef typename Grid::LevelIndexSet LevelIndexSet DUNE_UNUSED;
    typedef typename Grid::LeafIndexSet LeafIndexSet DUNE_UNUSED;
    typedef typename Grid::LocalIdSet LocalIdSet DUNE_UNUSED;
    typedef typename Grid::GlobalIdSet GlobalIdSet DUNE_UNUSED;

    g.levelIndexSet(0);

    if (g.levelGridView(0).template begin<0>() != g.levelGridView(0).template end<0>())
    {
      // Instantiate all methods of LevelIndexSet
      g.levelIndexSet(0).index(*g.levelGridView(0).template begin<0>());
      /** \todo Test for subindex is missing, because I don't know yet
          how to test for the existence of certain codims */
    }

    g.levelIndexSet(0).
    size(Dune::GeometryTypes::simplex(Grid::dimension));
    for( int codim = 0; codim < Grid::dimension; ++codim )
      g.levelIndexSet( 0 ).types( codim );

    /** \todo Test for subindex is missing, because I don't know yet
       how to test for the existence of certain codims */
    g.leafIndexSet().size(Dune::GeometryTypes::simplex(Grid::dimension));
    for( int codim = 0; codim < Grid::dimension; ++codim )
      g.leafIndexSet().types( codim );

    if ( g.levelGridView( 0 ).template begin< 0 >() != g.levelGridView( 0 ).template end< 0 >() )
    {
      // Instantiate all methods of LocalIdSet
      /** \todo Test for subindex is missing, because I don't know yet
          how to test for the existence of certain codims */
      g.localIdSet().id( *g.levelGridView( 0 ).template begin< 0 >() );
      // Instantiate all methods of GlobalIdSet
      /** \todo Test for subindex is missing, because I don't know yet
          how to test for the existence of certain codims */
      g.globalIdSet().id( *g.levelGridView( 0 ).template begin< 0 >() );
    }
    // recursively check entity-interface
    // ... we only allow grids with codim 0 zero entites
    static_assert((Dune::Capabilities::hasEntity<Grid, 0>::v),"Grid must have codim 0 entities");
    static_assert((Dune::Capabilities::hasEntity<const Grid, 0>::v),"Grid must have codim 0 entities");

    EntityInterface< Grid, 0, Grid::dimension, Dune::Capabilities::hasEntity< Grid, 0 >::v >();

    // !!! check for parallel grid?
    g.levelGridView(0).template begin<0, Dune::Ghost_Partition>();
    g.levelGridView(0).template end<0, Dune::Ghost_Partition>();
  }
  GridInterface()
  {
    c = check;
  }
  // member just to avoid "unused variable"-warning in constructor
  void (*c)(const Grid&);
};

#endif // #ifndef DUNE_GRID_TEST_STATICCHECK_HH
