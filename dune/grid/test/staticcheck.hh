// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id: gridcheck.cc 6156 2010-01-17 16:40:14Z dedner $
#ifndef DUNE_STATIC_CHECK_HH
#define DUNE_STATIC_CHECK_HH

/** \file
    \brief Implements static grid checks


   \todo check return types

 */

#include <dune/geometry/type.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridenums.hh>

template< class Geometry, int codim, int dim >
struct GeometryInterface
{
  static void check ( const Geometry &geo )
  {
    dune_static_assert( (Geometry::mydimension == dim-codim), "" );
    dune_static_assert( (Geometry::dimension == dim), "" );

    typedef typename Geometry::ctype ctype;

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
  typedef typename Entity::ctype ctype;

  // methods on each entity
  e.level();
  e.partitionType();
  e.geometry();

  // check interface of attached element-interface
  GeometryInterface< typename Entity::Geometry, Entity::codimension, Entity::dimension >();
}

// recursive check of codim-0-entity methods count(), entity()
template <class Grid, int cd, bool hasEntity>
struct ZeroEntityMethodCheck
{
  typedef typename Grid::template Codim<0>::Entity Entity;
  static void check(Entity &e)
  {
    // check types
    typedef typename Entity::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename Entity::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename Entity::HierarchicIterator HierarchicIterator;
    typedef typename Entity::EntityPointer EntityPointer;

    e.template count<cd>();
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
  static void check(Entity &e)
  {
    // check types
    typedef typename Entity::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename Entity::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename Entity::HierarchicIterator HierarchicIterator;
    typedef typename Entity::EntityPointer EntityPointer;

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
    typedef typename Entity::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename Entity::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename Entity::HierarchicIterator HierarchicIterator;
    typedef typename Entity::EntityPointer EntityPointer;

    e.template count<0>();
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
    typedef typename Entity::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename Entity::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename Entity::HierarchicIterator HierarchicIterator;
    typedef typename Entity::EntityPointer EntityPointer;

    e.template count<0>();
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
    i == j;
    i != j;
    j = i;

    // state
    typedef typename IntersectionIterator::Intersection Intersection;
    const Intersection &inter = *i;
    inter.boundary();
    inter.neighbor();

    inter.boundarySegmentIndex();
#if !DISABLE_DEPRECATED_METHOD_CHECK
    // id of boundary segment
    inter.boundaryId();
#endif // #if !DISABLE_DEPRECATED_METHOD_CHECK

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
    dune_static_assert( (Entity::codimension == codim), "" );
    dune_static_assert( (Entity::dimension == dim), "" );

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

  static void check (Entity &e)
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
    dune_static_assert( (Entity::codimension == 0), "" );
    dune_static_assert( (Entity::dimension == dim), "" );

    // do the common checking
    DoEntityInterfaceCheck(e);

    // special codim-0-entity methods which are parametrized by a codimension
    ZeroEntityMethodCheck
    <Grid, dim, Dune::Capabilities::hasEntity<Grid, dim>::v >();

    // grid hierarchy
    if ( e.hasFather() )
    {
      const typename Entity::EntityPointer fatherPtr = e.father();
      const Entity &father = *fatherPtr;
      father.hbegin(0);
      e.geometryInFather();
    }

    // intersection iterator
    if (checkLevelIter) {
#if !DISABLE_DEPRECATED_METHOD_CHECK
      e.ilevelbegin();
      e.ilevelend();
#endif // #if !DISABLE_DEPRECATED_METHOD_CHECK
      // #if 0 // WARNING must be updated to new interface
      IntersectionIteratorInterface<Grid,typename Grid::LevelIntersectionIterator>();
      // #endif
    }
#if !DISABLE_DEPRECATED_METHOD_CHECK
    e.ileafbegin();
    e.ileafend();
#endif // #if !DISABLE_DEPRECATED_METHOD_CHECK
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
    dune_static_assert( (Entity::codimension == dim), "" );
    dune_static_assert( (Entity::dimension == dim), "" );

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
    g.template leafbegin<0>();
    g.template leafend<0>();
  }
  LeafInterface()
  {
    c = check;
  }
  void (*c)(Grid&);
};

/** \brief Instantiate this class for a full static interface check */
template <class Grid>
struct GridInterface
{
  static void check (const Grid &g)
  {
    // check for exported types
    typedef typename Grid::ctype ctype;
    typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
    typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;
    typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;

    // check for member functions
    g.maxLevel();
    // number of grid entities of a given codim on a given level
    g.size(0,0);
    // number of leaf entities per codim in this process
    g.size(0);
    // number of entities per level and geometry type in this process
    g.size(0, Dune::GeometryType(Dune::GeometryType::cube,Grid::dimension));
    // number of leaf entities per geometry type in this process
    g.size(Dune::GeometryType(Dune::GeometryType::cube,Grid::dimension));

    // Check overlap and ghost size on level 0
    g.overlapSize(0,0);
    g.ghostSize(0,0);

    // Check overlap and ghost size on the leaf level
    g.overlapSize(0);
    g.ghostSize(0);

    // check for iterator functions
    g.template lbegin<0>(0);
    g.template lend<0>(0);

    LeafInterface< Grid>();

    // Check for index sets
    typedef typename Grid::LevelIndexSet LevelIndexSet;
    typedef typename Grid::LeafIndexSet LeafIndexSet;
    typedef typename Grid::LocalIdSet LocalIdSet;
    typedef typename Grid::GlobalIdSet GlobalIdSet;

    g.levelIndexSet(0);
    if (g.template lbegin<0>(0) !=g.template lend<0>(0) ) {
      // Instantiate all methods of LevelIndexSet
      g.levelIndexSet(0).index(*g.template lbegin<0>(0));
      /** \todo Test for subindex is missing, because I don't know yet
          how to test for the existence of certain codims */
    }
    g.levelIndexSet(0).
    size(Dune::GeometryType(Dune::GeometryType::simplex,Grid::dimension));
    for (int codim = 0; codim < Grid::dimension; codim++)
      g.levelIndexSet(0).geomTypes(codim);

    if (g.template leafbegin<0>() != g.template leafend<0>() )
    {
      // Instantiate all methods of LeafIndexSet
      g.leafIndexSet().index(*g.template leafbegin<0>());
    }
    /** \todo Test for subindex is missing, because I don't know yet
       how to test for the existence of certain codims */
    g.leafIndexSet().size(Dune::GeometryType(Dune::GeometryType::simplex,Grid::dimension));
    for (int codim = 0; codim < Grid::dimension; codim++)
      g.leafIndexSet().geomTypes(codim);

    if (g.template lbegin<0>(0) !=g.template lend<0>(0) ) {
      // Instantiate all methods of LocalIdSet
      /** \todo Test for subindex is missing, because I don't know yet
          how to test for the existence of certain codims */
      g.localIdSet().id(*g.template lbegin<0>(0));
      // Instantiate all methods of GlobalIdSet
      /** \todo Test for subindex is missing, because I don't know yet
          how to test for the existence of certain codims */
      g.globalIdSet().id(*g.template lbegin<0>(0));
    }
    // recursively check entity-interface
    // ... we only allow grids with codim 0 zero entites
    dune_static_assert((Dune::Capabilities::hasEntity<Grid, 0>::v),"Grid must have codim 0 entities");
    dune_static_assert((Dune::Capabilities::hasEntity<const Grid, 0>::v),"Grid must have codim 0 entities");

    EntityInterface< Grid, 0, Grid::dimension, Dune::Capabilities::hasEntity< Grid, 0 >::v >();

    // !!! check for parallel grid?
    g.template lbegin<0, Dune::Ghost_Partition>(0);
    g.template lend<0, Dune::Ghost_Partition>(0);
  }
  GridInterface()
  {
    c = check;
  }
  // member just to avoid "unused variable"-warning in constructor
  void (*c)(const Grid&);
};

#endif // #ifndef DUNE_STATIC_CHECK_HH
