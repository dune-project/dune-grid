// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CHECK_ENTITYSEED_CC
#define DUNE_CHECK_ENTITYSEED_CC

//- C++ includes
#include <cassert>
#include <ostream>

//- dune-common includes
#include <dune/common/forloop.hh>

//- dune-grid includes
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/gridview.hh>

/**
   @file
   @author Christoph Gersbacher
   @brief  Provides a check for EntitySeeds of all codimensions.
 */


namespace CheckEntitySeed // don't blur namespace Dune
{

  // Capability hasEntitySeed
  // ------------------------

  // allow for switching off the tests below by overloading this capability
  template< class Grid, int codim >
  struct hasEntitySeed
  {
    static const bool v = Dune::Capabilities::hasEntity< Grid, codim >::v;
  };



  // Equals
  // ------

  template< class T >
  struct Equals
  {
    typedef T Type;

    static bool apply ( const T &t1, const T &t2 )
    {
      return ( t1 == t2 );
    }
  };


  // Template specialization for const
  // ---------------------------------

  template< class T >
  struct Equals< const T >
    : public Equals< T >
  { };


  // Template specialization for Dune::Geometry
  // ------------------------------------------

  template< int mydim, int cdim, class GridImp, template< int, int, class > class GeometryImp >
  struct Equals< Dune::Geometry< mydim, cdim, GridImp, GeometryImp > >
  {
    typedef typename Dune::Geometry< mydim, cdim, GridImp, GeometryImp > Type;

    static bool apply ( const Type &t1, const Type &t2, const double eps = 1e-10 )
    {
      typedef typename Type::GlobalCoordinate GlobalCooridinate;

      // check geometry type
      if( t1.type() != t2.type() )
        return false;

      // compare corners
      assert( t1.corners() == t2.corners() );
      const int corners = t1.corners();
      for( int i = 0; i < corners; ++i )
      {
        if( (t1.corner( i ) - t2.corner( i )).two_norm() > eps )
          return false;
      }

      return true;
    }
  };



  // Check
  // -----

  template< int codim, class GridView,
      bool check = hasEntitySeed< typename GridView::Grid, codim >::v
      >
  class Check;


  // Template specialization for hasEntitySeed< Grid, codim >::v = false
  // -------------------------------------------------------------------

  template< int codim, class GridView >
  class Check< codim, GridView, false >
  {
  public:
    static void apply ( const GridView &gridView, std::ostream &output )
    { };
  };


  // Template specialization for hasEntitySeed< Grid, codim >::v = true
  // ------------------------------------------------------------------

  template< int codim, class GridView >
  class Check< codim, GridView, true >
  {
  public:
    // entity pointer type
    typedef typename GridView::template Codim< codim >::EntityPointer EntityPointer;
    // iterator type
    typedef typename GridView::template Codim< codim >::Iterator Iterator;
    // entity type
    typedef typename EntityPointer::Entity Entity;
    // geometry type
    typedef typename Entity::Geometry Geometry;

    // grid type
    typedef typename GridView::Grid Grid;
    // type of entity seed (not available on GridView)
    typedef typename Grid::template Codim< codim >::EntitySeed EntitySeed;

    static void apply ( const GridView &gridView, std::ostream &output )
    {
      // get grid, as method entityPointer() is missing on GridViews
      const Grid &grid = gridView.grid();

      const Iterator end = gridView.template end< codim >();
      for( Iterator it = gridView.template begin< codim >(); it != end; ++it )
      {
        // get entity and entity seed
        const Entity &entity = *it;
        EntitySeed seed = entity.seed();

        // get entity pointer from seed
        EntityPointer entityPointer = grid.entityPointer( seed );
        compare( entityPointer, EntityPointer( it ), output );

        // create copy of seed and compare again
        // we might like to check the assignment operator as well
        EntitySeed seed2( seed );
        compare( entityPointer, grid.entityPointer( seed2 ), output );
      }
    }

  private:
    // compare two entity pointers for equality
    static void compare ( const EntityPointer &ep1, const EntityPointer &ep2, std::ostream &output )
    {
      // compare entity pointers
      if( !Equals< EntityPointer >::apply( ep1, ep2 ) )
        output << "Warning: EntityPointers do not conincide" << std::endl;

      // compare geometries
      const double eps = 1e-10;
      if( !Equals< Geometry >::apply( ep1->geometry(), ep2->geometry(), eps ) )
        output << "Warning: Geometries do not conincide" << std::endl;
    }
  };



  // IfHasEntitySeed
  // ---------------

  template< int codim >
  struct IfHasEntitySeed
  {
    template< class GridView >
    static void apply ( const GridView &gridView, std::ostream &output )
    {
      Check< codim, GridView >::apply( gridView, output );
    }
  };

} // namespace CheckEntitySeed



namespace Dune
{

  /** \brief Check for EntitySeeds for all codimensions available in a grid implementation
      \param  gridView  Grid view for which EntitySeeds will be tested
   */
  template< class VT >
  void checkEntitySeed ( const GridView< VT > &gridView, std::ostream &output = std::cerr )
  {
    const int dimension = GridView< VT >::dimension;
    ForLoop< CheckEntitySeed::IfHasEntitySeed, 0, dimension >::apply( gridView, output );
  };

} // namespace Dune

#endif // #ifndef DUNE_CHECK_ENTITYSEED_CC
