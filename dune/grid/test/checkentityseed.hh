// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CHECK_ENTITYSEED_HH
#define DUNE_GRID_CHECK_ENTITYSEED_HH

//- C++ includes
#include <cassert>
#include <ostream>
#include <utility>

//- dune-common includes
#include <dune/common/hybridutilities.hh>

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


  // Equals for Dune::Geometry
  // ------------------------------------------

  template<class Geometry>
  struct GeometryEquals
  {

    static bool apply ( const Geometry &t1, const Geometry &t2, const double eps = 1e-10 )
    {
      //typedef typename Type::GlobalCoordinate GlobalCooridinate;

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
    static void apply ( const GridView &, std::ostream & )
    { };
  };


  // Template specialization for hasEntitySeed< Grid, codim >::v = true
  // ------------------------------------------------------------------

  template< int codim, class GridView >
  class Check< codim, GridView, true >
  {
  public:
    // iterator type
    typedef typename GridView::template Codim< codim >::Iterator Iterator;
    // entity type
    typedef typename GridView::template Codim< codim >::Entity Entity;
    // geometry type
    typedef typename Entity::Geometry Geometry;

    // grid type
    typedef typename GridView::Grid Grid;
    // type of entity seed (not available on GridView)
    typedef typename Grid::template Codim< codim >::EntitySeed EntitySeed;

    // Check whether EntitySeed reports the correct codimension
    static_assert(EntitySeed::codimension==codim, "Codimension exported by EntitySeed is incorrect!");

    static void apply ( const GridView &gridView, std::ostream &output )
    {
      // get grid, as method entity() is missing on GridViews
      const Grid &grid = gridView.grid();

      const Iterator end = gridView.template end< codim >();
      for( Iterator it = gridView.template begin< codim >(); it != end; ++it )
      {
        // get entity
        const Entity &entity = *it;
        EntitySeed seed = entity.seed();

        // get entity from seed
        Entity entity2 = grid.entity( seed );
        compare( entity, *it, output );

        // test default constructor
        EntitySeed seed2;
        assert(! seed2.isValid());

        // create copy of seed and compare again
        seed2 = seed;
        assert( seed2.isValid());

        // we might like to check the assignment operator as well
        compare( entity2, grid.entity( seed2 ), output );
      }
    }

  private:

    // compare two entities for equality
    static void compare ( const Entity &e1, const Entity &e2, std::ostream &output )
    {
      // compare entities
      if( !Equals< Entity >::apply( e1, e2 ) )
        output << "Warning: Entities do not conincide" << std::endl;

      // compare geometries
      const double eps = 1e-10;
      if( !GeometryEquals< Geometry >::apply( e1.geometry(), e2.geometry(), eps ) )
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
      if constexpr (Dune::Capabilities::hasEntityIterator<typename GridView::Grid, codim>::v)
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
    Hybrid::forEach( std::make_index_sequence< dimension+1 >{},
      [ & ]( auto i ){ CheckEntitySeed::IfHasEntitySeed< i >::apply( gridView, output ); } );
  }

} // namespace Dune

#endif // #ifndef DUNE_GRID_CHECK_ENTITYSEED_HH
