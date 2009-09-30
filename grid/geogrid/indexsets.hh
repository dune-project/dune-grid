// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INDEXSETS_HH
#define DUNE_GEOGRID_INDEXSETS_HH

#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/indexidset.hh>

#include <dune/grid/geogrid/capabilities.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class HostGrid, class CoordFunction,
        bool hasHierarchicIndexSet = Capabilities::hasHierarchicIndexSet< HostGrid >::v >
    class HierarchicIndexSetProvider;



    // IndexSet
    // --------

    template< class Grid, class HostIndexSet >
    class IndexSet
      : public Dune::IndexSet< Grid, IndexSet< Grid, HostIndexSet >, typename HostIndexSet::IndexType >
    {
      typedef IndexSet< Grid, HostIndexSet > This;

      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::HostGrid HostGrid;

    public:
      typedef Dune::IndexSet< Grid, This, typename HostIndexSet::IndexType > Base;

      static const int dimension = Grid::dimension;

      typedef typename Base::IndexType IndexType;

    private:
      const HostIndexSet *hostIndexSet_;

    public:
      using Base::subIndex;

      IndexSet ( const HostIndexSet &hostIndexSet )
        : hostIndexSet_( &hostIndexSet )
      {}

      template< int codim >
      IndexType index ( const typename Grid::template Codim< codim >::Entity &entity ) const
      {
        return Grid::getRealImplementation( entity ).index( hostIndexSet() );
      }

      template< class Entity >
      IndexType index ( const Entity &entity ) const
      {
        return index< Entity::codimension >( entity );
      }

      template< int codim >
      IndexType subIndex ( const typename Grid::template Codim< codim >::Entity &entity, int i, unsigned int subcodim ) const
      {
        return Grid::getRealImplementation( entity ).template subIndex( hostIndexSet(), i, subcodim );
      }

      template< class Entity >
      IndexType subIndex ( const Entity &entity, int i, unsigned int subcodim ) const
      {
        return subIndex< Entity::codimension >( entity, i, subcodim );
      }

      // This one is only necessary due to the using directive
      IndexType subIndex ( const typename Grid::template Codim< 0 >::Entity &entity, int i, unsigned int subcodim ) const
      {
        return subIndex< 0 >( entity, i, subcodim );
      }

      IndexType size ( GeometryType type ) const
      {
        return hostIndexSet().size( type );
      }

      int size ( int codim ) const
      {
        return hostIndexSet().size( codim );
      }

      template< int codim >
      bool contains ( const typename Grid::template Codim< codim >::Entity &entity ) const
      {
        return Grid::getRealImplementation( entity ).isContained( hostIndexSet() );
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return contains< Entity::codimension >( entity );
      }

      const std::vector< GeometryType > &geomTypes ( int codim ) const
      {
        return hostIndexSet().geomTypes( codim );
      }

    private:
      const HostIndexSet &hostIndexSet () const
      {
        assert( hostIndexSet_ != 0 );
        return *hostIndexSet_;
      }
    };



    // HierarchicIndexSetProvider
    // --------------------------

    template< class HostGrid, class CoordFunction >
    class HierarchicIndexSetProvider< HostGrid, CoordFunction, false >
    {};

    template< class HostGrid, class CoordFunction >
    class HierarchicIndexSetProvider< HostGrid, CoordFunction, true >
    {
      typedef HierarchicIndexSetProvider< HostGrid, CoordFunction, true > This;

      typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    public:
      typedef IndexSet< const Grid, typename HostGrid::HierarchicIndexSet >
      HierarchicIndexSet;

    private:
      const Grid &grid_;
      mutable HierarchicIndexSet *hierarchicIndexSet_;

    public:
      HierarchicIndexSetProvider ( const Grid &grid )
        : grid_( grid ),
          hierarchicIndexSet_( 0 )
      {}

      HierarchicIndexSetProvider ( const This &other )
        : grid_( other.grid_ ),
          hierarchicIndexSet_( 0 )
      {}

      ~HierarchicIndexSetProvider ()
      {
        if( hierarchicIndexSet_ != 0 )
          delete hierarchicIndexSet_;
      }

      const HierarchicIndexSet &hierarchicIndexSet () const
      {
        if( hierarchicIndexSet_ == 0 )
        {
          hierarchicIndexSet_
            = new HierarchicIndexSet( grid_.hostGrid().hierarchicIndexSet() );
        }
        assert( hierarchicIndexSet_ != 0 );
        return *hierarchicIndexSet_;
      }
    };

  }

}

#endif
