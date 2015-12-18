// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INDEXSETS_HH
#define DUNE_GEOGRID_INDEXSETS_HH

#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/indexidset.hh>

#include <dune/grid/geometrygrid/declaration.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // IndexSet
    // --------

    template< class Grid, class HostIndexSet >
    class IndexSet
      : public Dune::IndexSet< Grid, IndexSet< Grid, HostIndexSet >, typename HostIndexSet::IndexType, typename HostIndexSet::Types >
    {
      typedef IndexSet< Grid, HostIndexSet > This;
      typedef Dune::IndexSet< Grid, This, typename HostIndexSet::IndexType, typename HostIndexSet::Types > Base;

      typedef typename std::remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::HostGrid HostGrid;

    public:
      static const int dimension = Traits::dimension;

      typedef typename Base::IndexType IndexType;

      typedef typename Base::Types Types;

      IndexSet ()
        : hostIndexSet_( 0 )
      {}

      explicit IndexSet ( const HostIndexSet &hostIndexSet )
        : hostIndexSet_( &hostIndexSet )
      {}

      IndexSet ( const This &other )
        : hostIndexSet_( other.hostIndexSet_ )
      {}

      const This &operator= ( const This &other )
      {
        hostIndexSet_ = other.hostIndexSet_;
        return *this;
      }

      using Base::index;
      using Base::subIndex;

      template< int cc >
      IndexType index ( const typename Traits::template Codim< cc >::Entity &entity ) const
      {
        return Grid::getRealImplementation( entity ).index( hostIndexSet() );
      }

      template< int cc >
      IndexType subIndex ( const typename Traits::template Codim< cc >::Entity &entity, int i, unsigned int codim ) const
      {
        return Grid::getRealImplementation( entity ).subIndex( hostIndexSet(), i, codim );
      }

      IndexType size ( GeometryType type ) const
      {
        return hostIndexSet().size( type );
      }

      int size ( int codim ) const
      {
        return hostIndexSet().size( codim );
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return Grid::getRealImplementation( entity ).isContained( hostIndexSet() );
      }

      Types types ( int codim ) const { return hostIndexSet().types( codim ); }

      const std::vector< GeometryType > &geomTypes ( int codim ) const
      {
        return hostIndexSet().geomTypes( codim );
      }

      operator bool () const { return bool( hostIndexSet_ ); }

    private:
      const HostIndexSet &hostIndexSet () const
      {
        assert( *this );
        return *hostIndexSet_;
      }

      const HostIndexSet *hostIndexSet_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_INDEXSETS_HH
