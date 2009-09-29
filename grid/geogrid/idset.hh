// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_IDSET_HH
#define DUNE_GEOGRID_IDSET_HH

#include <dune/grid/common/indexidset.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // IdSet
    // -----

    template< class Grid, class HostIdSet >
    class IdSet
      : public Dune :: IdSet
        < Grid, IdSet< Grid, HostIdSet >, typename HostIdSet :: IdType >
    {
      typedef typename remove_const< Grid > :: type :: Traits Traits;

      const HostIdSet &hostIdSet_;

    public:
      typedef typename HostIdSet :: IdType IdType;

      IdSet ( const HostIdSet &hostIdSet )
        : hostIdSet_( hostIdSet )
      {}

      template< int codim >
      IdType id ( const typename Traits :: template Codim< codim > :: Entity &entity ) const
      {
        return Grid :: getRealImplementation( entity ).id( hostIdSet_ );
      }

      template< class Entity >
      IdType id ( const Entity &entity ) const
      {
        return id< Entity :: codimension >( entity );
      }

      template< int codim >
      IdType subId ( const typename Traits :: template Codim< 0 > :: Entity &entity, int i) const
      {
        return hostIdSet_.template subId< codim >( Grid :: template getHostEntity< 0 >( entity ), i );
      }

    private:
      IdSet ( const IdSet & );
      IdSet &operator= ( const IdSet & );
    };

  }

}

#endif
