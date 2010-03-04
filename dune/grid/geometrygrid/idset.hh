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
      : public Dune::IdSet< Grid, IdSet< Grid, HostIdSet >, typename HostIdSet::IdType >
    {
      typedef Dune::IdSet< Grid, IdSet< Grid, HostIdSet >, typename HostIdSet::IdType > Base;

      typedef typename remove_const< Grid >::type::Traits Traits;

    public:
      typedef typename HostIdSet::IdType IdType;

      using Base::subId;

      IdSet ( const HostIdSet &hostIdSet )
        : hostIdSet_( hostIdSet )
      {}

      template< int codim >
      IdType id ( const typename Traits::template Codim< codim >::Entity &entity ) const
      {
        return Grid::getRealImplementation( entity ).id( hostIdSet_ );
      }

      template< class Entity >
      IdType id ( const Entity &entity ) const
      {
        return id< Entity :: codimension >( entity );
      }

      IdType subId ( const typename Traits::template Codim< 0 >::Entity &entity, int i, unsigned int codim ) const
      {
        return hostIdSet_.subId( Grid::template getHostEntity< 0 >( entity ), i, codim );
      }

    private:
      IdSet ( const IdSet & );
      IdSet &operator= ( const IdSet & );

      const HostIdSet &hostIdSet_;
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_IDSET_HH
