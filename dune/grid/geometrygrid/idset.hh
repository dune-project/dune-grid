// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
      typedef IdSet< Grid, HostIdSet > This;
      typedef Dune::IdSet< Grid, This, typename HostIdSet::IdType > Base;

      typedef typename std::remove_const< Grid >::type::Traits Traits;

    public:
      typedef typename HostIdSet::IdType IdType;

      using Base::subId;

      IdSet ()
        : hostIdSet_( 0 )
      {}

      explicit IdSet ( const HostIdSet &hostIdSet )
        : hostIdSet_( &hostIdSet )
      {}

      IdSet ( const This &other )
        : hostIdSet_( other.hostIdSet_ )
      {}

      const This &operator= ( const This &other )
      {
        hostIdSet_ = other.hostIdSet_;
        return *this;
      }

      template< int codim >
      IdType id ( const typename Traits::template Codim< codim >::Entity &entity ) const
      {
        return entity.impl().id( hostIdSet() );
      }

      template< class Entity >
      IdType id ( const Entity &entity ) const
      {
        return id< Entity::codimension >( entity );
      }

      IdType subId ( const typename Traits::template Codim< 0 >::Entity &entity, int i, unsigned int codim ) const
      {
        return hostIdSet().subId( Grid::template getHostEntity< 0 >( entity ), i, codim );
      }

      explicit operator bool () const { return bool( hostIdSet_ ); }

    private:
      const HostIdSet &hostIdSet () const
      {
        assert( *this );
        return *hostIdSet_;
      }

      const HostIdSet *hostIdSet_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_IDSET_HH
