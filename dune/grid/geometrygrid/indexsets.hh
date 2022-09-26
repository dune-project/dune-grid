// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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

      IndexSet () = default;

      explicit IndexSet ( const HostIndexSet &hostIndexSet )
        : hostIndexSet_( &hostIndexSet )
      {}

      // The index set contains a pointer to the host index set, so copying or assigning this can be dangerous.
      IndexSet ( const This & ) = delete;
      IndexSet ( This && ) = delete;

      IndexSet &operator= ( const This & ) = delete;
      IndexSet &operator= ( This && ) = delete;

      using Base::index;
      using Base::subIndex;

      template< int cc >
      IndexType index ( const typename Traits::template Codim< cc >::Entity &entity ) const
      {
        return entity.impl().index( hostIndexSet() );
      }

      template< int cc >
      IndexType subIndex ( const typename Traits::template Codim< cc >::Entity &entity, int i, unsigned int codim ) const
      {
        return entity.impl().subIndex( hostIndexSet(), i, codim );
      }

      std::size_t size ( GeometryType type ) const
      {
        return hostIndexSet().size( type );
      }

      std::size_t size ( int codim ) const
      {
        return hostIndexSet().size( codim );
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return entity.impl().isContained( hostIndexSet() );
      }

      Types types ( int codim ) const { return hostIndexSet().types( codim ); }

      const std::vector< GeometryType > &geomTypes ( int codim ) const
      {
        return hostIndexSet().geomTypes( codim );
      }

      explicit operator bool () const { return bool( hostIndexSet_ ); }

      void reset () { hostIndexSet_ = nullptr; }
      void reset ( const HostIndexSet &hostIndexSet ) { hostIndexSet_ = &hostIndexSet; }

    private:
      const HostIndexSet &hostIndexSet () const
      {
        assert( *this );
        return *hostIndexSet_;
      }

      const HostIndexSet *hostIndexSet_ = nullptr;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_INDEXSETS_HH
