// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_DATAHANDLE_HH
#define DUNE_GEOGRID_DATAHANDLE_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/geometrygrid/capabilities.hh>
#include <dune/grid/geometrygrid/entity.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // GeometryGridDataHandle
    // ----------------------

    template< class Grid, class WrappedHandle >
    class CommDataHandle
      : public CommDataHandleIF< CommDataHandle< Grid, WrappedHandle >, typename WrappedHandle::DataType >
    {
      typedef typename std::remove_const< Grid >::type::Traits Traits;

    public:
      CommDataHandle ( const Grid &grid, WrappedHandle &handle )
        : grid_( grid ),
          wrappedHandle_( handle )
      {}

      bool contains ( int dim, int codim ) const
      {
        const bool contains = wrappedHandle_.contains( dim, codim );
        if( contains )
          assertHostEntity( dim, codim );
        return contains;
      }

      bool fixedSize ( int dim, int codim ) const
      {
        return wrappedHandle_.fixedSize( dim, codim );
      }

      template< class HostEntity >
      size_t size ( const HostEntity &hostEntity ) const
      {
        typedef typename Grid::Traits::template Codim< HostEntity::codimension >::Entity Entity;
        typedef typename Grid::Traits::template Codim< HostEntity::codimension >::EntityImpl EntityImpl;
        Entity entity( EntityImpl( grid_, hostEntity ) );
        return wrappedHandle_.size( entity );
      }

      template< class MessageBuffer, class HostEntity >
      void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
      {
        typedef typename Grid::Traits::template Codim< HostEntity::codimension >::Entity Entity;
        typedef typename Grid::Traits::template Codim< HostEntity::codimension >::EntityImpl EntityImpl;
        Entity entity( EntityImpl( grid_, hostEntity ) );
        wrappedHandle_.gather( buffer, entity );
      }

      template< class MessageBuffer, class HostEntity >
      void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size_ )
      {
        typedef typename Grid::Traits::template Codim< HostEntity::codimension >::Entity Entity;
        typedef typename Grid::Traits::template Codim< HostEntity::codimension >::EntityImpl EntityImpl;
        Entity entity( EntityImpl( grid_, hostEntity ) );
        wrappedHandle_.scatter( buffer, entity, size_ );
      }

    private:
      static void assertHostEntity ( int , int codim )
      {
        if( !Capabilities::CodimCache< Grid >::hasHostEntity( codim ) )
          DUNE_THROW( NotImplemented, "Host grid has no entities for codimension " << codim << "." );
      }

      const Grid &grid_;
      WrappedHandle &wrappedHandle_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_DATAHANDLE_HH
