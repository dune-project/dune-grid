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

    template< int codim, class Grid >
    class EntityProxy
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::template Codim< codim >::Entity Entity;
      typedef GeoGrid::Entity< codim, Traits::dimension, const Grid > EntityImpl;
      typedef typename EntityImpl::HostEntity HostEntity;

      template< bool >
      struct InitReal
      {
        static void apply ( EntityImpl &entityImpl, const HostEntity &hostEntity )
        {
          entityImpl.initialize( hostEntity );
        }
      };

      template< bool >
      struct InitFake
      {
        static void apply ( EntityImpl &entityImpl, const HostEntity &hostEntity )
        {
          DUNE_THROW( NotImplemented, "Host grid has no entities for codimension " << codim << "." );
        }
      };

      static const bool hasHostEntity = Capabilities::hasHostEntity< Grid, codim >::v;
      typedef typename SelectType< hasHostEntity, InitReal< true >, InitFake< false > >::Type Init;

    public:
      EntityProxy ( const Grid &grid, const HostEntity &hostEntity )
        : entity_( EntityImpl( grid ) )
      {
        Init::apply( Grid::getRealImplementation( entity_ ), hostEntity );
      }

      const Entity &operator* () const
      {
        return entity_;
      }

    private:
      Entity entity_;
    };



    // GeometryGridDataHandle
    // ----------------------

    template< class Grid, class WrappedHandle >
    class CommDataHandle
      : public CommDataHandleIF< CommDataHandle< Grid, WrappedHandle >, typename WrappedHandle::DataType >
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

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

      bool fixedsize ( int dim, int codim ) const
      {
        return wrappedHandle_.fixedsize( dim, codim );
      }

      template< class HostEntity >
      size_t size ( const HostEntity &hostEntity ) const
      {
        EntityProxy< HostEntity::codimension, Grid > proxy( grid_, hostEntity );
        return wrappedHandle_.size( *proxy );
      }

      template< class MessageBuffer, class HostEntity >
      void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
      {
        EntityProxy< HostEntity::codimension, Grid > proxy( grid_, hostEntity );
        wrappedHandle_.gather( buffer, *proxy );
      }

      template< class MessageBuffer, class HostEntity >
      void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size )
      {
        EntityProxy< HostEntity::codimension, Grid > proxy( grid_, hostEntity );
        wrappedHandle_.scatter( buffer, *proxy, size );
      }

    private:
      static void assertHostEntity ( int dim, int codim )
      {
        if( !Capabilities::CodimCache< Grid >::hasHostEntity( codim ) )
          noEntity( codim );
      }

      static void noEntity ( int codim )
      {
        DUNE_THROW( NotImplemented, "Host grid has no entities for codimension " << codim << "." );
      }

      const Grid &grid_;
      WrappedHandle &wrappedHandle_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_DATAHANDLE_HH
