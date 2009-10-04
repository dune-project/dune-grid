// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_DATAHANDLE_HH
#define DUNE_GEOGRID_DATAHANDLE_HH

#include <dune/common/typetraits.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid/geogrid/capabilities.hh>
#include <dune/grid/geogrid/entity.hh>
#include <dune/grid/geogrid/storage.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // External Forward Declarations
    // -----------------------------

    template< class Entity >
    class GeometryGridEntityWrapper;



    // GeometryGridDataHandle
    // ----------------------

    template< class Grid, class WrappedHandle >
    class CommDataHandle
      : public CommDataHandleIF< CommDataHandle< Grid, WrappedHandle >, typename WrappedHandle::DataType >
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

      template< class HostEntity >
      class EntityProxy;

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
        EntityProxy< HostEntity > proxy( grid_, hostEntity );
        return wrappedHandle_.size( *proxy );
      }

      template< class MessageBuffer, class HostEntity >
      void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
      {
        EntityProxy< HostEntity > proxy( grid_, hostEntity );
        wrappedHandle_.gather( buffer, *proxy );
      }

      template< class MessageBuffer, class HostEntity >
      void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size )
      {
        EntityProxy< HostEntity > proxy( grid_, hostEntity );
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
        DUNE_THROW( NotImplemented, "Host grid has no entities for codimension "
                    << codim << "." );
      }

      const Grid &grid_;
      WrappedHandle &wrappedHandle_;
    };



    template< class Grid, class WrappedHandle >
    template< class HostEntity >
    class CommDataHandle< Grid, WrappedHandle >::EntityProxy
    {
      static const int codimension = HostEntity::codimension;
      typedef typename Traits::template Codim< codimension >::Entity Entity;
      typedef GeoGrid::EntityWrapper< Entity > EntityWrapper;
      typedef GeoGrid::Storage< EntityWrapper > EntityStorage;

      template< bool >
      struct InitializeReal
      {
        static void
        apply ( EntityWrapper &entity, const Grid &grid, const HostEntity &hostEntity )
        {
          entity.initialize( grid, hostEntity );
        }
      };

      template< bool >
      struct InitializeFake
      {
        static void
        apply ( EntityWrapper &entity, const Grid &grid, const HostEntity &hostEntity )
        {
          noEntity( codimension );
        }
      };

      static const bool hasHostEntity = Capabilities::hasHostEntity< Grid, codimension >::v;
      typedef typename SelectType< hasHostEntity, InitializeReal<true>, InitializeFake<false> >::Type Initialize;

    public:
      EntityProxy ( const Grid &grid, const HostEntity &hostEntity )
        : entity_( EntityStorage :: alloc() )
      {
        Initialize::apply( *entity_, grid, hostEntity );
      }

      ~EntityProxy ()
      {
        EntityStorage::free( entity_ );
      }

      const Entity &operator* () const
      {
        return *entity_;
      }

    private:
      EntityWrapper *entity_;
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_DATAHANDLE_HH
