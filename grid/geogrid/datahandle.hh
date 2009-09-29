// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_DATAHANDLE_HH
#define DUNE_GEOGRID_DATAHANDLE_HH

#include <dune/common/typetraits.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid/geogrid/capabilities.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class Entity >
  class GeometryGridEntityWrapper;



  // GeometryGridDataHandle
  // ----------------------

  template< class Grid, class WrappedCommDataHandle >
  class GeometryGridCommDataHandle
    : public CommDataHandleIF
      < GeometryGridCommDataHandle< Grid, WrappedCommDataHandle >,
          typename WrappedCommDataHandle :: DataType >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

    const Grid &grid_;
    WrappedCommDataHandle &wrappedHandle_;

  public:
    GeometryGridCommDataHandle ( const Grid &grid,
                                 WrappedCommDataHandle &handle )
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
      const int codimension = HostEntity :: codimension;
      typedef typename Traits :: template Codim< codimension > :: Entity Entity;
      typedef GeometryGridEntityWrapper< Entity > EntityWrapper;
      typedef GeometryGridStorage< EntityWrapper > EntityStorage;

      EntityWrapper *entity = EntityStorage :: alloc();
      entity->initialize( grid_, hostEntity );
      const size_t size = wrappedHandle_.size( (const Entity &)(*entity) );
      EntityStorage :: free( entity );
      return size;
    }

    template< class MessageBuffer, class HostEntity >
    void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
    {
      const int codimension = HostEntity :: codimension;
      typedef typename Traits :: template Codim< codimension > :: Entity Entity;
      typedef GeometryGridEntityWrapper< Entity > EntityWrapper;
      typedef GeometryGridStorage< EntityWrapper > EntityStorage;

      EntityWrapper *entity = EntityStorage :: alloc();
      entity->initialize( grid_, hostEntity );
      wrappedHandle_.gather( buffer, (const Entity &)(*entity) );
      EntityStorage :: free( entity );
    }

    template< class MessageBuffer, class HostEntity >
    void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size )
    {
      const int codimension = HostEntity :: codimension;
      typedef typename Traits :: template Codim< codimension > :: Entity Entity;
      typedef GeometryGridEntityWrapper< Entity > EntityWrapper;
      typedef GeometryGridStorage< EntityWrapper > EntityStorage;

      EntityWrapper *entity = EntityStorage :: alloc();
      entity->initialize( grid_, hostEntity );
      wrappedHandle_.scatter( buffer, (const Entity &)(*entity), size );
      EntityStorage :: free( entity );
    }

  private:
    void assertHostEntity ( int dim, int codim ) const
    {
      if( !Capabilities :: CodimCache< Grid > :: hasHostEntity( codim ) )
      {
        DUNE_THROW( NotImplemented, "Host grid has no entities for codimension "
                    << codim << "." );
      }
    }
  };

}

#endif
