// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_PERSISTENTCONTAINER_HH
#define DUNE_GEOGRID_PERSISTENTCONTAINER_HH

#include <dune/common/typetraits.hh>
#include <dune/grid/utility/persistentcontainer.hh>

namespace Dune
{
  // PersistentContainer for GeometryGrid
  // ------------------------------------

  template< class HostGrid, class CoordFunction, class CoordAllocator, class Data, class Allocator >
  class PersistentContainer< GeometryGrid< HostGrid, CoordFunction, CoordAllocator >, Data, Allocator >
  {
    typedef PersistentContainer< HostGrid, Data, Allocator > HostContainer;

  public:
    typedef GeometryGrid< HostGrid, CoordFunction, CoordAllocator > GridType;

    typedef typename HostContainer::ConstIterator ConstIterator;
    typedef typename HostContainer::Iterator Iterator;

    typedef typename GridType::template Codim< 0 >::Entity ElementType;

    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
      : hostContainer_( grid.hostGrid(), codim, allocator ), grid_(grid)
    {}

    template< class Entity >
    const Data &operator[] ( const Entity &entity ) const
    {
      static const bool fake = ! Capabilities::hasHostEntity< GridType, Entity::codimension >::v;
      integral_constant<bool,fake> a;
      return data(entity,a);
    }

    template< class Entity >
    Data &operator[] ( const Entity &entity )
    {
      static const bool fake = ! Capabilities::hasHostEntity< GridType, Entity::codimension >::v;
      integral_constant<bool,fake> a;
      return data(entity,a);
    }

    const Data &operator() ( const ElementType &element, const int subEntity ) const
    {
      return hostContainer_( GridType::template getHostEntity<0>( element ), subEntity );
    }

    Data &operator() ( const ElementType &element, const int subEntity )
    {
      return hostContainer_( GridType::template getHostEntity<0>( element ), subEntity );
    }

    ConstIterator begin () const { return hostContainer_.begin(); }
    Iterator begin () { return hostContainer_.begin(); }

    ConstIterator end () const { return hostContainer_.end(); }
    Iterator end () { return hostContainer_.end(); }

    void resize ( )
    {
      hostContainer_.resize( );
    }

    void compress ( )
    {
      hostContainer_.compress( );
    }

  protected:
    template< class EntityImpl >
    const Data &data ( const EntityImpl &entity, integral_constant<bool, false > ) const
    {
      return hostContainer_[ GridType::template getHostEntity<EntityImpl::codimension>( entity) ];
    }

    template< class EntityImpl >
    Data &data ( const EntityImpl &entity, integral_constant<bool, false > )
    {
      return hostContainer_[ GridType::template getHostEntity<EntityImpl::codimension>( entity) ];
    }

    template< class EntityImpl >
    const Data &data ( const EntityImpl &entity, integral_constant<bool, true > ) const
    {
      return hostContainer_( entity.hostElement(). entity.subEntity() );
    }

    template< class EntityImpl >
    Data &data ( const EntityImpl &entity, integral_constant<bool, true > )
    {
      return hostContainer_( entity.hostElement(). entity.subEntity() );
    }

  private:
    HostContainer hostContainer_;
    const GridType &grid_;
  };
} // end namespace Dune

#endif // end DUNE_PERSISTENTCONTAINER_HH
