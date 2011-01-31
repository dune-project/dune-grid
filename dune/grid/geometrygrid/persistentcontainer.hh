// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_PERSISTENTCONTAINER_HH
#define DUNE_GEOGRID_PERSISTENTCONTAINER_HH

#include <dune/common/typetraits.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/geometrygrid.hh>

namespace Dune
{
  // PersistentContainer for GeometryGrid
  // ------------------------------------

  template< class HostGrid, class CoordFunction, class Data, class Allocator >
  class PersistentContainer< GeometryGrid< HostGrid, CoordFunction >, Data, Allocator >
    : public PersistentContainer< HostGrid, Data, Allocator >
  {
    typedef PersistentContainer< HostGrid, Data, Allocator > Base;

  public:
    typedef GeometryGrid< HostGrid, CoordFunction > GridType;

    typedef typename GridType::template Codim< 0 >::Entity ElementType;

    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
      : Base( grid.hostGrid(), codim, allocator ),
        grid_(grid)
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
      return Base::operator()( GridType::template getHostEntity<0>( element ), subEntity );
    }

    Data &operator() ( const ElementType &element, const int subEntity )
    {
      return Base::operator()( GridType::template getHostEntity<0>( element ), subEntity );
    }

  protected:
    template< class EntityImpl >
    const Data &data ( const EntityImpl &entity, integral_constant<bool, false > ) const
    {
      return Base::operator[]( GridType::template getHostEntity<EntityImpl::codimension>( entity) );
    }

    template< class EntityImpl >
    Data &data ( const EntityImpl &entity, integral_constant<bool, false > )
    {
      return Base::operator[]( GridType::template getHostEntity<EntityImpl::codimension>( entity) );
    }

    template< class EntityImpl >
    const Data &data ( const EntityImpl &entity, integral_constant<bool, true > ) const
    {
      return Base::operator()( entity.hostElement(). entity.subEntity() );
    }

    template< class EntityImpl >
    Data &data ( const EntityImpl &entity, integral_constant<bool, true > )
    {
      return Base::operator()( entity.hostElement(). entity.subEntity() );
    }

  private:
    const GridType &grid_;
  };
} // end namespace Dune

#endif // end DUNE_PERSISTENTCONTAINER_HH
