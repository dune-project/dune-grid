// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ENTITYPOINTER_CC
#define DUNE_ALBERTA_ENTITYPOINTER_CC

#include <dune/grid/albertagrid/entitypointer.hh>

namespace Dune
{

  template< int codim, class GridImp >
  inline AlbertaGridEntityPointer< codim, GridImp >
  ::AlbertaGridEntityPointer ( const GridImp &grid,
                               int level,
                               const ElementInfo &elementInfo,
                               int subEntity )
    : entity_( grid.template getNewEntity< codim >() )
  {
    assert( entity_ );
    entityImp().setElement( elementInfo, subEntity );
  }


  template<int codim, class GridImp >
  inline AlbertaGridEntityPointer< codim, GridImp >
  ::AlbertaGridEntityPointer ( const GridImp &grid, int level, bool end )
    : entity_( grid.template getNewEntity< codim >() )
  {
    if( end )
      done();
  }

  template< int codim, class GridImp >
  inline AlbertaGridEntityPointer< codim, GridImp >
  ::AlbertaGridEntityPointer ( const EntityImp &entity )
    : entity_( entity.grid().template getNewEntity< codim >() )
  {
    entityImp().setEntity( entity );
  }


  template< int codim, class GridImp >
  inline AlbertaGridEntityPointer< codim, GridImp >
  ::AlbertaGridEntityPointer ( const GridImp &grid, const EntityImp  &entity )
    : entity_( grid.template getNewEntity< codim >() )
  {
    entityImp().setEntity( entity );
  }


  template< int codim, class GridImp >
  inline AlbertaGridEntityPointer< codim, GridImp >
  ::AlbertaGridEntityPointer ( const This &other )
    : entity_( other.grid().template getNewEntity< codim >() )
  {
    entityImp().setEntity( other.entityImp() );
  }


  template< int codim, class GridImp >
  inline typename AlbertaGridEntityPointer< codim, GridImp >::This &
  AlbertaGridEntityPointer< codim, GridImp >::operator= ( const This &other )
  {
    entityImp().setEntity( other.entityImp() );
    return *this;
  }


  template< int codim, class GridImp >
  inline typename AlbertaGridEntityPointer<codim,GridImp>::EntityImp &
  AlbertaGridEntityPointer< codim, GridImp >::entityImp ()
  {
    assert( entity_ != 0 );
    return GridImp::getRealImplementation( *entity_ );
  }


  template< int codim, class GridImp >
  inline const typename AlbertaGridEntityPointer< codim, GridImp >::EntityImp &
  AlbertaGridEntityPointer< codim, GridImp >::entityImp () const
  {
    assert( entity_ != 0 );
    return GridImp::getRealImplementation( *entity_ );
  }


  template< int codim, class GridImp >
  inline const GridImp &AlbertaGridEntityPointer< codim, GridImp >::grid () const
  {
    return entityImp().grid();
  }



  template<int codim, class GridImp >
  inline AlbertaGridEntityPointer< codim, GridImp >::~AlbertaGridEntityPointer ()
  {
    done();
    grid().template freeEntity< codim >( entity_ );
    entity_ = 0;
  }


  template< int codim, class GridImp >
  inline void AlbertaGridEntityPointer< codim, GridImp >::done ()
  {
    entityImp().clearElement();
  }


  template<int codim, class GridImp >
  inline bool
  AlbertaGridEntityPointer< codim, GridImp >::equals ( const This &other ) const
  {
    return entityImp().equals( other.entityImp() );
  }


  template<int codim, class GridImp >
  inline typename AlbertaGridEntityPointer< codim, GridImp >::Entity &
  AlbertaGridEntityPointer< codim, GridImp >::dereference () const
  {
    assert( entity_ != 0 );
    return *entity_;
  }


  template< int codim, class GridImp >
  inline int AlbertaGridEntityPointer< codim, GridImp >::level () const
  {
    return entityImp().level();
  }

}

#endif
