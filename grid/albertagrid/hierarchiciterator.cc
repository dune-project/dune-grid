// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_HIERARCHICITERATOR_CC
#define DUNE_ALBERTA_HIERARCHICITERATOR_CC

#include <dune/grid/albertagrid/hierarchiciterator.hh>

namespace Dune
{

  template< class GridImp >
  inline AlbertaGridHierarchicIterator< GridImp >
  ::AlbertaGridHierarchicIterator( const GridImp &grid, int actLevel, int maxLevel )
    : Base( grid, actLevel, true, true ),
      startLevel_( actLevel ),
      maxlevel_( maxLevel )
  {}


  template< class GridImp >
  inline AlbertaGridHierarchicIterator< GridImp >
  ::AlbertaGridHierarchicIterator ( const GridImp &grid,
                                    const ElementInfo &elementInfo,
                                    int actLevel, int maxLevel, bool leafIt )
    : Base( grid, elementInfo.level(), leafIt, false ),
      startLevel_( elementInfo.level() ),
      maxlevel_( maxLevel )
  {
    increment( elementInfo );
  }


  template< class GridImp >
  inline AlbertaGridHierarchicIterator< GridImp >
  ::AlbertaGridHierarchicIterator( const This &other )
    : Base( other ),
      startLevel_( other.startLevel_ ),
      maxlevel_( other.maxlevel_ )
  {}


  template< class GridImp >
  inline typename AlbertaGridHierarchicIterator< GridImp >::This &
  AlbertaGridHierarchicIterator< GridImp >::operator= ( const This &other )
  {
    Base::operator=( other );

    startLevel_ = other.startLevel_;
    maxlevel_ = other.maxlevel_;
    return *this;
  }


  template< class GridImp >
  inline void AlbertaGridHierarchicIterator< GridImp >::increment ()
  {
    // note: since we are not the end iterator, we point to a valid entity
    increment( entityImp().elementInfo_ );
  }

  template< class GridImp >
  inline void AlbertaGridHierarchicIterator< GridImp >
  ::increment ( ElementInfo elementInfo )
  {
    if( (elementInfo.level() >= maxlevel_) || elementInfo.isLeaf() )
    {
      while( (elementInfo.level() > startLevel_) && (elementInfo.indexInFather() == 1) )
        elementInfo = elementInfo.father();
      if( elementInfo.level() > startLevel_ )
        entityImp().setElement( elementInfo.father().child( 1 ), 0 );
      else
        this->done();
    }
    else
      entityImp().setElement( elementInfo.child( 0 ), 0 );
  }

}

#endif
