// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/alugrid/2d/geometry.hh>
#include <dune/grid/alugrid/2d/intersection.hh>

namespace Dune
{

  // Implementation of ALU2DIntersectionGeometryStorage
  // --------------------------------------------------

  template< class LocalGeometry, class LocalGeometryImp >
  ALU2DIntersectionGeometryStorage< LocalGeometry, LocalGeometryImp >::ALU2DIntersectionGeometryStorage ()
  {
    for( int i = 0; i < 4; ++i )
    {
      for( int j = 0; j < 2; ++j )
      {
        LocalGeometryImp geo;
        if ( i < 3 )
        {
          // build geometry
          geo.buildLocalGeometry( i, j, 3 );
          // create dune geoemtry
          geoms_[ 0 ][ i ][ j ] = new LocalGeometry( geo );
        }
        else
          geoms_[ 0 ][ 3 ][ j ] = 0;

        // build geometry
        geo.buildLocalGeometry( i, j, 4 );
        // create dune geoemtry
        geoms_[ 1 ][ i ][ j ] = new LocalGeometry( geo );
      }
    }
  }


  template< class LocalGeometry, class LocalGeometryImp >
  ALU2DIntersectionGeometryStorage< LocalGeometry, LocalGeometryImp >::~ALU2DIntersectionGeometryStorage ()
  {
    for( int k = 0; k < 2; ++k )
      for( int i = 0; i < 4; ++i )
        for( int j = 0; j < 2; ++j )
          delete geoms_[ k ][ i ][ j ];
  }



  // Implementation of ALU2dGridLevelIntersectionIterator
  // ----------------------------------------------------

  template< class GridImp >
  inline int ALU2dGridLevelIntersectionIterator< GridImp >
  ::getOppositeInFather ( const int nrInChild, const int nrOfChild )
  {
    int ret = (nrInChild==0) ? (2-nrOfChild)
              : ((nrInChild-nrOfChild==2 || nrInChild-nrOfChild==0) ? -1 : 0);
    assert(ret >= -1 && ret < 3);
    return ret;
  }


  template< class GridImp >
  inline int ALU2dGridLevelIntersectionIterator< GridImp >
  ::getOppositeInChild ( const int nrInFather, const int nrOfChild )
  {
    int ret = (nrInFather==0) ? (nrOfChild+1) : ((nrInFather-nrOfChild==1) ? -1 : 0);
    assert( ret >= -1 && ret < 3 );
    return ret;
  }


  template< class GridImp >
  void ALU2dGridLevelIntersectionIterator< GridImp >::addNeighboursToStack ()
  {
    assert( current.index_ < current.nFaces() );

    ThinelementType *neighbor = current.inside()->neighbour( current.index_ );
    assert( neighbor );

    IntersectionInfo info;
    if( neighbor->thinis( ThinelementType::bndel_like ) )
    {
      HBndElType *bndel = (HBndElType *)neighbor;
      if( bndel->type() != HBndElType::periodic )
        return;

      PeriodicBndElType *bndnb = ((PeriodicBndElType *)bndel)->periodic_nb;
      assert( bndnb && bndnb->neighbour( 0 ) && bndnb->neighbour( 0 )->thinis( ThinelementType::element_like ) );
      info.first = (HElementType *)bndnb->neighbour( 0 );
      info.second = bndnb->opposite( 0 );
    }
    else
    {
      assert( neighbor->thinis( ThinelementType::element_like ) );
      info.first = (HElementType *)neighbor;
      info.second = current.inside()->opposite( current.index_ );
    }
    assert( info.first );

    while( info.first->level() > walkLevel_ )
    {
      info.second = getOppositeInFather( info.second, info.first->childNr() );
      assert( (info.second >= 0) && (info.second < current.nFaces()) );
      info.first = info.first->father();
    }

    if( info.first->level() >= walkLevel_ )
    {
      nbStack_.push( info );
      return;
    }

    // why should we go up, here?
    while( info.first && (info.first->level() < walkLevel_ - 1) )
    {
      info.second = getOppositeInFather( info.second, info.first->childNr() );
      assert( (info.second >= 0) && (info.second < current.nFaces()) );
      info.first = info.first->father();
    }

    if( info.first )
    {
      assert( info.first->level() == walkLevel_ - 1 );

      const int opposite = info.second;
      for( info.first = info.first->down(); info.first; info.first = info.first->next() )
      {
        info.second = getOppositeInChild( opposite, info.first->childNr() );
        if( info.second != -1 )
          nbStack_.push( info );
      }
    }
  }


  template< class GridImp >
  void ALU2dGridLevelIntersectionIterator< GridImp >::doIncrement ()
  {
    assert( current.index_ < current.nFaces() );

    this->unsetUp2Date();

    if( nbStack_.empty() )
    {
      ++current.index_;
      if( current.index_ >= current.nFaces() )
      {
        assert( current.index_ == current.nFaces() );
        return;
      }

      addNeighboursToStack();
      // if more then one element in stack we have non-conform intersection
      current.useOutside_ = (nbStack_.size() > 1);

      if( nbStack_.empty() )
        return current.setOutside( 0, -1 );
    }

    setupIntersection();

    assert( !current.outside() || (current.outside()->level() == walkLevel_) );
  }


  template< class GridImp >
  void ALU2dGridLevelIntersectionIterator< GridImp >::setFirstItem ( const HElementType &elem, int wLevel )
  {
    // empty stack first
    while( !nbStack_.empty() )
      nbStack_.pop();

    current.setInside( const_cast< HElementType * >( &elem ) );
    current.index_ = -1;
    current.setOutside( 0, -1 );

    walkLevel_ = wLevel;

    assert( current.inside() );

    increment();
  }


  template< class GridImp >
  inline void ALU2dGridLevelIntersectionIterator< GridImp >::setupIntersection ()
  {
    assert( !nbStack_.empty() );

    IntersectionInfo &info = nbStack_.top();
    current.setOutside( info.first, info.second );
    nbStack_.pop();
  }

}



// Template Instantiation
// ----------------------

template class Dune::ALU2DIntersectionGeometryStorage
< Dune::Geometry< 1, 2, const Dune::ALU2dGrid< 2, 2, ALU2DSPACE triangle >, Dune::ALU2dGridGeometry >,
    Dune::ALU2dGridGeometry< 1, 2, const Dune::ALU2dGrid< 2, 2, ALU2DSPACE triangle > > >;
#ifdef ALUGRID_SURFACE_2D
template class Dune::ALU2DIntersectionGeometryStorage
< Dune::Geometry< 1, 2, const Dune::ALU2dGrid< 2, 2, ALU2DSPACE quadrilateral >, Dune::ALU2dGridGeometry >,
    Dune::ALU2dGridGeometry< 1, 2, const Dune::ALU2dGrid< 2, 2, ALU2DSPACE quadrilateral > > >;

template class Dune::ALU2DIntersectionGeometryStorage
< Dune::Geometry< 1, 2, const Dune::ALU2dGrid< 2, 3, ALU2DSPACE triangle >, Dune::ALU2dGridGeometry >,
    Dune::ALU2dGridGeometry< 1, 2, const Dune::ALU2dGrid< 2, 3, ALU2DSPACE triangle > > >;
template class Dune::ALU2DIntersectionGeometryStorage
< Dune::Geometry< 1, 2, const Dune::ALU2dGrid< 2, 3, ALU2DSPACE quadrilateral >, Dune::ALU2dGridGeometry >,
    Dune::ALU2dGridGeometry< 1, 2, const Dune::ALU2dGrid< 2, 3, ALU2DSPACE quadrilateral > > >;
#endif // #ifdef ALUGRID_SURFACE_2D


template class Dune::ALU2dGridLevelIntersectionIterator< const Dune::ALU2dGrid< 2, 2, ALU2DSPACE triangle > >;
#ifdef ALUGRID_SURFACE_2D
template class Dune::ALU2dGridLevelIntersectionIterator< const Dune::ALU2dGrid< 2, 2, ALU2DSPACE quadrilateral > >;

template class Dune::ALU2dGridLevelIntersectionIterator< const Dune::ALU2dGrid< 2, 3, ALU2DSPACE triangle > >;
template class Dune::ALU2dGridLevelIntersectionIterator< const Dune::ALU2dGrid< 2, 3, ALU2DSPACE quadrilateral > >;
#endif // #ifdef ALUGRID_SURFACE_2D
