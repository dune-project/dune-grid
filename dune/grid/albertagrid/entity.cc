// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ENTITY_CC
#define DUNE_ALBERTA_ENTITY_CC

namespace Dune
{

  // AlbertaGridEntity (for codim > 0)
  // ---------------------------------

  template<int codim, int dim, class Grid>
  inline AlbertaGridEntity< codim, dim, Grid >
    ::AlbertaGridEntity ( const Grid &grid, const ElementInfo &elementInfo, int subEntity )
  : grid_( &grid ),
    elementInfo_( elementInfo ),
    subEntity_( subEntity )
  {}

  template<int codim, int dim, class Grid>
  inline AlbertaGridEntity< codim, dim, Grid >
    ::AlbertaGridEntity ( const Grid &grid )
  : grid_( &grid ),
    elementInfo_(),
    subEntity_( -1 )
  {}

  template<int codim, int dim, class Grid>
  inline AlbertaGridEntity< codim, dim, Grid >
    ::AlbertaGridEntity ()
  : grid_( NULL ),
    elementInfo_(),
    subEntity_( -1 )
  {}

  template< int codim, int dim, class Grid >
  inline PartitionType
  AlbertaGridEntity< codim,dim,Grid >::partitionType () const
  {
    return InteriorEntity;
  }


  template< int codim, int dim, class Grid >
  inline bool
  AlbertaGridEntity< codim, dim, Grid >::equals ( const This &other ) const
  {
    const Alberta::Element *e1 = elementInfo().el();
    const Alberta::Element *e2 = other.elementInfo().el();

    // if both element null then they are equal
    if( (e1 == NULL) && (e2 == NULL) )
      return true;
    return ((e1 == e2) && (subEntity_ == other.subEntity_));
  }


  template< int codim, int dim, class Grid >
  inline ALBERTA EL_INFO *
  AlbertaGridEntity< codim, dim, Grid >::getElInfo () const
  {
    return &(elementInfo_.elInfo());
  }


  template< int codim, int dim, class Grid >
  inline void
  AlbertaGridEntity< codim, dim, Grid >::clearElement ()
  {
    elementInfo_ = ElementInfo();
  }


  template< int codim, int dim, class Grid >
  inline void AlbertaGridEntity< codim, dim, Grid >
    ::setElement ( const ElementInfo &elementInfo, int subEntity )
  {
    elementInfo_ = elementInfo;
    subEntity_ = subEntity;
  }


  template< int codim, int dim, class Grid >
  inline void
  AlbertaGridEntity< codim, dim, Grid >::setEntity ( const This &other )
  {
    setElement( other.elementInfo_, other.subEntity_ );
  }


  template< int codim, int dim, class Grid >
  inline int AlbertaGridEntity< codim, dim, Grid >::level() const
  {
    assert( elementInfo_.level() == grid().levelProvider() ( elementInfo_ ) );
    return elementInfo_.level();
  }


  template< int codim, int dim, class Grid >
  inline typename AlbertaGridEntity< codim, dim, Grid >::Geometry
  AlbertaGridEntity< codim, dim, Grid >::geometry () const
  {
    typedef AlbertaGridCoordinateReader< codim, Grid > CoordReader;

    assert( elementInfo_ );
    const CoordReader coordReader( grid(), elementInfo_, subEntity_ );
    return Geometry( GeometryImpl( coordReader ) );
  }


  template< int codim, int dim, class Grid >
  inline GeometryType AlbertaGridEntity< codim, dim, Grid >::type () const
  {
    return GeometryTypes::simplex( mydimension );
  }



  // AlbertaGridEntity (for codim = 0)
  // ---------------------------------

  template< int dim, class Grid >
  inline AlbertaGridEntity< 0, dim, Grid >
    ::AlbertaGridEntity ( const Grid &grid, const ElementInfo &elementInfo, int subEntity )
  : grid_( &grid ),
    elementInfo_( elementInfo )
  {
    assert( subEntity == 0 );
  }

  template< int dim, class Grid >
  inline AlbertaGridEntity< 0, dim, Grid >
    ::AlbertaGridEntity( const Grid &grid )
  : grid_( &grid ),
    elementInfo_()
  {}

  template< int dim, class Grid >
  inline AlbertaGridEntity< 0, dim, Grid >
    ::AlbertaGridEntity ()
  : grid_( NULL ),
    elementInfo_()
  {}


  template< int dim, class Grid >
  inline int AlbertaGridEntity< 0, dim, Grid >::boundaryId() const
  {
    // elements are always inside of our Domain
    return 0;
  }


  template< int dim, class Grid >
  inline bool AlbertaGridEntity< 0, dim, Grid >::isNew () const
  {
    return grid().levelProvider().isNew( elementInfo_ );
  }


  template< int dim, class Grid >
  inline bool AlbertaGridEntity< 0, dim, Grid >::mightVanish () const
  {
    return elementInfo_.mightVanish();
  }


  template< int dim, class Grid >
  inline bool
  AlbertaGridEntity< 0, dim, Grid >::hasBoundaryIntersections () const
  {
    assert( elementInfo_ );
    bool isBoundary = false;
    for( int i = 0; i < dim+1; ++i )
      isBoundary |= elementInfo_.isBoundary( i );
    return isBoundary;
  }


  template< int dim, class Grid >
  inline PartitionType
  AlbertaGridEntity< 0, dim, Grid >::partitionType () const
  {
    return InteriorEntity;
  }


  template< int dim, class Grid >
  inline bool AlbertaGridEntity< 0, dim,Grid >::isLeaf () const
  {
    return elementInfo_.isLeaf();
  }


  template< int dim, class Grid >
  inline bool
  AlbertaGridEntity< 0, dim, Grid >::equals ( const This &other ) const
  {
    // element pointers are unique
    return (elementInfo().el() == other.elementInfo().el());
  }


  template< int dim, class Grid >
  template< int codim >
  inline typename Grid::template Codim< codim >::Entity
  AlbertaGridEntity< 0, dim, Grid >::subEntity ( int i ) const
  {
    typedef AlbertaGridEntity< codim, dim, Grid > EntityImpl;
    return EntityImpl( grid(), elementInfo_, grid().generic2alberta( codim, i ) );
  }


  template< int dim, class Grid >
  inline ALBERTA EL_INFO *
  AlbertaGridEntity< 0, dim, Grid >::getElInfo () const
  {
    return &(elementInfo_.elInfo());
  }


  template< int dim, class Grid >
  inline int AlbertaGridEntity< 0, dim, Grid >::level () const
  {
    assert( elementInfo_.level() == grid().levelProvider() ( elementInfo_ ) );
    return elementInfo_.level();
  }


  template< int dim, class Grid >
  inline void
  AlbertaGridEntity< 0, dim, Grid >::clearElement ()
  {
    elementInfo_ = ElementInfo();
  }


  template< int dim, class Grid >
  inline void AlbertaGridEntity< 0, dim, Grid >
  ::setElement ( const ElementInfo &elementInfo, int /* subEntity */ )
  {
    elementInfo_ = elementInfo;
  }


  template< int dim, class Grid >
  inline void
  AlbertaGridEntity< 0, dim, Grid >::setEntity( const This &other )
  {
    setElement( other.elementInfo_, 0 );
  }


  template< int dim, class Grid >
  inline typename AlbertaGridEntity< 0, dim, Grid >::Geometry
  AlbertaGridEntity< 0, dim, Grid >::geometry () const
  {
    typedef AlbertaGridCoordinateReader< 0, Grid > CoordReader;

    assert( elementInfo_ );
    const CoordReader coordReader( grid(), elementInfo_, 0 );
    return Geometry( GeometryImpl( coordReader ) );
  }


  template< int dim, class Grid >
  inline GeometryType AlbertaGridEntity< 0, dim, Grid>::type () const
  {
    return GeometryTypes::simplex( mydimension );
  }


  template< int dim, class Grid >
  inline typename Grid::template Codim< 0 >::Entity
  AlbertaGridEntity< 0, dim, Grid >::father () const
  {
    typedef AlbertaGridEntity< 0, dim, Grid > EntityImpl;

    assert( elementInfo_ );
    const ElementInfo fatherInfo = elementInfo_.father();

    return EntityImpl( grid(), fatherInfo, 0 );
  }


  template< int dim, class Grid >
  inline int AlbertaGridEntity< 0, dim, Grid >::nChild () const
  {
    return elementInfo_.indexInFather();
  }


  template< int dim, class Grid >
  inline typename AlbertaGridEntity< 0, dim, Grid >::LocalGeometry
  AlbertaGridEntity< 0, dim, Grid >::geometryInFather() const
  {
    typedef AlbertaGridLocalGeometryProvider< Grid > LocalGeoProvider;
    const int indexInFather = elementInfo_.indexInFather();
    const int orientation = (elementInfo_.type() == 1 ? -1 : 1);
    return LocalGeometry( LocalGeoProvider::instance().geometryInFather( indexInFather, orientation ) );
  }


  template< int dim, class Grid >
  inline typename AlbertaGridEntity< 0, dim, Grid >::HierarchicIterator
  AlbertaGridEntity< 0, dim, Grid >::hbegin( int maxlevel ) const
  {
    assert( elementInfo_ );
    typedef AlbertaGridHierarchicIterator< Grid > IteratorImp;
    return IteratorImp( grid(), elementInfo_, maxlevel );
  }


  template< int dim, class Grid >
  inline typename AlbertaGridEntity< 0, dim, Grid >::HierarchicIterator
  AlbertaGridEntity< 0, dim, Grid>::hend( int maxlevel ) const
  {
    assert( elementInfo_ );
    typedef AlbertaGridHierarchicIterator< Grid > IteratorImp;
    return IteratorImp( grid(), level(), maxlevel );
  }


  template< int dim, class Grid >
  inline typename AlbertaGridEntity< 0, dim, Grid >::AlbertaGridLeafIntersectionIterator
  AlbertaGridEntity< 0, dim, Grid >::ileafbegin() const
  {
    assert( elementInfo_ );

#ifndef NDEBUG
    for( int i = 0; i <= dimension; ++i )
    {
      // std::cout << "Opposite vertex " << i << ": "
      //           << (int)(getElInfo()->opp_vertex[ i ]) << std::endl;
      if( getElInfo()->opp_vertex[ i ] == 127 )
      {
        assert( false );
        DUNE_THROW( NotImplemented, "AlbertaGrid: Intersections on outside "
                    "entities are not fully implemented, yet." );
      }
    }
#endif // #ifndef NDEBUG

    typename AlbertaGridLeafIntersectionIterator::Begin begin;
    return AlbertaGridLeafIntersectionIterator( *this, begin );
  }


  template< int dim, class Grid >
  inline typename AlbertaGridEntity< 0, dim, Grid >::AlbertaGridLeafIntersectionIterator
  AlbertaGridEntity< 0, dim, Grid >::ileafend() const
  {
    assert( elementInfo_ );
    typename AlbertaGridLeafIntersectionIterator::End end;
    return AlbertaGridLeafIntersectionIterator( *this, end );
  }

} // namespace Dune

#endif // #ifndef DUNE_ALBERTA_ENTITY_CC
