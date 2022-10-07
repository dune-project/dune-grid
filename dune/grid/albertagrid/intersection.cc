// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_INTERSECTION_CC
#define DUNE_ALBERTA_INTERSECTION_CC

#include <dune/grid/albertagrid/intersection.hh>

namespace Dune
{

  // AlbertaGridIntersectionBase
  // ---------------------------

  template< class Grid >
  inline AlbertaGridIntersectionBase< Grid >
  ::AlbertaGridIntersectionBase ()
    : grid_( nullptr ),
      elementInfo_(),
      oppVertex_( -1 ) // mark invalid intersection
  {}

  template< class Grid >
  inline AlbertaGridIntersectionBase< Grid >
  ::AlbertaGridIntersectionBase ( const EntityImp &entity, const int oppVertex )
    : grid_( &entity.grid() ),
      elementInfo_( entity.elementInfo() ),
      oppVertex_( oppVertex )
  {}


  template< class Grid >
  inline typename Grid::template Codim< 0 >::Entity
  AlbertaGridIntersectionBase< Grid >::inside () const
  {
    typedef AlbertaGridEntity< 0, Grid::dimension, Grid > EntityImp;
    return EntityImp( grid(), elementInfo(), 0 );
  }


  template< class Grid >
  inline bool AlbertaGridIntersectionBase< Grid >::boundary () const
  {
    return elementInfo().isBoundary( oppVertex_ );
  }


  template< class Grid >
  inline int AlbertaGridIntersectionBase< Grid >::boundaryId () const
  {
    if( boundary() )
    {
      const int id = elementInfo().boundaryId( oppVertex_ );
      assert( id != 0 );
      return id;
    }
    else
      return 0;
  }


  template< class Grid >
  inline size_t AlbertaGridIntersectionBase< Grid >::boundarySegmentIndex () const
  {
    assert( boundary() );
    const Alberta::BasicNodeProjection *projection = elementInfo().boundaryProjection( oppVertex_ );
    assert( projection );
    return projection->boundaryIndex();
  }


  template< class Grid >
  inline int AlbertaGridIntersectionBase< Grid >::indexInInside () const
  {
    const int face = (dimension > 1 ? oppVertex_ : 1-oppVertex_);
    return grid().alberta2generic( 1, face );
  }


  template< class Grid >
  inline GeometryType AlbertaGridIntersectionBase< Grid >::type () const
  {
    return GeometryTypes::simplex( dimension-1 );
  }


  template< class Grid >
  inline typename AlbertaGridIntersectionBase< Grid >::NormalVector
  AlbertaGridIntersectionBase< Grid >::centerIntegrationOuterNormal () const
  {
    const typename Entity::Geometry geoInside = inside().geometry();

    const int face = indexInInside();
    auto refSimplex = ReferenceElements< ctype, dimension >::simplex();
    FieldVector< ctype, dimension > refNormal = refSimplex.integrationOuterNormal( face );

    const typename Entity::Geometry::JacobianInverseTransposed &jInvT
      = geoInside.impl().jacobianInverseTransposed();

    NormalVector normal;
    jInvT.mv( refNormal, normal );
    normal *= geoInside.impl().integrationElement();

    return normal;
  }

  template<>
  inline AlbertaGridIntersectionBase< const AlbertaGrid< 1, 1 > >::NormalVector
  AlbertaGridIntersectionBase< const AlbertaGrid< 1, 1 > >::centerIntegrationOuterNormal () const
  {
    const Alberta::GlobalVector &oppCoord = grid().getCoord( elementInfo(), oppVertex_ );
    const Alberta::GlobalVector &myCoord = grid().getCoord( elementInfo(), 1-oppVertex_ );
    NormalVector n;
    n[ 0 ] = (myCoord[ 0 ] > oppCoord[ 0 ] ? ctype( 1 ) : -ctype( 1 ));
    return n;
  }

#if defined GRIDDIM && GRIDDIM > 1
  template<>
  inline AlbertaGridIntersectionBase< const AlbertaGrid< 2, 2 > >::NormalVector
  AlbertaGridIntersectionBase< const AlbertaGrid< 2, 2 > >::centerIntegrationOuterNormal () const
  {
    const Alberta::GlobalVector &coordOne = grid().getCoord( elementInfo(), (oppVertex_+1)%3 );
    const Alberta::GlobalVector &coordTwo = grid().getCoord( elementInfo(), (oppVertex_+2)%3 );

    NormalVector n;
    n[ 0 ] = -(coordOne[ 1 ] - coordTwo[ 1 ]);
    n[ 1 ] =   coordOne[ 0 ] - coordTwo[ 0 ];
    return n;
  }
#endif // defined GRIDDIM && GRIDDIM > 1

  template<>
  inline AlbertaGridIntersectionBase< const AlbertaGrid< 3, 3 > >::NormalVector
  AlbertaGridIntersectionBase< const AlbertaGrid< 3, 3 > >::centerIntegrationOuterNormal () const
  {
    // in this case the orientation is negative, multiply by -1
    const ALBERTA EL_INFO &elInfo = elementInfo().elInfo();
    const ctype val = (elInfo.orientation > 0) ? 1.0 : -1.0;

    static const int faceVertices[ 4 ][ 3 ]
      = { {1,3,2}, {0,2,3}, {0,3,1}, {0,1,2} };
    const int *localFaces = faceVertices[ oppVertex_ ];

    const Alberta::GlobalVector &coord0 = grid().getCoord( elementInfo(), localFaces[ 0 ] );
    const Alberta::GlobalVector &coord1 = grid().getCoord( elementInfo(), localFaces[ 1 ] );
    const Alberta::GlobalVector &coord2 = grid().getCoord( elementInfo(), localFaces[ 2 ] );

    FieldVector< ctype, dimensionworld > u;
    FieldVector< ctype, dimensionworld > v;
    for( int i = 0; i < dimension; ++i )
    {
      v[ i ] = coord1[ i ] - coord0[ i ];
      u[ i ] = coord2[ i ] - coord1[ i ];
    }

    NormalVector n;
    for( int i = 0; i < dimension; ++i )
    {
      const int j = (i+1)%dimension;
      const int k = (i+2)%dimension;
      n[ i ] = val * (u[ j ] * v[ k ] - u[ k ] * v[ j ]);
    }
    return n;
  }


  template< class Grid >
  inline typename AlbertaGridIntersectionBase< Grid >::NormalVector
  AlbertaGridIntersectionBase< Grid >::centerOuterNormal() const
  {
    return centerIntegrationOuterNormal();
  }


  template< class Grid >
  inline typename AlbertaGridIntersectionBase< Grid >::NormalVector
  AlbertaGridIntersectionBase< Grid >::centerUnitOuterNormal () const
  {
    NormalVector normal = centerOuterNormal();
    normal *= (1.0 / normal.two_norm());
    return normal;
  }


  template< class Grid >
  inline typename AlbertaGridIntersectionBase< Grid >::NormalVector
  AlbertaGridIntersectionBase< Grid >::integrationOuterNormal ( const LocalCoordType & ) const
  {
    return centerIntegrationOuterNormal();
  }


  template< class Grid >
  inline typename AlbertaGridIntersectionBase< Grid >::NormalVector
  AlbertaGridIntersectionBase< Grid >::outerNormal( const LocalCoordType & ) const
  {
    return centerOuterNormal();
  }


  template< class Grid >
  inline typename AlbertaGridIntersectionBase< Grid >::NormalVector
  AlbertaGridIntersectionBase< Grid >::unitOuterNormal ( const LocalCoordType & ) const
  {
    return centerUnitOuterNormal();
  }


  template< class Grid >
  inline AlbertaTransformation
  AlbertaGridIntersectionBase< Grid >::transformation () const
  {
    return AlbertaTransformation( elementInfo().transformation( oppVertex_ ) );
  }


  template< class Grid >
  inline const Grid &AlbertaGridIntersectionBase< Grid >::grid () const
  {
    return *grid_;
  }


  template< class Grid >
  inline const typename AlbertaGridIntersectionBase< Grid >::ElementInfo &
  AlbertaGridIntersectionBase< Grid >::elementInfo () const
  {
    assert( !!elementInfo_ );
    return elementInfo_;
  }



  // AlbertaGridIntersectionBase::GlobalCoordReader
  // ----------------------------------------------

  template< class GridImp >
  struct AlbertaGridIntersectionBase< GridImp >::GlobalCoordReader
  {
    typedef typename std::remove_const< GridImp >::type Grid;

    static const int dimension = Grid::dimension;
    static const int codimension = 1;
    static const int mydimension = dimension - codimension;
    static const int coorddimension = Grid::dimensionworld;

    typedef Alberta::Real ctype;

    typedef Alberta::ElementInfo< dimension > ElementInfo;
    typedef FieldVector< ctype, coorddimension > Coordinate;

  private:
    const Grid &grid_;
    const ElementInfo &elementInfo_;
    const int subEntity_;
    const int twist_;

  public:
    GlobalCoordReader ( const GridImp &grid,
                        const ElementInfo &elementInfo,
                        int subEntity )
      : grid_( grid ),
        elementInfo_( elementInfo ),
        subEntity_( subEntity ),
        twist_( elementInfo.template twist< codimension >( subEntity ) )
    {}

    void coordinate ( int i, Coordinate &x ) const
    {
      assert( !elementInfo_ == false );
      assert( (i >= 0) && (i <= mydimension) );

      const int ti = Alberta::applyInverseTwist< mydimension >( twist_, i );
      const int k = mapVertices( subEntity_, ti );
      const Alberta::GlobalVector &coord = grid_.getCoord( elementInfo_, k );
      for( int j = 0; j < coorddimension; ++j )
        x[ j ] = coord[ j ];
    }

    bool hasDeterminant () const
    {
      return false;
    }

    ctype determinant () const
    {
      assert( false );
      return ctype( 0 );
    }

  private:
    static int mapVertices ( int subEntity, int i )
    {
      return Alberta::MapVertices< dimension, codimension >::apply( subEntity, i );
    }
  };




  // AlbertaGridIntersectionBase::LocalCoordReader
  // ---------------------------------------------

  template< class GridImp >
  struct AlbertaGridIntersectionBase< GridImp >::LocalCoordReader
  {
    typedef typename std::remove_const< GridImp >::type Grid;

    static const int dimension = Grid::dimension;
    static const int codimension = 1;
    static const int mydimension = dimension - codimension;
    static const int coorddimension = dimension;

    typedef Alberta::Real ctype;

    typedef FieldVector< ctype, coorddimension > Coordinate;

    typedef typename Grid::template Codim< 0 >::Geometry ElementGeometry;
    typedef typename Grid::template Codim< 1 >::Geometry FaceGeometry;

  private:
    const ElementGeometry &elementGeometry_;
    const FaceGeometry &faceGeometry_;

  public:
    LocalCoordReader ( const ElementGeometry &elementGeometry,
                       const FaceGeometry &faceGeometry )
      : elementGeometry_( elementGeometry ),
        faceGeometry_( faceGeometry )
    {}

    void coordinate ( int i, Coordinate &x ) const
    {
      x = elementGeometry_.local( faceGeometry_.corner( i ) );
    }

    bool hasDeterminant () const
    {
      return false;
    }

    ctype determinant () const
    {
      return ctype( 0 );
    }
  };



  // AlbertaGridLeafIntersection
  // ---------------------------

  template< class GridImp >
  inline void AlbertaGridLeafIntersection<GridImp>::next ()
  {
    assert( oppVertex_ <= dimension );
    ++oppVertex_;
    neighborInfo_ = ElementInfo();
  }

  template< class GridImp >
  inline typename GridImp::template Codim< 0 >::Entity
  AlbertaGridLeafIntersection< GridImp >::outside () const
  {
    typedef AlbertaGridEntity< 0, GridImp::dimension, GridImp > EntityImp;

    if( !neighborInfo_ )
    {
      assert( neighbor() );

      neighborInfo_ = elementInfo().leafNeighbor( oppVertex_ );
    }

    assert( !neighborInfo_ == false );
    assert( neighborInfo_.el() != NULL );
    return EntityImp( grid(), neighborInfo_, 0 );
  }

  template< class GridImp >
  inline bool AlbertaGridLeafIntersection< GridImp >::neighbor () const
  {
    assert( oppVertex_ <= dimension );
    return elementInfo().hasLeafNeighbor( oppVertex_ );
  }


  template< class GridImp >
  inline typename AlbertaGridLeafIntersection< GridImp >::LocalGeometry
  AlbertaGridLeafIntersection< GridImp >::geometryInInside () const
  {
    typedef AlbertaGridLocalGeometryProvider< GridImp > LocalGeoProvider;
    const int twist = elementInfo().template twist< 1 >( oppVertex_ );
    const int face = (dimension > 1 ? oppVertex_ : 1-oppVertex_);
    return LocalGeometry( LocalGeoProvider::instance().faceGeometry( face, twist ) );
  }


  template< class GridImp >
  inline typename AlbertaGridLeafIntersection< GridImp >::LocalGeometry
  AlbertaGridLeafIntersection< GridImp >::geometryInOutside () const
  {
    assert( neighbor() );

    typedef AlbertaGridLocalGeometryProvider< GridImp > LocalGeoProvider;
    const ALBERTA EL_INFO &elInfo = elementInfo().elInfo();
    const int oppVertex = elInfo.opp_vertex[ oppVertex_ ];
    const int twist = elementInfo().twistInNeighbor( oppVertex_ );
    const int face = (dimension > 1 ? oppVertex : 1-oppVertex);
    return LocalGeometry( LocalGeoProvider::instance().faceGeometry( face, twist ) );
  }


  template< class GridImp >
  inline typename AlbertaGridLeafIntersection< GridImp >::Geometry
  AlbertaGridLeafIntersection< GridImp >::geometry () const
  {
    const int face = (dimension > 1 ? oppVertex_ : 1-oppVertex_);
    const GlobalCoordReader coordReader( grid(), elementInfo(), face );
    return Geometry( GeometryImpl( coordReader ) );
  }


  template< class GridImp >
  inline int AlbertaGridLeafIntersection< GridImp >::indexInOutside () const
  {
    const ALBERTA EL_INFO &elInfo = elementInfo().elInfo();
    const int oppVertex = elInfo.opp_vertex[ oppVertex_ ];
    const int face = (dimension > 1 ? oppVertex : 1-oppVertex);
    return grid().alberta2generic( 1, face );
  }

} // namespace Dune

#endif // #ifndef DUNE_ALBERTA_INTERSECTION_CC
