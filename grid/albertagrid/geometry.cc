// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GEOMETRY_CC
#define DUNE_ALBERTA_GEOMETRY_CC

#include <dune/grid/albertagrid/algebra.hh>
#include <dune/grid/albertagrid/geometry.hh>
#include <dune/grid/albertagrid/refinement.hh>

namespace Dune
{

  // AlbertaGridGeometry
  // -------------------

#if !USE_GENERICGEOMETRY
  template< int mydim, int cdim, class GridImp >
  inline AlbertaGridGeometry< mydim, cdim, GridImp >::AlbertaGridGeometry ()
  {
    invalidate();
  }


  template< int mydim, int cdim, class GridImp >
  inline AlbertaGridGeometry< mydim, cdim, GridImp >
  ::AlbertaGridGeometry ( const This &other )
    : coord_( other.coord_ ),
      jT_( other.jT_ ),
      jTInv_( other.jTInv_ ),
      builtJT_( other.builtJT_ ),
      builtJTInv_( other.builtJTInv_ ),
      calcedDet_( other.calcedDet_ ),
      elDet_( other.elDet_ )
  {}


  template< int mydim, int cdim, class GridImp >
  template< class CoordReader >
  inline AlbertaGridGeometry< mydim, cdim, GridImp >
  ::AlbertaGridGeometry ( const CoordReader &coordReader )
  {
    build( coordReader );
  }


  // print the GeometryInformation
  template <int mydim, int cdim, class GridImp>
  inline void AlbertaGridGeometry<mydim,cdim,GridImp>::print (std::ostream& ss) const
  {
    ss << "AlbertaGridGeometry<" << mydim << "," << cdim << "> = { \n";
    for(int i=0; i<corners(); i++)
    {
      ss << " corner " << i << " = ";
      ss << ((*this)[i]); ss << "\n";
    }
    ss << "} \n";
  }

  template< int mydim, int cdim, class GridImp >
  inline GeometryType AlbertaGridGeometry< mydim, cdim, GridImp >::type () const
  {
    return GeometryType( GeometryType::simplex, mydimension );
  }

  template< int mydim, int cdim, class GridImp >
  inline int AlbertaGridGeometry< mydim, cdim, GridImp >::corners () const
  {
    return numCorners;
  }


  template< int mydim, int cdim, class GridImp >
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::GlobalVector
  AlbertaGridGeometry< mydim, cdim, GridImp >::corner ( const int i ) const
  {
    // for simplices old and new corner numbering is identical
    assert( (i >= 0) && (i < numCorners) );
    return coord_[ i ];
  }


  template< int mydim, int cdim, class GridImp >
  inline const typename AlbertaGridGeometry< mydim, cdim, GridImp >::GlobalVector &
  AlbertaGridGeometry< mydim, cdim, GridImp >::operator[] ( int i ) const
  {
    assert( (i >= 0) && (i < numCorners) );
    return coord_[ i ];
  }


  template< int mydim, int cdim, class GridImp >
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::GlobalVector
  AlbertaGridGeometry< mydim, cdim, GridImp >::global ( const LocalVector &local ) const
  {
    GlobalVector y = coord_[ 0 ];
    jacobianTransposed().umtv( local, y );
    return y;
  }

  //local implementation for mydim < cdim
  template< int mydim, int cdim, class GridImp >
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::LocalVector
  AlbertaGridGeometry< mydim, cdim, GridImp>::local ( const GlobalVector &global ) const
  {
    GlobalVector y = global;
    y -= coord_[ 0 ];

    LocalVector x;
    FMatrixHelp::multAssignTransposed( jacobianInverseTransposed(), y, x );
    return x;
  }


  template< int mydim, int cdim, class GridImp >
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::ctype
  AlbertaGridGeometry< mydim, cdim, GridImp >::elDeterminant () const
  {
    return std::abs( Alberta::determinant( jacobianTransposed() ) );
  }



  template< int mydim, int cdim, class GridImp >
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::ctype
  AlbertaGridGeometry< mydim, cdim, GridImp >::volume () const
  {
    assert( calcedDet_ );
    static const ctype refVolume
      = ctype( 1 ) / ctype( Factorial< mydimension >::factorial );
    return refVolume * elDet_;
  }


  template< int mydim, int cdim, class GridImp >
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::ctype
  AlbertaGridGeometry< mydim, cdim, GridImp >::integrationElement () const
  {
    assert( calcedDet_ );
    return elDet_;
  }


  template< int mydim, int cdim, class GridImp >
  inline const typename AlbertaGridGeometry< mydim, cdim, GridImp >::JacobianTransposed &
  AlbertaGridGeometry< mydim, cdim, GridImp >::jacobianTransposed () const
  {
    if( !builtJT_ )
    {
      const FieldVector< ctype, coorddimension > &origin = coord_[ 0 ];
      for( int i = 0; i < mydimension; ++i )
      {
        jT_[ i ] = coord_[ i+1 ];
        jT_[ i ] -= origin;
      }
      builtJT_ = true;
    }
    return jT_;
  }


  template< int mydim, int cdim, class GridImp >
  inline const typename AlbertaGridGeometry< mydim, cdim, GridImp >::JacobianInverseTransposed &
  AlbertaGridGeometry< mydim, cdim, GridImp >::jacobianInverseTransposed () const
  {
    if( !builtJTInv_ )
    {
      elDet_ = std::abs( Alberta::invert( jacobianTransposed(), jTInv_ ) );
      assert( elDet_ > 1.0e-25 );
      calcedDet_ = true;
      builtJTInv_ = true;
    }
    return jTInv_;
  }


  template< int mydim, int cdim, class GridImp >
  inline void AlbertaGridGeometry< mydim, cdim, GridImp >::invalidate ()
  {
    builtJT_ = false;
    builtJTInv_ = false;
    calcedDet_ = false;
  }


  // built Geometry
  template< int mydim, int cdim, class GridImp >
  template< class CoordReader >
  inline void AlbertaGridGeometry< mydim, cdim, GridImp >
  ::build ( const CoordReader &coordReader )
  {
    builtJT_ = false;
    builtJTInv_ = false;

    for( int i = 0; i <= mydimension; ++i )
      coordReader.coordinate( i, coord_[ i ] );

    elDet_ = (coordReader.hasDeterminant() ? coordReader.determinant() : elDeterminant());
    assert( std::abs( elDet_ ) > 0.0 );
    calcedDet_ = true;
  }
#endif // #if !USE_GENERICGEOMETRY



  // AlbertaGridLocalGeometryProvider
  // --------------------------------

  template< class Grid >
  void AlbertaGridLocalGeometryProvider< Grid >::buildGeometryInFather ()
  {
    typedef MakeableInterfaceObject< LocalElementGeometry > LocalGeoObject;
    typedef typename LocalGeoObject::ImplementationType LocalGeoImp;

    for( int child = 0; child < numChildren; ++child )
    {
      for( int orientation = 0; orientation < 2; ++orientation )
      {
        const GeoInFatherCoordReader coordReader( child, orientation );
        geometryInFather_[ child ][ orientation ]
          = new LocalGeoObject( LocalGeoImp( coordReader ) );
      }
    }
  }


  template< class Grid >
  void AlbertaGridLocalGeometryProvider< Grid >::buildFaceGeometry ()
  {
    typedef MakeableInterfaceObject< LocalFaceGeometry > LocalGeoObject;
    typedef typename LocalGeoObject::ImplementationType LocalGeoImp;

    for( int face = 0; face < numFaces; ++face )
    {
      for( int twist = minFaceTwist; twist <= maxFaceTwist; ++twist )
      {
        const FaceCoordReader coordReader( face, twist );
        faceGeometry_[ face ][ twist - minFaceTwist ] = new LocalGeoObject( LocalGeoImp( coordReader ) );
      }
    }
  }



  // AlbertaGridLocalGeometryProvider::GeoInFatherCoordReader
  // --------------------------------------------------------

  template< class Grid >
  struct AlbertaGridLocalGeometryProvider< Grid >::GeoInFatherCoordReader
  {
    typedef Alberta::Real ctype;

    typedef FieldVector< ctype, dimension > Coordinate;

  private:
    typedef Alberta::GeometryInFather< dimension > GeoInFather;

    const int child_;
    const int orientation_;

  public:
    GeoInFatherCoordReader ( int child, int orientation )
      : child_( child ),
        orientation_( orientation )
    {}

    void coordinate ( int i, Coordinate &x ) const
    {
      const typename GeoInFather::LocalVector &coord
        = GeoInFather::coordinate( child_, orientation_, i );
      for( int j = 0; j < dimension; ++j )
        x[ j ] = coord[ j ];
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



  // AlbertaGridLocalGeometryProvider::FaceCoordReader
  // -------------------------------------------------

  template< class Grid >
  struct AlbertaGridLocalGeometryProvider< Grid >::FaceCoordReader
  {
    typedef Alberta::Real ctype;

    typedef FieldVector< ctype, dimension > Coordinate;

  private:
    int face_;
    int twist_;

  public:
    FaceCoordReader ( const int face, const int twist = 0 )
      : face_( face ),
        twist_( twist )
    {}

    void coordinate ( const int i, Coordinate &x ) const
    {
      const int ti = Alberta::applyInverseTwist< dimension-1 >( twist_, i );
      const int j = mapVertices< 1 >( face_, ti );
      refCorner( j, x );
    }

    bool hasDeterminant () const
    {
      return false;
    }

    ctype determinant () const
    {
      return ctype( 0 );
    }

  private:
    static void refCorner ( const int i, Coordinate &x )
    {
      x = ctype( 0 );
      if( i > 0 )
        x[ i-1 ] = ctype( 1 );
    }
  };

}

#endif // #ifndef DUNE_ALBERTA_GEOMETRY_CC
