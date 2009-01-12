// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GEOMETRY_CC
#define DUNE_ALBERTA_GEOMETRY_CC

#include <dune/grid/albertagrid/geometry.hh>
#include <dune/grid/albertagrid/refinement.hh>

namespace Dune
{

  namespace Alberta
  {

    template< class ctype, int cdim, int mydim >
    struct MatrixOps
    {
      static ctype invert ( const FieldMatrix< ctype, cdim, mydim > &A,
                            FieldMatrix< ctype, cdim, mydim > &invAT )
      {
        // calc A^T*A
        FieldMatrix< ctype, mydim, mydim > ATA;
        FMatrixHelp::multTransposedMatrix( A, ATA );

        // calc invA = A (A^T*A)^-1
        FieldMatrix< ctype, mydim, mydim > invATA;
        const ctype det = FMatrixHelp::invertMatrix( ATA, invATA );
        FMatrixHelp::multMatrix( A, invATA, invAT );
        return sqrt( det );
      }
    };

    template< class ctype, int cdim >
    struct MatrixOps< ctype, cdim, cdim >
    {
      static ctype invert ( const FieldMatrix< ctype, cdim, cdim > &A,
                            FieldMatrix< ctype, cdim, cdim > &invAT )
      {
        const ctype det = FMatrixHelp::invertMatrix_retTransposed( A, invAT );
        return std::abs( det );
      }
    };

  }



  // AlbertaGridGeometry
  // -------------------

  template< int mydim, int cdim, class GridImp >
  inline AlbertaGridGeometry< mydim, cdim, GridImp > :: AlbertaGridGeometry ()
  {
    invalidate();
  }


  template< int mydim, int cdim, class GridImp >
  inline AlbertaGridGeometry< mydim, cdim, GridImp >
  ::AlbertaGridGeometry ( const CoordMatrix &coords )
    : coord_( coords ),
      builtElMat_( false ),
      builtinverse_( false )
  {
    elDet_ = elDeterminant ();
    calcedDet_ = true;

    if( !(elDet_ > 0.0) )
      DUNE_THROW( AlbertaError, "Degenerate Geometry.");
  }


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
  inline const typename AlbertaGridGeometry< mydim, cdim, GridImp >::GlobalVector &
  AlbertaGridGeometry< mydim, cdim, GridImp >::operator[] ( int i ) const
  {
    assert( (i >= 0) && (i < numCorners) );
    return coord_[ i ];
  }


  template <int mydim, int cdim, class GridImp>
  inline void AlbertaGridGeometry<mydim,cdim,GridImp>::calcElMatrix () const
  {
    if( (mydimension == 0) || builtElMat_ )
      return;

    for( int i = 0; i < coorddimension; ++i )
    {
      for( int j = 0; j < mydimension; ++j )
        elMat_[ i ][ j ] = coord_[ j+1 ][ i ] - coord_[ 0 ][ i ];
    }
    builtElMat_ = true;
  }


  template< int mydim, int cdim, class GridImp >
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::GlobalVector
  AlbertaGridGeometry< mydim, cdim, GridImp >::global ( const LocalVector &local ) const
  {
    calcElMatrix();

    GlobalVector y = coord_[ 0 ];
    elMat_.umv( local, y );
    return y;
  }

  //local implementation for mydim < cdim
  template< int mydim, int cdim, class GridImp >
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::LocalVector
  AlbertaGridGeometry< mydim, cdim, GridImp>::local ( const GlobalVector &global ) const
  {
    if( !builtinverse_ )
      buildJacobianInverseTransposed();

    FieldVector< ctype, coorddimension > y = global;
    y -= coord_[ 0 ];

    FieldVector< ctype, mydimension > x;
    FMatrixHelp::multAssignTransposed( Jinv_, y, x );
    return x;
  }

  // determinant of one Geometry, here line
  template<>
  inline AlbertaGridGeometry< 1, 2, const AlbertaGrid< 2, 2 > >::ctype
  AlbertaGridGeometry< 1, 2, const AlbertaGrid< 2, 2 > >::elDeterminant () const
  {
    // volume is length of edge
    FieldVector< ctype, coorddimension > z = coord_[ 0 ] - coord_[ 1 ];
    return z.two_norm();
  }

  // determinant of one Geometry, here line
  template <>
  inline AlbertaGridGeometry< 1, 3, const AlbertaGrid< 3, 3 > >::ctype
  AlbertaGridGeometry<1,3,const AlbertaGrid<3,3> >::elDeterminant () const
  {
    // volume is length of edge
    FieldVector< ctype, coorddimension > z = coord_[0] - coord_[1];
    return z.two_norm();
  }

  // determinant of one Geometry, here triangle
  template <>
  inline AlbertaGridGeometry< 2, 2, const AlbertaGrid< 2, 2 > >::ctype
  AlbertaGridGeometry<2,2,const AlbertaGrid<2,2> >::elDeterminant () const
  {
    calcElMatrix();
    return std::abs ( elMat_.determinant () );
  }

  // determinant of one Geometry, here triangle in 3d
  template <>
  inline AlbertaGridGeometry< 2, 3, const AlbertaGrid< 3, 3 > >::ctype
  AlbertaGridGeometry<2,3,const AlbertaGrid<3,3> >::elDeterminant () const
  {
    enum { dim = 3 };

    // create vectors of face
    FieldVector< ctype, coorddimension > v = coord_[1] - coord_[0];
    FieldVector< ctype, coorddimension > u = coord_[2] - coord_[1];

    // calculate scaled outer normal
    FieldVector< ctype, coorddimension > z;
    for( int i = 0; i < dim; ++i )
      z[i] = u[(i+1)%dim] * v[(i+2)%dim] - u[(i+2)%dim] * v[(i+1)%dim];

    // z is the same as 2.0 * the outer normal
    return z.two_norm();
  }

  // volume of one Geometry, here therahedron
  template <>
  inline AlbertaGridGeometry< 3, 3, const AlbertaGrid< 3, 3 > >::ctype
  AlbertaGridGeometry<3,3,const AlbertaGrid<3,3> >::elDeterminant () const
  {
    calcElMatrix();
    return std::abs(elMat_.determinant ());
  }

  // volume of one Geometry, here point
  template <>
  inline AlbertaGridGeometry< 0, 2, const AlbertaGrid< 2, 2 > >::ctype
  AlbertaGridGeometry<0,2,const AlbertaGrid<2,2> >::elDeterminant () const
  {
    return 1.0;
  }
  // volume of one Geometry, here point
  template <>
  inline AlbertaGridGeometry< 0, 3, const AlbertaGrid< 3, 3 > >::ctype
  AlbertaGridGeometry<0,3,const AlbertaGrid<3,3> >::elDeterminant () const
  {
    return 1.0;
  }

  // build jacboianInverseTransposed for (mydim < cdim)
  template< int mydim, int cdim, class GridImp >
  inline void
  AlbertaGridGeometry< mydim, cdim, GridImp >
  :: buildJacobianInverseTransposed () const
  {
    // calc A and stores it in elMat_
    calcElMatrix();
    assert( builtElMat_ == true );

    // Jinv = A^-1^T
    elDet_ = Alberta::MatrixOps< ctype, cdim, mydim >::invert( elMat_, Jinv_ );
    assert( elDet_ > 1.0e-25 );

    calcedDet_ = true;
    builtinverse_ = true;
  }


  template <int mydim, int cdim, class GridImp>
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::ctype
  AlbertaGridGeometry<mydim,cdim,GridImp>::volume () const
  {
    assert( calcedDet_ );
    const ctype refVolume
      = ctype( 1 ) / ctype( Factorial< mydimension > :: factorial );
    return refVolume * elDet_;
  }


  template <int mydim, int cdim, class GridImp >
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::ctype
  AlbertaGridGeometry< mydim, cdim, GridImp>
  ::integrationElement ( const LocalVector &local ) const
  {
    assert( calcedDet_ );
    return elDet_;
  }


  template< int mydim, int cdim, class GridImp >
  inline const FieldMatrix< Alberta::Real, cdim, mydim > &
  AlbertaGridGeometry< mydim, cdim, GridImp >
  ::jacobianInverseTransposed ( const LocalVector &local ) const
  {
    if(builtinverse_)
      return Jinv_;

    // builds the jacobian inverse and calculates the volume
    buildJacobianInverseTransposed();
    return Jinv_;
  }


  template< int mydim, int cdim, class GridImp >
  inline bool AlbertaGridGeometry< mydim, cdim, GridImp >
  ::checkInside ( const LocalVector &local ) const
  {
    albertCtype sum = 0.0;

    for(int i=0; i<mydim; i++)
    {
      sum += local[i];
      if(local[i] < 0.0)
      {
        if(std::abs(local[i]) > 1e-13)
        {
          return false;
        }
      }
    }

    if( sum > 1.0 )
    {
      if(sum > (1.0 + 1e-13))
        return false;
    }

    return true;
  }


  template< int mydim, int cdim, class GridImp >
  inline void AlbertaGridGeometry< mydim, cdim, GridImp >::invalidate ()
  {
    builtinverse_ = false;
    builtElMat_ = false;
    calcedDet_ = false;
  }


  // built Geometry
  template< int mydim, int cdim, class GridImp >
  template< class CoordReader >
  inline void AlbertaGridGeometry< mydim, cdim, GridImp >
  ::build ( const CoordReader &coordReader )
  {
    builtinverse_ = false;
    builtElMat_ = false;

    for( int i = 0; i <= mydimension; ++i )
      coordReader.coordinate( i, coord_[ i ] );

    elDet_ = (coordReader.hasDeterminant() ? coordReader.determinant() : elDeterminant());
    assert( std::abs( elDet_ ) > 0.0 );
    calcedDet_ = true;
  }



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

    FieldMatrix< ctype, dimension+1, dimension > refCorners = 0.0;
    for( int i = 0; i < dimension; ++i )
      refCorners[ i+1 ][ i ] = 1.0;

    for( int face = 0; face < numFaces; ++face )
    {
      // use reference element here!
      FieldMatrix< ctype, dimension, dimension > faceCorners;
      for( int i = 0; i < face; ++i )
        faceCorners[ i ] = refCorners[ i ];
      for( int i = face; i < dimension; ++i )
        faceCorners[ i ] = refCorners[ i+1 ];
      faceGeometry_[ face ] = new LocalGeoObject( LocalGeoImp( faceCorners ) );
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


}

#endif
