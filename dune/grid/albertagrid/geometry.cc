// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::GlobalCoordinate
  AlbertaGridGeometry< mydim, cdim, GridImp >::global ( const LocalCoordinate &local ) const
  {
    GlobalCoordinate y = corner( 0 );
    jacobianTransposed().umtv( local, y );
    return y;
  }


  //local implementation for mydim < cdim
  template< int mydim, int cdim, class GridImp >
  inline typename AlbertaGridGeometry< mydim, cdim, GridImp >::LocalCoordinate
  AlbertaGridGeometry< mydim, cdim, GridImp >::local ( const GlobalCoordinate &global ) const
  {
    LocalCoordinate x;
    jacobianInverseTransposed().mtv( global - corner( 0 ), x );
    return x;
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


  // built Geometry
  template< int mydim, int cdim, class GridImp >
  template< class CoordReader >
  inline void AlbertaGridGeometry< mydim, cdim, GridImp >
  ::build ( const CoordReader &coordReader )
  {
    builtJT_ = false;
    builtJTInv_ = false;

    // copy corners
    for( int i = 0; i < numCorners; ++i )
      coordReader.coordinate( i, coord_[ i ] );

    // calculate centroid
    centroid_ = coord_[ 0 ];
    for( int i = 1; i < numCorners; ++i )
      centroid_ += coord_[ i ];
    centroid_ *= 1.0 / numCorners;

    elDet_ = (coordReader.hasDeterminant() ? coordReader.determinant() : elDeterminant());
    assert( std::abs( elDet_ ) > 0.0 );
    calcedDet_ = true;
  }


#if !DUNE_ALBERTA_CACHE_COORDINATES
  template< int dim, int cdim >
  inline typename AlbertaGridGlobalGeometry< dim, cdim, const AlbertaGrid< dim, cdim > >::GlobalCoordinate
  AlbertaGridGlobalGeometry< dim, cdim, const AlbertaGrid< dim, cdim > >::global ( const LocalCoordinate &local ) const
  {
    GlobalCoordinate y = corner( 0 );
    jacobianTransposed().umtv( local, y );
    return y;
  }


  //local implementation for mydim < cdim
  template< int dim, int cdim >
  inline typename AlbertaGridGlobalGeometry< dim, cdim, const AlbertaGrid< dim, cdim > >::LocalCoordinate
  AlbertaGridGlobalGeometry< dim, cdim, const AlbertaGrid< dim, cdim > >::local ( const GlobalCoordinate &global ) const
  {
    LocalCoordinate x;
    jacobianInverseTransposed().mtv( global - corner( 0 ), x );
    return x;
  }
#endif // #if !DUNE_ALBERTA_CACHE_COORDINATES



  // AlbertaGridLocalGeometryProvider
  // --------------------------------

  template< class Grid >
  void AlbertaGridLocalGeometryProvider< Grid >::buildGeometryInFather ()
  {
    for( int child = 0; child < numChildren; ++child )
    {
      for( int orientation = 0; orientation < 2; ++orientation )
      {
        const GeoInFatherCoordReader coordReader( child, orientation );
        geometryInFather_[ child ][ orientation ] = new LocalElementGeometry( coordReader );
      }
    }
  }


  template< class Grid >
  void AlbertaGridLocalGeometryProvider< Grid >::buildFaceGeometry ()
  {
    for( int face = 0; face < numFaces; ++face )
    {
      for( int twist = minFaceTwist; twist <= maxFaceTwist; ++twist )
      {
        const FaceCoordReader coordReader( face, twist );
        faceGeometry_[ face ][ twist - minFaceTwist ] = new LocalFaceGeometry( coordReader );
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

  private:
    const int child_;
    const int orientation_;
  };



  // AlbertaGridLocalGeometryProvider::FaceCoordReader
  // -------------------------------------------------

  template< class Grid >
  struct AlbertaGridLocalGeometryProvider< Grid >::FaceCoordReader
  {
    typedef Alberta::Real ctype;

    typedef FieldVector< ctype, dimension > Coordinate;

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

    int face_;
    int twist_;
  };

} // namespace Dune

#endif // #ifndef DUNE_ALBERTA_GEOMETRY_CC
