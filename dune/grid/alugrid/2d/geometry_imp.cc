// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_GEOMETRYIMP_CC
#define DUNE_ALU2DGRID_GEOMETRYIMP_CC

#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/common/genericreferenceelements.hh>

namespace Dune
{

  // implementation of methods

  //**********************************************************************
  //
  // --ALU2dGridGeometry
  // --Geometry
  //**********************************************************************


  template< int mydim, int cdim, class GridImp >
  inline ALU2dGridGeometry< mydim, cdim, GridImp >::ALU2dGridGeometry ()
    : geoImpl_(),
      det_( 1.0 ),
      up2Date_( false )
#ifndef NDEBUG
      , haveProjection_( false )
#endif
  {}


  //! print the GeometryInformation
  template< int mydim, int cdim, class GridImp >
  inline void ALU2dGridGeometry< mydim, cdim, GridImp >::print ( std::ostream &out ) const
  {
    out << "ALU2dGridGeometry< " << mydim << ", " << cdim << " > = ";
    char c = '{';
    for( int i = 0; i < corners(); ++i )
    {
      out << c << " (" << corner( i ) << ")";
      c = ',';
    }
    out << " }" << std::endl;
  }


  ///////////////////////////////////////////////////////////////////////


  //! access to coordinates of corners. Index is the number of the corner
  template <int mydim, int cdim, class GridImp>
  inline const typename ALU2dGridGeometry< mydim, cdim, GridImp >::GlobalCoordinate &
  ALU2dGridGeometry< mydim, cdim, GridImp >::operator[] ( int i ) const
  {
    typedef GenericGeometry::MapNumberingProvider< mydim > Numbering;
    const unsigned int tid = GenericGeometry::topologyId( type() );
    static GlobalCoordinate c;
    return (c = corner( Numbering::template dune2generic< mydim >( tid, i ) ));
  }


  //! access to coordinates of corners. Index is the number of the corner
  template< int mydim, int cdim, class GridImp >
  inline typename ALU2dGridGeometry< mydim, cdim, GridImp >::GlobalCoordinate
  ALU2dGridGeometry< mydim, cdim, GridImp >::corner ( int i ) const
  {
    const GenericReferenceElement< alu2d_ctype, mydim > &refElement
      = GenericReferenceElements< alu2d_ctype, mydim >::general( type() );
    return global( refElement.position( i, mydim ) );
  }


  //! maps a local coordinate within reference element to
  //! global coordinate in element
  template< int mydim, int cdim, class GridImp >
  inline typename ALU2dGridGeometry< mydim, cdim, GridImp >::GlobalCoordinate
  ALU2dGridGeometry< mydim, cdim, GridImp >::global ( const LocalCoordinate &local ) const
  {
    GlobalCoordinate global;
    geoImpl_.map2world( local, global );
    return global;
  }


  //! maps a global coordinate within the element to a
  //! local coordinate in its reference element
  template< int mydim, int cdim, class GridImp >
  inline typename ALU2dGridGeometry< mydim, cdim, GridImp >::LocalCoordinate
  ALU2dGridGeometry< mydim, cdim, GridImp >::local ( const GlobalCoordinate &global ) const
  {
    if( mydim == 0 )
      return LocalCoordinate( 1 );

    LocalCoordinate local;
    geoImpl_.world2map( global, local );
    return local;
  }


  template< int mydim, int cdim, class GridImp >
  inline alu2d_ctype
  ALU2dGridGeometry< mydim, cdim, GridImp >::integrationElement ( const LocalCoordinate &local ) const
  {
    return (mydim == 0 ? 1.0 : det_);
  }


  template< int mydim, int cdim, class GridImp >
  inline alu2d_ctype ALU2dGridGeometry< mydim, cdim, GridImp >::volume () const
  {
    if( mydim == 2 )
    {
      switch( GridImp::elementType )
      {
      case ALU2DSPACE triangle :
        return 0.5 * det_;

      case ALU2DSPACE quadrilateral :
        return det_;

      case ALU2DSPACE mixed :
        DUNE_THROW( NotImplemented, "Geometry::volume() not implemented for ElementType mixed." );
      }
    }
    else
      return (mydim == 0 ? 1.0 : det_);
  }


  template< int mydim, int cdim, class GridImp >
  inline const FieldMatrix<alu2d_ctype,mydim,cdim>&
  ALU2dGridGeometry< mydim, cdim, GridImp>::jacobianTransposed ( const LocalCoordinate &local ) const
  {
    return geoImpl_.jacobianTransposed( local );
  }


  template< int mydim, int cdim, class GridImp >
  inline const FieldMatrix<alu2d_ctype,cdim,mydim>&
  ALU2dGridGeometry< mydim, cdim, GridImp >::jacobianInverseTransposed ( const LocalCoordinate &local ) const
  {
    return geoImpl_.jacobianInverseTransposed( local );
  }


  //! built Geometry for triangles
  template< int mydim, int cdim, class GridImp >
  inline bool
  ALU2dGridGeometry< mydim, cdim, GridImp >::buildGeom( const HElementType &item )
  {
    // check item
    assert( &item );
    assert( mydim == 2 );

    // update geometry impl
    geoImpl_.update( item );

    // store volume
    det_ = (geoImpl_.corners() == 3 ? 2 * item.area() : item.area());

    // geom is up2date
    up2Date_ = true;

    // geometry built
    return true;
  }


  //! built Geometry for edges
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  buildGeom(const HElementType & item, const int aluFace)
  {
    // face 1 is twisted
    const int twist = (aluFace % 2);

    // check item
    assert( &item );

    // update geometry impl
    geoImpl_.update( item.getVertex( (aluFace + (1 + twist)) % 3 )->coord() ,
                     item.getVertex( (aluFace + (2 - twist)) % 3 )->coord() );

    // store volume
    det_ = item.sidelength( aluFace );

    // geom is up2date
    up2Date_ = true;

    // geometry built
    return true;
  }

  //! built Geometry for vertices
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  buildGeom(const VertexType & item , const int )
  {
    assert( &item );
    // update geometry impl
    geoImpl_.update( item.coord() );

    // volume is already 1.0

    // geom is up2date
    up2Date_ = true;

    return true;
  }

  // built Geometry
  template <int mydim, int cdim, class GridImp>
  template <class GeometryType, class LocalGeometryType >
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  buildLocalGeom(const GeometryType &geo, const LocalGeometryType & localGeom)
  {
    // update geometry
    geoImpl_.updateLocal( geo, localGeom );

    // calculate volume
    LocalCoordinate local( 0.25 );
    det_ = geoImpl_.det( local );

    // geom is up2date
    up2Date_ = true;

    // geometry built
    return true;
  }

  // built Geometry (faceNumber is in generic numbering)
  template <int mydim, int cdim, class GridImp>
  inline FieldMatrix<alu2d_ctype, 3 , 3> ALU2dGridGeometry<mydim,cdim,GridImp>::
  calculateReferenceCoords() const
  {
    // calculate reference coordinates of aLUGrid reference triangle

    // set all to zero
    FieldMatrix<alu2d_ctype, 3 , 3> refCoord (0.0);

    // point 1
    refCoord[1][0] = 1.0;
    // point 2
    refCoord[2][1] = 1.0;

    // length of faces
    refCoord[0][2] = M_SQRT2;
    refCoord[1][2] = 1.0;
    refCoord[2][2] = 1.0;

    return refCoord;
  }

  // built Geometry (faceNumber is in generic numbering)
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  buildLocalGeometry(const int aluFace, const int twist)
  {
    assert( twist == 0 || twist == 1 );
    assert( mydim == 1 );

    // get coordinates of reference element
    static const FieldMatrix<alu2d_ctype, 3, 3> refCoord ( calculateReferenceCoords() );

    // just map the point of the global intersection to the local
    // coordinates , this is the default procedure
    // for simplices this is not so bad
    geoImpl_.update( refCoord[( aluFace + ( twist   %2) + 1) % 3],
                     refCoord[( aluFace + ((twist+1)%2) + 1) % 3] );

    // get length of faces
    det_ = refCoord[ aluFace ][2];

    // geom is up2date
    up2Date_ = true;

    // geometry built
    return true;
  }

  // built Geometry
  template <int mydim, int cdim, class GridImp >
  inline bool ALU2dGridGeometry<mydim, cdim,GridImp>::
  buildGeomInFather(const Geometry & fatherGeom ,
                    const Geometry & myGeom,
                    const bool hasBndProjection )
  {
    // update geometry
    geoImpl_.updateLocal( fatherGeom, myGeom );

    // store volume which is a part of one
    det_ = myGeom.volume() / fatherGeom.volume();
    assert( (det_ > 0.0) && (det_ < 1.0) );

    // geom is up2date
    up2Date_ = true;

#ifndef NDEBUG
    // set projection information
    haveProjection_ = hasBndProjection;
#endif

    return true;
  }

} //end namespace Dune

#endif // #ifndef DUNE_ALU2DGRID_GEOMETRYIMP_CC
