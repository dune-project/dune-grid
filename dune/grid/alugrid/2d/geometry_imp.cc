// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_GEOMETRYIMP_CC
#define DUNE_ALU2DGRID_GEOMETRYIMP_CC

#include <dune/geometry/referenceelements.hh>

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
  {
    getObject();
  }

  template< int mydim, int cdim, class GridImp >
  inline ALU2dGridGeometry< mydim, cdim, GridImp >&
  ALU2dGridGeometry< mydim, cdim, GridImp >::operator = ( const ALU2dGridGeometry& other )
  {
    if( &other != this )
    {
      removeObj();
      assign( other );
    }
    return *this;
  }

  template< int mydim, int cdim, class GridImp >
  inline ALU2dGridGeometry< mydim, cdim, GridImp >::ALU2dGridGeometry ( const ALU2dGridGeometry& other )
  {
    assign( other );
  }

  template< int mydim, int cdim, class GridImp >
  inline ALU2dGridGeometry< mydim, cdim, GridImp >::~ALU2dGridGeometry ()
  {
    removeObj();
  }

  template< int mydim, int cdim, class GridImp >
  inline void ALU2dGridGeometry< mydim, cdim, GridImp >::getObject()
  {
    geoImpl_ = geoProvider().getEmptyObject();
    geoImpl().reset();
  }

  template< int mydim, int cdim, class GridImp >
  inline void ALU2dGridGeometry< mydim, cdim, GridImp >::assign( const ALU2dGridGeometry& other )
  {
    // copy pointer
    geoImpl_ = other.geoImpl_ ;

    // increase reference count
    ++ geoImpl();
  }

  template< int mydim, int cdim, class GridImp >
  inline void ALU2dGridGeometry< mydim, cdim, GridImp >::removeObj()
  {
    // decrease reference count
    -- geoImpl();

    // if reference count is zero free the object
    if( ! geoImpl() )
    {
      geoProvider().freeObject( geoImpl_ );
    }

    // reset pointer
    geoImpl_ = 0;
  }

  template< int mydim, int cdim, class GridImp >
  inline void ALU2dGridGeometry< mydim, cdim, GridImp >::invalidate()
  {
    // if geometry is used elsewhere remove the pointer
    // and get a new one
    if( geoImpl().stillUsed() )
    {
      // remove old object
      removeObj();

      // get new object
      getObject();
    }
    else
    {
      // otherwise invalidate object
      geoImpl().invalidate();
    }
  }

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
  template< int mydim, int cdim, class GridImp >
  inline typename ALU2dGridGeometry< mydim, cdim, GridImp >::GlobalCoordinate
  ALU2dGridGeometry< mydim, cdim, GridImp >::corner ( int i ) const
  {
    return geoImpl().corner( i );
  }


  //! maps a local coordinate within reference element to
  //! global coordinate in element
  template< int mydim, int cdim, class GridImp >
  inline typename ALU2dGridGeometry< mydim, cdim, GridImp >::GlobalCoordinate
  ALU2dGridGeometry< mydim, cdim, GridImp >::global ( const LocalCoordinate &local ) const
  {
    GlobalCoordinate global;
    geoImpl().map2world( local, global );
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
    geoImpl().world2map( global, local );
    return local;
  }


  template< int mydim, int cdim, class GridImp >
  inline alu2d_ctype
  ALU2dGridGeometry< mydim, cdim, GridImp >::integrationElement ( const LocalCoordinate &local ) const
  {
    assert( geoImpl().valid() );
    if( mydim == cdim )
    {
      if( std::abs( geoImpl().det( local ) - jacobianTransposed( local
                                                                 ).determinant() ) > 1e-8 )
        std::cout << "det = " << geoImpl().det( local ) << "  inv.det " <<
        jacobianTransposed( local ).determinant() << std::endl;
      assert( std::abs( geoImpl().det( local ) - jacobianTransposed( local
                                                                     ).determinant() ) < 1e-8 );
    }
    return geoImpl().det( local );
  }


  template< int mydim, int cdim, class GridImp >
  inline alu2d_ctype ALU2dGridGeometry< mydim, cdim, GridImp >::volume () const
  {
    assert( geoImpl().valid() );
    return geoImpl().volume();
  }


  template< int mydim, int cdim, class GridImp >
  inline const FieldMatrix<alu2d_ctype,mydim,cdim>&
  ALU2dGridGeometry< mydim, cdim, GridImp>::jacobianTransposed ( const LocalCoordinate &local ) const
  {
    return geoImpl().jacobianTransposed( local );
  }


  template< int mydim, int cdim, class GridImp >
  inline const FieldMatrix<alu2d_ctype,cdim,mydim>&
  ALU2dGridGeometry< mydim, cdim, GridImp >::jacobianInverseTransposed ( const LocalCoordinate &local ) const
  {
    return geoImpl().jacobianInverseTransposed( local );
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
    geoImpl().update( item );

    // geometry built
    return true;
  }


  //! built Geometry for edges
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  buildGeom(const HElementType & item, const int aluFace)
  {
    assert( mydim == 1 );
    // face 1 is twisted
    const int nf = item.numfaces();

    // for triangles face 1 is twisted and for quatrilaterals faces 1 and 2
    const int twist = (nf == 3) ? (aluFace % 2) : (aluFace>>1)^(aluFace&1);
    // std::cout << nf << " " << aluFace << " " << twist << std::endl;
    //           << " ( " item.getVertex( (aluFace + 1 + twist ) % nf )->coord() ,
    //           << " , " item.getVertex( (aluFace + 2 - twist ) % nf )->coord()
    //           << " ) " << std::endl;

    // check item
    assert( &item );


    // update geometry impl
    geoImpl().update( item.getVertex( (aluFace + 1 + twist ) % nf )->coord() ,
                      item.getVertex( (aluFace + 2 - twist ) % nf )->coord() ,
                      item.sidelength( aluFace )
                      );

    // geometry built
    return true;
  }

  //! built Geometry for vertices
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  buildGeom(const VertexType & item , const int )
  {
    assert( mydim == 0 );

    assert( &item );
    // update geometry impl
    // volume is already 1.0
    geoImpl().update( item.coord() );

    // geometry built
    return true;
  }

  // built Geometry
  template <int mydim, int cdim, class GridImp>
  template <class GeometryType, class LocalGeometryType >
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  buildLocalGeom(const GeometryType &geo, const LocalGeometryType & localGeom)
  {
    // update geometry
    geoImpl().updateLocal( geo, localGeom );
    // geometry built
    return true;
  }

  // built Geometry (faceNumber is in generic numbering)
  template< int mydim, int cdim, class GridImp >
  inline std::pair< FieldMatrix< alu2d_ctype, 4, 2 >, FieldVector< alu2d_ctype, 4 > >
  ALU2dGridGeometry< mydim, cdim, GridImp >::calculateReferenceCoords ( const int corners )
  {
    // calculate reference coordinates of aLUGrid reference triangle

    FieldMatrix< alu2d_ctype, 4, 2 > refCoord( 0. );
    FieldVector< alu2d_ctype, 4 > lengths( 1. );

    // point 1
    refCoord[1][0] = 1.0;
    // point (corners-1)
    refCoord[corners-1][1] = 1.0;
    if( corners == 3 )
    {
      lengths[0] = M_SQRT2;
    }
    else
    {
      // point 2
      refCoord[2][0] = 1.0;
      refCoord[2][1] = 1.0;
    }
    return std::make_pair( refCoord, lengths );
  }

  // built Geometry (faceNumber is in generic numbering)
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  buildLocalGeometry(const int aluFace, const int twist, const int corners)
  {
    assert( twist == 0 || twist == 1 );
    assert( mydim == 1 );

    // get coordinates of reference element
    typedef std::pair< FieldMatrix< alu2d_ctype, 4, 2 >, FieldVector< alu2d_ctype, 4 > > RefCoord;
    RefCoord refCoord( calculateReferenceCoords( corners ) );

    geoImpl().update( refCoord.first[ ( aluFace + 1+twist ) % corners ],
                      refCoord.first[ ( aluFace + 2-twist ) % corners ],
                      refCoord.second[ aluFace ]
                      );

    // geometry built
    return true;
  }

  // built Geometry
  template <int mydim, int cdim, class GridImp >
  inline bool ALU2dGridGeometry<mydim, cdim,GridImp>::
  buildGeomInFather(const Geometry & fatherGeom ,
                    const Geometry & myGeom)
  {
    // update geometry
    geoImpl().updateLocal( fatherGeom, myGeom );
    // geometry built
    return true;
  }

} //end namespace Dune

#endif // #ifndef DUNE_ALU2DGRID_GEOMETRYIMP_CC
