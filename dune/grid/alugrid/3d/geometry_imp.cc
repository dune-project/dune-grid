// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_GEOMETRY_IMP_CC
#define DUNE_ALUGRID_GEOMETRY_IMP_CC

#include <dune/common/math.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

#include "grid.hh"
#include "mappings.hh"
#include "geometry.hh"

namespace Dune {
  // --Geometry

  template <int mydim, int cdim, class GridImp>
  inline ALU3dGridGeometry<mydim, cdim, GridImp> ::
  ALU3dGridGeometry()
  {
    getObject();
  }

  template <int mydim, int cdim, class GridImp>
  inline ALU3dGridGeometry<mydim, cdim, GridImp> ::
  ALU3dGridGeometry( const ALU3dGridGeometry& other )
  {
    assign( other );
  }

  template <int mydim, int cdim, class GridImp>
  inline void ALU3dGridGeometry<mydim, cdim, GridImp> ::
  getObject ( )
  {
    geoImpl_ = geoProvider().getEmptyObject();
    geoImpl().reset();
  }

  template <int mydim, int cdim, class GridImp>
  inline void ALU3dGridGeometry<mydim, cdim, GridImp> ::
  assign ( const ALU3dGridGeometry& other )
  {
    // copy pointer
    geoImpl_ = other.geoImpl_ ;

    // increase reference count
    ++ geoImpl();
  }

  template <int mydim, int cdim, class GridImp>
  inline void ALU3dGridGeometry<mydim, cdim, GridImp> ::
  removeObj ()
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

  template <int mydim, int cdim, class GridImp>
  inline ALU3dGridGeometry<mydim, cdim, GridImp>&
  ALU3dGridGeometry<mydim, cdim, GridImp> :: operator = (const ALU3dGridGeometry& other )
  {
    if( &other != this )
    {
      removeObj();
      assign( other );
    }
    return *this;
  }

  template <int mydim, int cdim, class GridImp>
  inline ALU3dGridGeometry<mydim, cdim, GridImp> ::
  ~ALU3dGridGeometry()
  {
    removeObj();
  }

  template< int mydim, int cdim, class GridImp>
  inline void
  ALU3dGridGeometry< mydim, cdim, GridImp > :: invalidate ()
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

  template< int mydim, int cdim, class GridImp>
  inline bool
  ALU3dGridGeometry< mydim, cdim, GridImp > :: valid () const
  {
    return geoImpl().valid();
  }

  template< int mydim, int cdim, class GridImp>
  inline GeometryType
  ALU3dGridGeometry< mydim, cdim, GridImp > :: type () const
  {
    return GeometryType( (elementType == tetra) ?
                         GenericGeometry :: SimplexTopology< mydim > :: type :: id :
                         GenericGeometry :: CubeTopology   < mydim > :: type :: id,
                         mydim );
  }

  template< int mydim, int cdim, class GridImp>
  inline int
  ALU3dGridGeometry<mydim, cdim, GridImp >::corners() const
  {
    return corners_;
  }

  template< int mydim, int cdim, class GridImp>
  inline typename ALU3dGridGeometry<mydim, cdim, GridImp >::GlobalCoordinate
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  corner (int i) const
  {
    return geoImpl()[ i ];
  }


  template< int mydim, int cdim, class GridImp>
  inline typename ALU3dGridGeometry<mydim, cdim, GridImp >::GlobalCoordinate
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  global (const LocalCoordinate& local) const
  {
    GlobalCoordinate global;
    geoImpl().mapping().map2world(local, global);
    return global;
  }

  template< int mydim, int cdim, class GridImp >
  inline typename ALU3dGridGeometry<mydim, cdim, GridImp >::LocalCoordinate
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  local (const GlobalCoordinate& global) const
  {
    LocalCoordinate local;
    geoImpl().mapping().world2map(global, local);
    return local;
  }

  template< int mydim, int cdim, class GridImp>
  inline typename ALU3dGridGeometry<mydim, cdim, GridImp >::ctype
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  integrationElement (const LocalCoordinate& local) const
  {
    // this is the only case we need to specialize
    if( mydim == cdim && elementType == tetra )
    {
      assert( geoImpl().valid() );
      return 6.0 * geoImpl().volume();
    }
    else
      return geoImpl().mapping().det( local );
  }

  template<int mydim, int cdim, class GridImp>
  inline typename ALU3dGridGeometry<mydim, cdim, GridImp >::ctype
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  volume () const
  {
    if( mydim == cdim )
    {
      assert( geoImpl().valid() );
      return geoImpl().volume() ;
    }
    else if ( mydim == cdim - 1 && elementType == tetra )
    {
      enum { factor = Factorial<mydim>::factorial };
      // local vector does not affect the result
      const LocalCoordinate dummy(0);
      return integrationElement( dummy ) / ((ctype) factor);
    }
    else
    {
      return integrationElement(LocalCoordinate(0.5));
    }
  }

  template< int mydim, int cdim, class GridImp>
  inline bool
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  affine() const
  {
    return geoImpl().mapping().affine();
  }

  template< int mydim, int cdim, class GridImp>
  inline const typename ALU3dGridGeometry<mydim, cdim, GridImp >::JacobianInverseTransposed&
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  jacobianInverseTransposed (const LocalCoordinate & local) const
  {
    return geoImpl().mapping().jacobianInverseTransposed( local );
  }

  template< int mydim, int cdim, class GridImp>
  inline const typename ALU3dGridGeometry<mydim, cdim, GridImp >::JacobianTransposed&
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  jacobianTransposed (const LocalCoordinate & local) const
  {
    return geoImpl().mapping().jacobianTransposed( local );
  }

  template <int mydim, int cdim, class GridImp>
  inline void
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  print (std::ostream& ss) const
  {
    const char* charElType = (elementType == tetra) ? "tetra" : "hexa";
    ss << "ALU3dGridGeometry<" << mydim << "," << cdim << ", " << charElType << "> = {\n";
    for(int i=0; i<corners(); i++)
    {
      ss << " corner " << i << " ";
      ss << "{" << ((*this)[i]) << "}"; ss << std::endl;
    }
    ss << "} \n";
  }

  // built Geometry
  template <int mydim, int cdim, class GridImp>
  template <class GeometryType>
  inline bool
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  buildGeomInFather(const GeometryType &fatherGeom , const GeometryType & myGeom)
  {
    // update geo impl
    geoImpl().updateInFather( fatherGeom, myGeom );

    // my volume is a part of 1 for hexas, for tetra adjust with factor
    double volume = myGeom.volume() / fatherGeom.volume() ;
    if( elementType == tetra )
    {
      volume /= 6.0;
      geoImpl().setVolume( volume );
#ifndef NDEBUG
      LocalCoordinate local( 0.0 );
      assert( std::abs( 6.0 * geoImpl().volume() - integrationElement( local ) ) < 1e-12 );
#endif
    }
    else
      geoImpl().setVolume( volume );

    return true;
  }

  //--hexaBuildGeom
  template <int mydim, int cdim, class GridImp>
  inline bool
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  buildGeom(const IMPLElementType& item)
  {
    if ( elementType == hexa )
    {
      // if this assertion is thrown, use ElementTopo::dune2aluVertex instead
      // of number when calling myvertex
      assert( ElementTopo::dune2aluVertex(0) == 0 );
      assert( ElementTopo::dune2aluVertex(1) == 1 );
      assert( ElementTopo::dune2aluVertex(2) == 3 );
      assert( ElementTopo::dune2aluVertex(3) == 2 );
      assert( ElementTopo::dune2aluVertex(4) == 4 );
      assert( ElementTopo::dune2aluVertex(5) == 5 );
      assert( ElementTopo::dune2aluVertex(6) == 7 );
      assert( ElementTopo::dune2aluVertex(7) == 6 );

      // update geo impl
      geoImpl().update( item.myvertex(0)->Point(),
                        item.myvertex(1)->Point(),
                        item.myvertex(3)->Point(),
                        item.myvertex(2)->Point(),
                        item.myvertex(4)->Point(),
                        item.myvertex(5)->Point(),
                        item.myvertex(7)->Point(),
                        item.myvertex(6)->Point() );
    }
    else if( elementType == tetra )
    {
      // if this assertion is thrown, use ElementTopo::dune2aluVertex instead
      // of number when calling myvertex
      assert( ElementTopo::dune2aluVertex(0) == 0 );
      assert( ElementTopo::dune2aluVertex(1) == 1 );
      assert( ElementTopo::dune2aluVertex(2) == 2 );
      assert( ElementTopo::dune2aluVertex(3) == 3 );

      // update geo impl
      geoImpl().update( item.myvertex(0)->Point(),
                        item.myvertex(1)->Point(),
                        item.myvertex(2)->Point(),
                        item.myvertex(3)->Point() );
    }

    // get volume of element
    geoImpl().setVolume( item.volume() );

    return true;
  }

  // buildFaceGeom
  template <int mydim, int cdim, class GridImp>
  inline bool
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  buildGeom(const HFaceType & item, int twist, int duneFace )
  {
    // get geo face
    const GEOFaceType& face = static_cast<const GEOFaceType&> (item);

    // if face was not set (when face comes from face iteration),
    // then set it to zero
    if( duneFace < 0 ) duneFace = 0;

    enum { numVertices = ElementTopo::numVerticesPerFace };
    // for all vertices of this face get rotatedIndex
    int rotatedALUIndex[ numVertices ];
    for (int i = 0; i < numVertices; ++i)
    {
      // Transform Dune index to ALU index and apply twist
      const int localALUIndex = ElementTopo::dune2aluFaceVertex(duneFace,i);
      rotatedALUIndex[ i ] = FaceTopo::twist(localALUIndex, twist);
    }

    if( elementType == hexa )
    {
      // update geometry implementation
      geoImpl().update( face.myvertex(rotatedALUIndex[0])->Point(),
                        face.myvertex(rotatedALUIndex[1])->Point(),
                        face.myvertex(rotatedALUIndex[2])->Point(),
                        face.myvertex(rotatedALUIndex[3])->Point() );
    }
    else if ( elementType == tetra )
    {
      // update geometry implementation
      geoImpl().update( face.myvertex(rotatedALUIndex[0])->Point(),
                        face.myvertex(rotatedALUIndex[1])->Point(),
                        face.myvertex(rotatedALUIndex[2])->Point());
    }

    return true;
  }

  // --buildFaceGeom
  template <int mydim, int cdim, class GridImp>
  template <class coord_t>
  inline bool
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  buildGeom(const coord_t& p0,
            const coord_t& p1,
            const coord_t& p2,
            const coord_t& p3)
  {
    // update geometry implementation
    geoImpl().update( p0, p1, p2, p3 );
    return true;
  }

  // --buildFaceGeom
  template <int mydim, int cdim, class GridImp>
  template <class coord_t>
  inline bool
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  buildGeom(const coord_t& p0,
            const coord_t& p1,
            const coord_t& p2)
  {
    // update geometry implementation
    geoImpl().update( p0, p1, p2 );
    return true;
  }

  template <int mydim, int cdim, class GridImp> // for faces
  inline bool
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  buildGeom(const FaceCoordinatesType& coords)
  {
    if ( elementType == hexa )
      return buildGeom( coords[0], coords[1], coords[2], coords[3] );
    else
    {
      assert( elementType == tetra );
      return buildGeom( coords[0], coords[1], coords[2] );
    }
  }

  template <int mydim, int cdim, class GridImp> // for edges
  inline bool
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  buildGeom(const HEdgeType & item, int twist, int)
  {
    const GEOEdgeType & edge = static_cast<const GEOEdgeType &> (item);
    // update geometry implementation
    geoImpl().update( edge.myvertex((twist)  %2)->Point(),
                      edge.myvertex((1+twist)%2)->Point() );
    return true;
  }

  template <int mydim, int cdim, class GridImp> // for Vertices ,i.e. Points
  inline bool
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  buildGeom(const VertexType & item, int twist, int)
  {
    // update geometry implementation
    geoImpl().update( static_cast<const GEOVertexType &> (item).Point() );
    return true;
  }

} // end namespace Dune
#endif // end DUNE_ALUGRID_GEOMETRY_IMP_CC
