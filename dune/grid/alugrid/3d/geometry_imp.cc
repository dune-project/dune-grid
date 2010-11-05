// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_GEOMETRY_IMP_CC
#define DUNE_ALUGRID_GEOMETRY_IMP_CC

#include <dune/grid/genericgeometry/conversion.hh>

#include "grid.hh"
#include "mappings.hh"
#include "geometry.hh"

namespace Dune {
  // --Geometry

  template <int mydim, int cdim, class GridImp>
  inline ALU3dGridGeometry<mydim, cdim, GridImp> ::
  ALU3dGridGeometry()
    : geoImpl_(),
      volume_(1.0)
  {}

  template< int mydim, int cdim, class GridImp>
  inline GeometryType
  ALU3dGridGeometry< mydim, cdim, GridImp > :: type () const
  {
    return GeometryType( (elementType == tetra) ?
                         GeometryType :: simplex : GeometryType :: cube, mydim );
  }

  template< int mydim, int cdim, class GridImp>
  inline int
  ALU3dGridGeometry<mydim, cdim, GridImp >::corners() const {
    return corners_;
  }

  template< int mydim, int cdim, class GridImp>
  inline const typename ALU3dGridGeometry<mydim, cdim, GridImp >::GlobalCoordinate&
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  operator[] (int i) const
  {
    typedef GenericGeometry::MapNumberingProvider< mydim > Numbering;
    const unsigned int tid = GenericGeometry::topologyId( type() );
    const int j = Numbering::template dune2generic< mydim >( tid, i );
    return geoImpl_[ j ];
  }

  template< int mydim, int cdim, class GridImp>
  inline typename ALU3dGridGeometry<mydim, cdim, GridImp >::GlobalCoordinate
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  corner (int i) const
  {
    return geoImpl_[ i ];
  }


  template< int mydim, int cdim, class GridImp>
  inline typename ALU3dGridGeometry<mydim, cdim, GridImp >::GlobalCoordinate
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  global (const LocalCoordinate& local) const
  {
    GlobalCoordinate global;
    geoImpl_.mapping().map2world(local, global);
    return global;
  }

  template< int mydim, int cdim, class GridImp >
  inline typename ALU3dGridGeometry<mydim, cdim, GridImp >::LocalCoordinate
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  local (const GlobalCoordinate& global) const
  {
    LocalCoordinate local;
    geoImpl_.mapping().world2map(global, local);
    return local;
  }

  template< int mydim, int cdim, class GridImp>
  inline typename ALU3dGridGeometry<mydim, cdim, GridImp >::ctype
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  integrationElement (const LocalCoordinate& local) const
  {
    if( mydim == cdim && elementType == tetra )
      return 6.0 * volume_;
    else if ( mydim == 0 )
      return 1.0;
    else
      return geoImpl_.mapping().det( local );
  }

  template<int mydim, int cdim, class GridImp>
  inline typename ALU3dGridGeometry<mydim, cdim, GridImp >::ctype
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  volume () const
  {
    if( mydim == cdim )
    {
      return volume_ ;
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
    return geoImpl_.mapping().affine();
  }

  template< int mydim, int cdim, class GridImp>
  inline const typename ALU3dGridGeometry<mydim, cdim, GridImp >::Jacobian&
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  jacobianInverseTransposed (const LocalCoordinate & local) const
  {
    return geoImpl_.mapping().jacobianInverseTransposed( local );
  }

  template< int mydim, int cdim, class GridImp>
  inline const typename ALU3dGridGeometry<mydim, cdim, GridImp >::JacobianTransposed&
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  jacobianTransposed (const LocalCoordinate & local) const
  {
    return geoImpl_.mapping().jacobianTransposed( local );
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
    geoImpl_.updateInFather( fatherGeom, myGeom );

    // my volume is a part of 1 for hexas, for tetra adjust with factor
    volume_ = myGeom.volume() / fatherGeom.volume();
    if( elementType == tetra )
    {
      volume_ /= 6.0;
#ifndef NDEBUG
      LocalCoordinate local( 0.0 );
      assert( std::abs( 6.0 * volume_ - integrationElement( local ) ) < 1e-12 );
#endif
    }

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
      geoImpl_.update( item.myvertex(0)->Point(),
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
      geoImpl_.update( item.myvertex(0)->Point(),
                       item.myvertex(1)->Point(),
                       item.myvertex(2)->Point(),
                       item.myvertex(3)->Point() );
    }

    // get volume of element
    volume_ = item.volume();

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

    // if face was not set, set it to zero
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
      geoImpl_.update( face.myvertex(rotatedALUIndex[0])->Point(),
                       face.myvertex(rotatedALUIndex[1])->Point(),
                       face.myvertex(rotatedALUIndex[2])->Point(),
                       face.myvertex(rotatedALUIndex[3])->Point() );
    }
    else if ( elementType == tetra )
    {
      // update geometry implementation
      geoImpl_.update( face.myvertex(rotatedALUIndex[0])->Point(),
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
    geoImpl_.update( p0, p1, p2, p3 );
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
    geoImpl_.update( p0, p1, p2 );
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
    geoImpl_.update( edge.myvertex((twist)  %2)->Point(),
                     edge.myvertex((1+twist)%2)->Point() );
    return true;
  }

  template <int mydim, int cdim, class GridImp> // for Vertices ,i.e. Points
  inline bool
  ALU3dGridGeometry<mydim, cdim, GridImp >::
  buildGeom(const VertexType & item, int twist, int)
  {
    // update geometry implementation
    geoImpl_.update( static_cast<const GEOVertexType &> (item).Point() );
    return true;
  }

} // end namespace Dune
#endif // end DUNE_ALUGRID_GEOMETRY_IMP_CC
