// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDFACEUTILITY_HH
#define DUNE_ALU3DGRIDFACEUTILITY_HH

#include <dune/common/misc.hh>
#include <dune/grid/common/genericreferenceelements.hh>

#include "mappings.hh"
#include "alu3dinclude.hh"
#include "topology.hh"

namespace Dune
{

  // convert FieldVectors to alu3dtypes
  // only used for calculating the normals because the method of the
  // mapping classes want double (&)[3] and we have FieldVectors which store an
  // double [3] this is why we can cast here
  // plz say notin' Adrian
  template< int dim >
  inline alu3d_ctype (&fieldVector2alu3d_ctype ( FieldVector< alu3d_ctype, dim > &val ))[ dim ]
  {
    return ((alu3d_ctype (&)[dim])(*( &(val[0])) ));
  }

  // convert const FieldVectors to const alu3dtypes
  template< int dim >
  inline const alu3d_ctype (&fieldVector2alu3d_ctype ( const FieldVector< alu3d_ctype, dim > &val ))[ dim ]
  {
    return ((const alu3d_ctype (&)[dim])(*( &(val[0])) ) );
  }


  // * Note: reconsider lazy evaluation of coordinates

  //- class ALU3dGridFaceInfo
  /* \brief Stores face and adjoining elements of the underlying ALU3dGrid
     The class has the same notion of inner and outer element as the
     intersection iterator.
   */
  template< ALU3dGridElementType type, class Comm >
  class ALU3dGridFaceInfo
  {
    typedef ALU3dImplTraits< type, Comm >  ImplTraits;
    //- private typedefs
    typedef typename ImplTraits::HasFaceType HasFaceType;
  public:
    enum ConformanceState {CONFORMING, REFINED_INNER, REFINED_OUTER, UNDEFINED };
    //- typedefs
    typedef typename ImplTraits::GEOFaceType GEOFaceType;
    typedef typename ImplTraits::GEOElementType GEOElementType;
    typedef typename ImplTraits::GEOPeriodicType GEOPeriodicType;
    typedef typename ImplTraits::IMPLElementType IMPLElementType;
    typedef typename ImplTraits::GhostPairType GhostPairType;
    typedef typename ImplTraits::BNDFaceType BNDFaceType;

  public:
    //! constructor creating empty face info
    ALU3dGridFaceInfo();
    void updateFaceInfo(const GEOFaceType& face, int innerLevel, int innerTwist);

    //- constructors and destructors
    //! Construct a connector from a face and the twist seen from the inner
    //! element
    //! \note: The user is responsible for the consistency of the input data
    //! as well as for choosing the appropriate (i.e. most refined) face
    ALU3dGridFaceInfo(const GEOFaceType& face, int innerTwist);
    //! Copy constructor
    ALU3dGridFaceInfo(const ALU3dGridFaceInfo &orig);
    //! Destructor
    ~ALU3dGridFaceInfo();

    //- queries
    //! Lies the face on an outer boundary?
    bool outerBoundary() const;
    //! Lies the face on an internal boundary
    bool boundary() const;

    bool periodicBoundary() const;
    //! is the neighbour element a ghost elemenet or a ghost face
    //! ic case of face true is returned
    bool ghostBoundary () const;

    //! Returns the ALU3dGrid face
    const GEOFaceType& face() const;
    //! Returns the inner element at that face
    const GEOElementType& innerEntity() const;
    //! Returns the outer element at that face
    //! \note This function is only meaningful in the interior
    const GEOElementType& outerEntity() const;
    //! Returns the inner element at that face
    //! \note This function is only meaningful at a boundary
    const BNDFaceType& innerFace() const;
    //! Returns the boundary (outer) element at that face
    //! \note This function is only meaningful at a boundary
    const BNDFaceType& boundaryFace() const;

    //! Twist of the face seen from the inner element
    int innerTwist() const;
    //! Twist of the face seen from the outer element
    int outerTwist() const;

    //! Twist of the face seen from the inner element
    int duneTwist(const int faceIdx, const int aluTwist) const;

    //! Local number of the face in inner element (ALU3dGrid reference element)
    int innerALUFaceIndex() const;
    //! Local number of the face in outer element (ALU3dGrid reference element)
    int outerALUFaceIndex() const;

    //! Description of conformance on the face
    ConformanceState conformanceState() const;

  private:
    //! Description of conformance on the face
    ConformanceState getConformanceState(const int innerLevel) const;

    //- forbidden methods
    const ALU3dGridFaceInfo &
    operator=(const ALU3dGridFaceInfo &orig);

  private:
    //- member data
    const GEOFaceType* face_;
    const HasFaceType* innerElement_;
    const HasFaceType* outerElement_;

    int innerFaceNumber_;
    int outerFaceNumber_;

    int innerTwist_;
    int outerTwist_;

    bool outerBoundary_;
    bool ghostBoundary_;
    bool innerBoundary_;
    bool periodicBoundary_;

    ConformanceState conformanceState_;
  };


  // ALU3dGridSurfaceMappingFactory
  // ------------------------------

  template< ALU3dGridElementType type, class Comm >
  struct ALU3dGridSurfaceMappingFactory;

  template< class Comm >
  struct ALU3dGridSurfaceMappingFactory< tetra, Comm >
  {
    // this is the original ALUGrid LinearSurfaceMapping,
    // see mapp_tetra_3d.* in ALUGrid code
    typedef ALU3DSPACE LinearSurfaceMapping SurfaceMappingType;
    typedef typename ALU3dGridFaceInfo< tetra, Comm >::GEOFaceType GEOFaceType;

    static const int numVerticesPerFace = EntityCount< tetra >::numVerticesPerFace;

    typedef FieldMatrix< alu3d_ctype, numVerticesPerFace, 3 > CoordinateType;

    // old method, copies values for tetra twice
    SurfaceMappingType *buildSurfaceMapping ( const CoordinateType &coords ) const;
    // get face but doesn't copy values twice
    SurfaceMappingType *buildSurfaceMapping ( const GEOFaceType &face ) const;
  };

  template< class Comm >
  struct ALU3dGridSurfaceMappingFactory< hexa, Comm >
  {
    typedef BilinearSurfaceMapping SurfaceMappingType;
    typedef typename ALU3dGridFaceInfo< hexa, Comm >::GEOFaceType GEOFaceType;

    static const int numVerticesPerFace = EntityCount< hexa >::numVerticesPerFace;

    typedef FieldMatrix< alu3d_ctype, numVerticesPerFace, 3 > CoordinateType;

    // old method, copies values for tetra twice
    SurfaceMappingType *buildSurfaceMapping ( const CoordinateType &coords ) const;
    // get face but doesn't copy values twice
    SurfaceMappingType *buildSurfaceMapping ( const GEOFaceType &face ) const;
  };



  // ALU3dGridGeometricFaceInfoBase
  // ------------------------------

  //! Helper class which provides geometric face information for the
  //! ALU3dGridIntersectionIterator
  template< ALU3dGridElementType type, class Comm >
  class ALU3dGridGeometricFaceInfoBase
    : public ALU3dGridSurfaceMappingFactory< type, Comm >
  {
    typedef ALU3dGridSurfaceMappingFactory< type, Comm > Base;

  public:
    typedef ElementTopologyMapping<type> ElementTopo;
    typedef FaceTopologyMapping<type> FaceTopo;
    typedef NonConformingFaceMapping< type, Comm > NonConformingMappingType;

    // type of container for reference elements
    typedef GenericReferenceElements< alu3d_ctype, 3 > ReferenceElementContainerType;
    // type of container for reference faces
    typedef GenericReferenceElements< alu3d_ctype, 2 > ReferenceFaceContainerType;

    // type of reference element
    typedef GenericReferenceElement<alu3d_ctype, 3> ReferenceElementType;
    // type of reference face
    typedef GenericReferenceElement<alu3d_ctype, 2> ReferenceFaceType;

    enum SideIdentifier { INNER, OUTER };
    enum { dimworld = 3 }; // ALU is a pure 3d grid
    enum { numVerticesPerFace =
             EntityCount<type>::numVerticesPerFace };

    //- public typedefs
    typedef FieldVector<alu3d_ctype, 3> NormalType;
    typedef FieldMatrix<alu3d_ctype,
        numVerticesPerFace,
        dimworld> CoordinateType;

    typedef typename ALU3dGridFaceInfo< type, Comm >::GEOFaceType GEOFaceType;

  public:
    typedef ALU3dGridFaceInfo< type, Comm > ConnectorType;

    //- constructors and destructors
    ALU3dGridGeometricFaceInfoBase(const ConnectorType &);
    ALU3dGridGeometricFaceInfoBase(const ALU3dGridGeometricFaceInfoBase &);

    //! reset status of faceGeomInfo
    void resetFaceGeom();

    //- functions
    const CoordinateType& intersectionSelfLocal() const;
    const CoordinateType& intersectionNeighborLocal() const;

  private:
    //- forbidden methods
    const ALU3dGridGeometricFaceInfoBase &operator=(const ALU3dGridGeometricFaceInfoBase &);

  private:
    //- private methods
    void generateLocalGeometries() const;

    int globalVertexIndex(int duneFaceIndex,
                          int faceTwist,
                          int duneFaceVertexIndex) const;

    void referenceElementCoordinatesRefined(SideIdentifier side,
                                            CoordinateType& result) const;
    void referenceElementCoordinatesUnrefined(SideIdentifier side,
                                              CoordinateType& result) const;

    void convert2CArray(const FieldVector<alu3d_ctype, 3>& in,
                        alu3d_ctype (&out)[3]) const;
    void convert2FieldVector(const alu3d_ctype (&in)[3],
                             FieldVector<alu3d_ctype, 3>& out) const;
  protected:
    //- private data
    const ConnectorType& connector_;

    mutable CoordinateType coordsSelfLocal_;
    mutable CoordinateType coordsNeighborLocal_;

    mutable bool generatedGlobal_;
    mutable bool generatedLocal_;

    inline static const ReferenceElementType& getReferenceElement()
    {
      return (type == tetra) ?
             ReferenceElementContainerType :: simplex() :
             ReferenceElementContainerType :: cube();
    }

    inline static const ReferenceFaceType& getReferenceFace()
    {
      return (type == tetra) ?
             ReferenceFaceContainerType :: simplex() :
             ReferenceFaceContainerType :: cube();
    }
  };

  //! Helper class which provides geometric face information for the
  //! ALU3dGridIntersectionIterator
  template< class Comm >
  class ALU3dGridGeometricFaceInfoTetra
    : public ALU3dGridGeometricFaceInfoBase< tetra, Comm >
  {
    typedef ALU3dGridGeometricFaceInfoBase< tetra, Comm > Base;

  public:
    //- public typedefs
    typedef FieldVector<alu3d_ctype, 3> NormalType;
    typedef typename Base::FaceTopo FaceTopo;
    typedef typename ALU3dGridFaceInfo< tetra, Comm >::GEOFaceType GEOFaceType;

    typedef ALU3dGridFaceInfo< tetra, Comm > ConnectorType;

    //- constructors and destructors
    ALU3dGridGeometricFaceInfoTetra(const ConnectorType& ctor);
    ALU3dGridGeometricFaceInfoTetra(const ALU3dGridGeometricFaceInfoTetra & orig);

    NormalType & outerNormal(const FieldVector<alu3d_ctype, 2>& local) const;

    //! reset status of faceGeomInfo
    void resetFaceGeom();

    //! update global geometry
    template <class GeometryImp>
    void buildGlobalGeom(GeometryImp& geo) const;

  private:
    //- forbidden methods
    const ALU3dGridGeometricFaceInfoTetra & operator=(const ALU3dGridGeometricFaceInfoTetra &);

  protected:
    using Base::connector_;

  private:
    //- private data
    mutable NormalType outerNormal_;

    // false if surface mapping needs a update
    mutable bool normalUp2Date_;
  };

  //! Helper class which provides geometric face information for the
  //! ALU3dGridIntersectionIterator
  template< class Comm >
  class ALU3dGridGeometricFaceInfoHexa
    : public ALU3dGridGeometricFaceInfoBase< hexa, Comm >
  {
    typedef ALU3dGridGeometricFaceInfoBase< hexa, Comm > Base;

  public:
    //- public typedefs
    typedef FieldVector<alu3d_ctype, 3> NormalType;
    typedef typename Base::FaceTopo FaceTopo;
    typedef typename ALU3dGridFaceInfo< hexa, Comm >::GEOFaceType GEOFaceType;
    typedef SurfaceNormalCalculator SurfaceMappingType;

    typedef ALU3dGridFaceInfo< hexa, Comm > ConnectorType;

    //- constructors and destructors
    ALU3dGridGeometricFaceInfoHexa(const ConnectorType &);
    ALU3dGridGeometricFaceInfoHexa(const ALU3dGridGeometricFaceInfoHexa &);

    NormalType & outerNormal(const FieldVector<alu3d_ctype, 2>& local) const;

    //! reset status of faceGeomInfo
    void resetFaceGeom();

    //! update global geometry
    template <class GeometryImp>
    void buildGlobalGeom(GeometryImp& geo) const;

  private:
    //- forbidden methods
    const ALU3dGridGeometricFaceInfoHexa & operator=(const ALU3dGridGeometricFaceInfoHexa &);

  protected:
    using Base::connector_;

  private:
    //- private data
    mutable NormalType outerNormal_;

    // surface mapping for calculating the outer normal
    mutable SurfaceMappingType mappingGlobal_;

    // false if surface mapping needs a update
    mutable bool mappingGlobalUp2Date_;
  };

} // end namespace Dune

#include "faceutility_imp.cc"

#endif
