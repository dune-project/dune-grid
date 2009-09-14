// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FACEUTILITY_IMP_HH
#define DUNE_FACEUTILITY_IMP_HH

namespace Dune {


  template <ALU3dGridElementType type>
  inline ALU3dGridFaceInfo<type>::ALU3dGridFaceInfo() :
    face_(0),
    innerElement_(0),
    outerElement_(0),
    innerFaceNumber_(-1),
    outerFaceNumber_(-1),
    innerTwist_(-665),
    outerTwist_(-665),
    outerBoundary_  ( false ),
    ghostBoundary_  ( false ),
    innerBoundary_  ( false ),
    conformanceState_(UNDEFINED)
  {}

  // points face from inner element away?
  template <ALU3dGridElementType type>
  inline void
  ALU3dGridFaceInfo<type>::updateFaceInfo(const GEOFaceType& face,
                                          int innerLevel,
                                          int innerTwist)
  {
    face_ = &face;
    innerElement_ = 0;
    outerElement_ = 0;
    innerFaceNumber_ = -1;
    outerFaceNumber_ = -1;
    outerBoundary_   = false;
    ghostBoundary_   = false;

    // points face from inner element away?
    if (innerTwist < 0)
    {
      innerElement_ = face.nb.rear().first;
      innerFaceNumber_ = face.nb.rear().second;
      outerElement_ = face.nb.front().first;
      outerFaceNumber_ = face.nb.front().second;
    }
    else
    {
      innerElement_ = face.nb.front().first;
      innerFaceNumber_ = face.nb.front().second;
      outerElement_ = face.nb.rear().first;
      outerFaceNumber_ = face.nb.rear().second;
    } // end if

    // if not true we are accessing a fake bnd
    assert( innerElement_->isRealObject() );
    // if not true we are accessing a fake bnd
    assert( outerElement_->isRealObject() );

#if ALU3DGRID_PARALLEL
    innerBoundary_ = innerElement_->isboundary();
    if( innerBoundary_ )
    {
      // check for ghosts
      // this check is only need in the parallel case
      const BndFaceType * bnd = dynamic_cast<const BndFaceType *> (innerElement_);

      if(bnd->bndtype() == ALU3DSPACE ProcessorBoundary_t)
      {
        // if nonconformity occurs then go up one level
        if( bnd->level () != bnd->ghostLevel() )
        {
          bnd = static_cast<const BndFaceType *>(bnd->up());
          assert( bnd );
          innerElement_ = dynamic_cast<const HasFaceType*> (bnd);
        }

        // get ghost and internal number
        GhostPairType p  = bnd->getGhost();

        innerFaceNumber_ = p.second;

        // this doesn't count as outer boundary
        const GEOElementType* ghost = static_cast<const GEOElementType*> (p.first);
        assert(ghost);

        innerTwist_ = ghost->twist(innerFaceNumber_);
      }
      else
      {
        innerTwist_ = innerFace().twist(innerALUFaceIndex());
      }
    }
    else
    {
      assert(innerTwist == innerEntity().twist(innerFaceNumber_));
      innerTwist_ = innerTwist;
    }

#endif

    outerBoundary_ = outerElement_->isboundary();
    if( outerBoundary_ )
    {
#if ALU3DGRID_PARALLEL
      // check for ghosts
      // this check is only need in the parallel case
      const BndFaceType * bnd = dynamic_cast<const BndFaceType *> (outerElement_);

      if(bnd->bndtype() == ALU3DSPACE ProcessorBoundary_t)
      {
        // if nonconformity occurs then go up one level
        if( bnd->level () != bnd->ghostLevel() )
        {
          bnd = static_cast<const BndFaceType *>(bnd->up());
          assert( bnd );
          outerElement_ = dynamic_cast<const HasFaceType*> (bnd);
        }

        // get ghost and internal number
        GhostPairType p  = bnd->getGhost();

        outerFaceNumber_ = p.second;
        ghostBoundary_   = true;

        // this doesn't count as outer boundary
        outerBoundary_ = false;
        const GEOElementType* ghost = static_cast<const GEOElementType*> (p.first);
        assert(ghost);

        outerTwist_ = ghost->twist(outerFaceNumber_);
      }
      else
#endif
      {
        outerTwist_ = boundaryFace().twist(outerALUFaceIndex());
      }
    }
    else
    {
      outerTwist_ = outerEntity().twist(outerALUFaceIndex());
    }

    // set conformance information
    conformanceState_ = getConformanceState(innerLevel);
  }

  // points face from inner element away?
  template <ALU3dGridElementType type>
  inline ALU3dGridFaceInfo<type>::ALU3dGridFaceInfo(const GEOFaceType& face,
                                                    int innerTwist)
  {
    updateFaceInfo(face,innerTwist);
  }

  template <ALU3dGridElementType type>
  inline ALU3dGridFaceInfo<type>::~ALU3dGridFaceInfo() {}

  template <ALU3dGridElementType type>
  ALU3dGridFaceInfo<type>::
  ALU3dGridFaceInfo(const ALU3dGridFaceInfo<type>& orig) :
    face_(orig.face_),
    innerElement_(orig.innerElement_),
    outerElement_(orig.outerElement_),
    innerFaceNumber_(orig.innerFaceNumber_),
    outerFaceNumber_(orig.outerFaceNumber_),
    innerTwist_(orig.innerTwist_),
    outerTwist_(orig.outerTwist_),
    outerBoundary_(orig.outerBoundary_),
    ghostBoundary_(orig.ghostBoundary_),
    innerBoundary_(orig.innerBoundary_),
    conformanceState_(orig.conformanceState_)
  {}

  template <ALU3dGridElementType type>
  inline bool ALU3dGridFaceInfo<type>::outerBoundary() const {
    return outerBoundary_;
  }

  template <ALU3dGridElementType type>
  inline bool ALU3dGridFaceInfo<type>::boundary() const {
    return outerBoundary_ || ghostBoundary_;
  }

  template <ALU3dGridElementType type>
  inline bool ALU3dGridFaceInfo<type>::ghostBoundary () const {
    return ghostBoundary_;
  }

  template <ALU3dGridElementType type>
  inline const typename ALU3dGridFaceInfo<type>::GEOFaceType&
  ALU3dGridFaceInfo<type>::face() const {
    return *face_;
  }

  template <ALU3dGridElementType type>
  inline const typename ALU3dGridFaceInfo<type>::GEOElementType&
  ALU3dGridFaceInfo<type>::innerEntity() const {
    assert(!innerElement_->isboundary());
    return static_cast<const GEOElementType&>(*innerElement_);
  }

  template <ALU3dGridElementType type>
  inline const typename ALU3dGridFaceInfo<type>::GEOElementType&
  ALU3dGridFaceInfo<type>::outerEntity() const {
    assert(!boundary());
    return static_cast<const GEOElementType&>(*outerElement_);
  }

  template <ALU3dGridElementType type>
  inline const typename ALU3dGridFaceInfo<type>::BndFaceType&
  ALU3dGridFaceInfo<type>::innerFace() const {
    assert( innerElement_->isboundary() );
    return static_cast<const BndFaceType&>(*innerElement_);
  }

  template <ALU3dGridElementType type>
  inline const typename ALU3dGridFaceInfo<type>::BndFaceType&
  ALU3dGridFaceInfo<type>::boundaryFace() const {
    assert(boundary());
    return static_cast<const BndFaceType&>(*outerElement_);
  }

  template <ALU3dGridElementType type>
  inline int ALU3dGridFaceInfo<type>::innerTwist() const
  {
    // don't check ghost boundaries here
    assert( ( ! innerBoundary_) ?
            innerEntity().twist(innerALUFaceIndex()) == innerTwist_ : true );
    return innerTwist_;
  }

  template <ALU3dGridElementType type>
  inline int ALU3dGridFaceInfo<type>::outerTwist() const
  {
    // don't check ghost boundaries here
    assert( (outerBoundary_) ?
            (outerTwist_ == boundaryFace().twist(0)) :
            (! ghostBoundary_) ?
            (outerTwist_ == outerEntity().twist(outerALUFaceIndex())) : true
            );
    return outerTwist_;
  }

  template <ALU3dGridElementType type>
  inline int ALU3dGridFaceInfo<type>::innerALUFaceIndex() const {
    return innerFaceNumber_;
  }

  template <ALU3dGridElementType type>
  inline int ALU3dGridFaceInfo<type>::outerALUFaceIndex() const {
    return outerFaceNumber_;
  }

  template <ALU3dGridElementType type>
  typename ALU3dGridFaceInfo<type>::ConformanceState
  inline ALU3dGridFaceInfo<type>::conformanceState() const
  {
    assert( conformanceState_ != UNDEFINED );
    return conformanceState_;
  }

  // calculate conformance state
  template <ALU3dGridElementType type>
  typename ALU3dGridFaceInfo<type>::ConformanceState
  inline ALU3dGridFaceInfo<type>::getConformanceState(const int innerLevel) const
  {
    ConformanceState result = CONFORMING;

    // A boundary is always unrefined
    int levelDifference = 0 ;
    if ( boundary() )
      levelDifference = innerLevel - boundaryFace().level();
    else
      levelDifference = innerLevel - outerEntity().level();

    if (levelDifference < 0) {
      result = REFINED_OUTER;
    }
    else if (levelDifference > 0) {
      result = REFINED_INNER;
    }

    return result;
  }

  template <ALU3dGridElementType type>
  inline ALU3dGridGeometricFaceInfoBase<type>::
  ALU3dGridGeometricFaceInfoBase(const ConnectorType& connector) :
    connector_(connector),
    generatedGlobal_(false),
    generatedLocal_(false),
    coordsSelfLocal_(-1.0),
    coordsNeighborLocal_(-1.0),
    refElem_( getReferenceElement() ),
    refFace_( getReferenceFace() )
  {}

  template <ALU3dGridElementType type>
  inline void
  ALU3dGridGeometricFaceInfoBase<type>::
  resetFaceGeom()
  {
    generatedGlobal_ = false;
    generatedLocal_  = false;
  }

  template <ALU3dGridElementType type>
  inline ALU3dGridGeometricFaceInfoBase<type>::
  ALU3dGridGeometricFaceInfoBase(const ALU3dGridGeometricFaceInfoBase<type>& orig) :
    connector_(orig.connector_),
    generatedGlobal_(orig.generatedGlobal_),
    generatedLocal_(orig.generatedLocal_),
    coordsSelfLocal_(orig.coordsSelfLocal_),
    coordsNeighborLocal_(orig.coordsNeighborLocal_),
    refElem_( getReferenceElement() ),
    refFace_( getReferenceFace() )
  {}

  template <ALU3dGridElementType type>
  inline const typename ALU3dGridGeometricFaceInfoBase<type>::CoordinateType&
  ALU3dGridGeometricFaceInfoBase<type>::intersectionSelfLocal() const {
    generateLocalGeometries();
    assert(generatedLocal_);
    return coordsSelfLocal_;
  }

  template <ALU3dGridElementType type>
  inline const typename ALU3dGridGeometricFaceInfoBase<type>::CoordinateType&
  ALU3dGridGeometricFaceInfoBase<type>::intersectionNeighborLocal() const {
    assert(!connector_.outerBoundary());
    generateLocalGeometries();
    assert(generatedLocal_);
    return coordsNeighborLocal_;
  }


  //sepcialisation for tetra and hexa
  inline ALU3dGridGeometricFaceInfoTetra::
  ALU3dGridGeometricFaceInfoTetra(const ConnectorType& connector)
    : ALU3dGridGeometricFaceInfoBase<tetra>(connector)
      , normalUp2Date_(false)
  {}

  inline void ALU3dGridGeometricFaceInfoTetra::
  resetFaceGeom()
  {
    ALU3dGridGeometricFaceInfoBase<tetra>::resetFaceGeom();
    normalUp2Date_ = false;
  }

  inline ALU3dGridGeometricFaceInfoTetra::
  ALU3dGridGeometricFaceInfoTetra(const ALU3dGridGeometricFaceInfoTetra& orig)
    : ALU3dGridGeometricFaceInfoBase<tetra>(orig)
      , normalUp2Date_(orig.normalUp2Date_)
  {}

  template <class GeometryImp>
  inline void
  ALU3dGridGeometricFaceInfoTetra::
  buildGlobalGeom(GeometryImp& geo) const
  {
    if (! this->generatedGlobal_)
    {
      // calculate the normal
      const GEOFaceType & face = this->connector_.face();

      geo.buildGeom( face.myvertex(FaceTopo::dune2aluVertex(0))->Point() ,
                     face.myvertex(FaceTopo::dune2aluVertex(1))->Point() ,
                     face.myvertex(FaceTopo::dune2aluVertex(2))->Point() );
      this->generatedGlobal_ = true ;
    }
  }

  inline FieldVector<alu3d_ctype, 3> &
  ALU3dGridGeometricFaceInfoTetra::
  outerNormal(const FieldVector<alu3d_ctype, 2>& local) const
  {
    // if geomInfo was not reseted then normal is still correct
    if(!normalUp2Date_)
    {
      // calculate the normal
      const GEOFaceType & face = this->connector_.face();
      const double (&_p0)[3] = face.myvertex(0)->Point();
      const double (&_p1)[3] = face.myvertex(1)->Point();
      const double (&_p2)[3] = face.myvertex(2)->Point();

      // change sign if face normal points into inner element
      // factor is 1.0 to get integration outer normal and not volume outer normal
      const double factor = (this->connector_.innerTwist() < 0) ? 1.0 : -1.0;

      // see mapp_tetra_3d.h for this piece of code
      outerNormal_[0] = factor * ((_p1[1]-_p0[1]) *(_p2[2]-_p1[2]) - (_p2[1]-_p1[1]) *(_p1[2]-_p0[2]));
      outerNormal_[1] = factor * ((_p1[2]-_p0[2]) *(_p2[0]-_p1[0]) - (_p2[2]-_p1[2]) *(_p1[0]-_p0[0]));
      outerNormal_[2] = factor * ((_p1[0]-_p0[0]) *(_p2[1]-_p1[1]) - (_p2[0]-_p1[0]) *(_p1[1]-_p0[1]));

      normalUp2Date_ = true;
    } // end if mapp ...

    return outerNormal_;
  }

  //-sepcialisation for and hexa
  inline ALU3dGridGeometricFaceInfoHexa::
  ALU3dGridGeometricFaceInfoHexa(const ConnectorType& connector)
    : ALU3dGridGeometricFaceInfoBase<hexa>(connector)
      , mappingGlobal_()
      , mappingGlobalUp2Date_(false)
  {}

  inline void ALU3dGridGeometricFaceInfoHexa::
  resetFaceGeom()
  {
    ALU3dGridGeometricFaceInfoBase<hexa>::resetFaceGeom();
    mappingGlobalUp2Date_ = false;
  }

  inline ALU3dGridGeometricFaceInfoHexa::
  ALU3dGridGeometricFaceInfoHexa(const ALU3dGridGeometricFaceInfoHexa& orig)
    : ALU3dGridGeometricFaceInfoBase<hexa>(orig)
      , mappingGlobal_(orig.mappingGlobal_)
      , mappingGlobalUp2Date_(orig.mappingGlobalUp2Date_)
  {}

  template <class GeometryImp>
  inline void
  ALU3dGridGeometricFaceInfoHexa::
  buildGlobalGeom(GeometryImp& geo) const
  {
    if (! this->generatedGlobal_)
    {
      // calculate the normal
      const GEOFaceType & face = this->connector_.face();

      geo.buildGeom( face.myvertex(FaceTopo::dune2aluVertex(0))->Point() ,
                     face.myvertex(FaceTopo::dune2aluVertex(1))->Point() ,
                     face.myvertex(FaceTopo::dune2aluVertex(2))->Point() ,
                     face.myvertex(FaceTopo::dune2aluVertex(3))->Point() );
      this->generatedGlobal_ = true ;
    }
  }

  inline FieldVector<alu3d_ctype, 3> &
  ALU3dGridGeometricFaceInfoHexa::
  outerNormal(const FieldVector<alu3d_ctype, 2>& local) const
  {
    // if mapping calculated and affine, nothing more to do
    if ( mappingGlobal_.affine () && mappingGlobalUp2Date_ )
      return outerNormal_ ;

    // update surface mapping
    if(! mappingGlobalUp2Date_ )
    {
      const GEOFaceType & face = connector_.face();
      // update mapping to actual face
      mappingGlobal_.buildMapping(
        face.myvertex( FaceTopo::dune2aluVertex(0) )->Point(),
        face.myvertex( FaceTopo::dune2aluVertex(1) )->Point(),
        face.myvertex( FaceTopo::dune2aluVertex(2) )->Point(),
        face.myvertex( FaceTopo::dune2aluVertex(3) )->Point()
        );
      mappingGlobalUp2Date_ = true;
    }

    // calculate the normal
    // has to be calculated every time normal called, because
    // depends on local
    if (connector_.innerTwist() < 0)
      mappingGlobal_.negativeNormal(local,outerNormal_);
    else
      mappingGlobal_.normal(local,outerNormal_);

    // end if
    return outerNormal_;
  }

  template <ALU3dGridElementType type>
  inline void ALU3dGridGeometricFaceInfoBase<type>::
  generateLocalGeometries() const
  {
    if (!generatedLocal_) {
      // Get the coordinates of the face in the reference element of the
      // adjoining inner and outer elements and initialise the respective
      // geometries
      switch (connector_.conformanceState()) {
      case (ConnectorType::CONFORMING) :
        referenceElementCoordinatesRefined(INNER, coordsSelfLocal_);
        // generate outer local geometry only when not at boundary
        // * in the parallel case, this needs to be altered for the ghost cells
        if (!connector_.outerBoundary()) {
          referenceElementCoordinatesRefined(OUTER, coordsNeighborLocal_);
        } // end if
        break;
      case (ConnectorType::REFINED_INNER) :
        referenceElementCoordinatesRefined(INNER, coordsSelfLocal_);
        referenceElementCoordinatesUnrefined(OUTER, coordsNeighborLocal_);
        break;
      case (ConnectorType::REFINED_OUTER) :
        referenceElementCoordinatesUnrefined(INNER, coordsSelfLocal_);
        referenceElementCoordinatesRefined(OUTER, coordsNeighborLocal_);
        break;
      default :
        std::cerr << "ERROR: Wrong conformanceState in generateLocalGeometries! in: " << __FILE__ << " line: " << __LINE__<< std::endl;
        assert(false);
        exit(1);
      } // end switch

      generatedLocal_ = true;
    } // end if
  }

  template <ALU3dGridElementType type>
  inline int ALU3dGridGeometricFaceInfoBase<type>::
  globalVertexIndex(int duneFaceIndex,
                    int faceTwist,
                    int duneFaceVertexIndex) const {

    int localALUIndex = FaceTopo::dune2aluVertex(duneFaceVertexIndex,
                                                 faceTwist);

    // get local ALU vertex number on the element's face
    int localDuneIndex = ElementTopo::
                         alu2duneFaceVertex(ElementTopo::dune2aluFace(duneFaceIndex),
                                            localALUIndex);

    return refElem_.subEntity(duneFaceIndex, 1, localDuneIndex, 3);
  }

  template <ALU3dGridElementType type>
  inline void ALU3dGridGeometricFaceInfoBase<type>::
  referenceElementCoordinatesRefined(SideIdentifier side,
                                     CoordinateType& result) const
  {
    // this is a dune face index
    int faceIndex =
      (side == INNER ?
       ElementTopo::alu2duneFace(connector_.innerALUFaceIndex()) :
       ElementTopo::alu2duneFace(connector_.outerALUFaceIndex()));
    int faceTwist =
      (side == INNER ?
       connector_.innerTwist() :
       connector_.outerTwist());

    for (int i = 0; i < numVerticesPerFace; ++i) {
      int duneVertexIndex = globalVertexIndex(faceIndex, faceTwist, i);
      result[i] = refElem_.position(duneVertexIndex, 3);
    }
  }

  template <ALU3dGridElementType type>
  inline void ALU3dGridGeometricFaceInfoBase<type>::
  convert2CArray(const FieldVector<alu3d_ctype, 3>& in,
                 alu3d_ctype (&out)[3]) const {
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
  }

  template <ALU3dGridElementType type>
  inline void ALU3dGridGeometricFaceInfoBase<type>::
  convert2FieldVector(const alu3d_ctype (&in)[3],
                      FieldVector<alu3d_ctype, 3>& out) const {
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
  }

} //end namespace Dune
#endif
