// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FACEUTILITY_IMP_HH
#define DUNE_FACEUTILITY_IMP_HH

namespace Dune
{


  template< ALU3dGridElementType type, class Comm >
  inline ALU3dGridFaceInfo< type, Comm >::ALU3dGridFaceInfo() :
    face_(0),
    innerElement_(0),
    outerElement_(0),
    innerFaceNumber_(-1),
    outerFaceNumber_(-1),
    innerTwist_(-665),
    outerTwist_(-665),
    segmentIndex_( -1 ),
    bndType_( noBoundary ),
    conformanceState_(UNDEFINED)
  {}

  // points face from inner element away?
  template< ALU3dGridElementType type, class Comm >
  inline void
  ALU3dGridFaceInfo< type, Comm >::
  updateFaceInfo(const GEOFaceType& face,
                 int innerLevel,
                 int innerTwist)
  {
    face_ = &face;
    innerElement_ = 0;
    outerElement_ = 0;
    innerFaceNumber_ = -1;
    outerFaceNumber_ = -1;
    bndType_ = noBoundary;
    segmentIndex_ = -1;

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

    // we only have to do this in parallel runs
    if( parallel() && innerElement_->isboundary() )
    {
      bndType_ = innerGhostBoundary;
    }

    if( parallel() && innerBoundary() )
    {
      // check for ghosts
      // this check is only need in the parallel case
      const BNDFaceType * bnd = static_cast<const BNDFaceType *> (innerElement_);

      if(bnd->bndtype() == ALU3DSPACE ProcessorBoundary_t)
      {
        // if nonconformity occurs then go up one level
        if( bnd->level () != bnd->ghostLevel() )
        {
          bnd = static_cast<const BNDFaceType *>(bnd->up());
          assert( bnd );
          innerElement_ = static_cast<const HasFaceType*> (bnd);
        }

        // get ghost and internal number
        GhostPairType p  = bnd->getGhost();

        // get face number
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
      // set inner twist
      assert(innerTwist == innerEntity().twist(innerFaceNumber_));
      innerTwist_ = innerTwist;
    }

    if( outerElement_->isboundary() )
    {
      assert( ! innerBoundary() );
      // set to default boundary (with domain boundary)
      bndType_ = domainBoundary ;

      // check for ghosts
      // this check is only need in the parallel case
      // if this cast fails we have a periodic element
      const BNDFaceType * bnd = dynamic_cast<const BNDFaceType *> (outerElement_);

      if( ! bnd ) // the periodic case
      {
        bndType_ = periodicBoundary ;
        assert( dynamic_cast< const GEOPeriodicType* > ( outerElement_ ) );
        const GEOPeriodicType* periodicClosure = static_cast< const GEOPeriodicType* > ( outerElement_ ) ;

#ifdef ALUGRID_PERIODIC_BOUNDARY
        // previously, the segmentIndex( 1 - outerFaceNumber_ ) was used, why?
        segmentIndex_ = periodicClosure->segmentIndex( outerFaceNumber_ );
#else
        // set to zero (grid test will fail)
        segmentIndex_ = 0 ;
#endif

        const GEOFaceType* face = ImplTraits::getFace( *periodicClosure, 1 - outerFaceNumber_ );
        if( dynamic_cast< const GEOPeriodicType * >( face->nb.rear().first ) )
        //if( face->nb.rear().first->isboundary() )
        {
          //assert( !face->nb.front().first->isboundary() );
          outerElement_    = face->nb.front().first ;
          outerFaceNumber_ = face->nb.front().second ;
        }
        else
        {
          //assert( face->nb.front().first->isboundary() );
          outerElement_    = face->nb.rear().first ;
          outerFaceNumber_ = face->nb.rear().second ;
        }

        assert( outerElement_->isRealObject() );
        if( outerElement_->isboundary() )
        {
          assert( dynamic_cast< const BNDFaceType * >( outerElement_ ) );
          bnd = static_cast< const BNDFaceType * >( outerElement_ );
        }
        else
          outerTwist_ = outerEntity().twist( outerFaceNumber_ );
      }

      if ( bnd ) // the boundary case
      {
        assert( bnd );

        // if this cast is valid we have either
        // a boundary or a ghost element
        // the ghost element case
        if( parallel() && bnd->bndtype() == ALU3DSPACE ProcessorBoundary_t)
        {
          // if nonconformity occurs then go up one level
          if( bnd->level () != bnd->ghostLevel() )
          {
            bnd = static_cast<const BNDFaceType *>(bnd->up());
            assert( bnd );
            outerElement_ = static_cast<const HasFaceType*> (bnd);
          }

          // get ghost and internal number
          GhostPairType p  = bnd->getGhost();
          outerFaceNumber_ = p.second;

          // set boundary type to ghost boundary
          bndType_ = outerGhostBoundary ;

          const GEOElementType* ghost = static_cast<const GEOElementType*> (p.first);
          assert(ghost);

          outerTwist_ = ghost->twist(outerFaceNumber_);
        }
        else // the normal boundary case
        {
          // get outer twist
          outerTwist_ = boundaryFace().twist(outerALUFaceIndex());
          // store segment index
          segmentIndex_ = boundaryFace().segmentIndex();
        }
      }
    }
    else
    {
      // get outer twist
      outerTwist_ = outerEntity().twist(outerALUFaceIndex());
    }

    // set conformance information
    conformanceState_ = getConformanceState(innerLevel);
  }

  // points face from inner element away?
  template< ALU3dGridElementType type, class Comm >
  inline ALU3dGridFaceInfo< type, Comm >::
  ALU3dGridFaceInfo(const GEOFaceType& face,
                    int innerTwist)
  {
    updateFaceInfo(face,innerTwist);
  }

  template< ALU3dGridElementType type, class Comm >
  inline ALU3dGridFaceInfo< type, Comm >::~ALU3dGridFaceInfo() {}

  template< ALU3dGridElementType type, class Comm >
  ALU3dGridFaceInfo< type, Comm >::
  ALU3dGridFaceInfo ( const ALU3dGridFaceInfo &orig )
    : face_(orig.face_),
      innerElement_(orig.innerElement_),
      outerElement_(orig.outerElement_),
      innerFaceNumber_(orig.innerFaceNumber_),
      outerFaceNumber_(orig.outerFaceNumber_),
      innerTwist_(orig.innerTwist_),
      outerTwist_(orig.outerTwist_),
      segmentIndex_( orig.segmentIndex_ ),
      bndType_( orig.bndType_ ),
      conformanceState_(orig.conformanceState_)
  {}

  template< ALU3dGridElementType type, class Comm >
  inline bool ALU3dGridFaceInfo< type, Comm >::isElementLike() const {
    return bndType_ < domainBoundary;
  }

  template< ALU3dGridElementType type, class Comm >
  inline bool ALU3dGridFaceInfo< type, Comm >::innerBoundary() const {
    return bndType_ == innerGhostBoundary;
  }

  template< ALU3dGridElementType type, class Comm >
  inline bool ALU3dGridFaceInfo< type, Comm >::outerBoundary() const {
    return bndType_ == domainBoundary;
  }

  template< ALU3dGridElementType type, class Comm >
  inline bool ALU3dGridFaceInfo< type, Comm >::boundary() const {
    return outerBoundary() || (bndType_ == periodicBoundary);
  }

  template< ALU3dGridElementType type, class Comm >
  inline bool ALU3dGridFaceInfo< type, Comm >::neighbor() const
  {
    return isElementLike() || ghostBoundary();
  }

  template< ALU3dGridElementType type, class Comm >
  inline bool ALU3dGridFaceInfo< type, Comm >::ghostBoundary () const
  {
    // when communicator is No_Comm there is no ghost boundary
    return parallel() ? ( bndType_ == outerGhostBoundary ) : false ;
  }

  template< ALU3dGridElementType type, class Comm >
  inline const typename ALU3dGridFaceInfo< type, Comm >::GEOFaceType&
  ALU3dGridFaceInfo< type, Comm >::face() const
  {
    return *face_;
  }

  template< ALU3dGridElementType type, class Comm >
  inline const typename ALU3dGridFaceInfo< type, Comm >::GEOElementType&
  ALU3dGridFaceInfo< type, Comm >::innerEntity() const
  {
    assert( ! innerElement_->isboundary() );
    return static_cast<const GEOElementType&>(*innerElement_);
  }

  template< ALU3dGridElementType type, class Comm >
  inline const typename ALU3dGridFaceInfo< type, Comm >::GEOElementType&
  ALU3dGridFaceInfo< type, Comm >::outerEntity() const
  {
    assert( isElementLike() );
    return static_cast<const GEOElementType&>(*outerElement_);
  }

  template< ALU3dGridElementType type, class Comm >
  inline const typename ALU3dGridFaceInfo< type, Comm >::BNDFaceType&
  ALU3dGridFaceInfo< type, Comm >::innerFace() const
  {
    assert( innerElement_->isboundary() );
    return static_cast<const BNDFaceType&>(*innerElement_);
  }

  template< ALU3dGridElementType type, class Comm >
  inline const typename ALU3dGridFaceInfo< type, Comm >::BNDFaceType&
  ALU3dGridFaceInfo< type, Comm >::boundaryFace() const {
    assert( ! isElementLike() );
    return static_cast<const BNDFaceType&>(*outerElement_);
  }

  template< ALU3dGridElementType type, class Comm >
  inline int ALU3dGridFaceInfo< type, Comm >::outsideLevel() const
  {
    assert( outerElement_ );
    assert( !isElementLike() || outerEntity().level() == outerElement_->nbLevel() );
    assert( isElementLike() || boundaryFace().level() == outerElement_->nbLevel() );
    return outerElement_->nbLevel();
  }

  template< ALU3dGridElementType type, class Comm >
  inline int ALU3dGridFaceInfo< type, Comm >::segmentIndex() const
  {
    assert( segmentIndex_ >= 0 );
    return segmentIndex_;
  }

  template< ALU3dGridElementType type, class Comm >
  inline int ALU3dGridFaceInfo< type, Comm >::boundaryId() const
  {
    return ( outerBoundary() ) ? boundaryFace().bndtype() : 20 ;
  }
  template< ALU3dGridElementType type, class Comm >
  inline int ALU3dGridFaceInfo< type, Comm >::innerTwist() const
  {
    // don't check ghost boundaries here
    assert( ( ! innerBoundary() ) ?
            innerEntity().twist(innerALUFaceIndex()) == innerTwist_ : true );
    return innerTwist_;
  }

  template< ALU3dGridElementType type, class Comm >
  inline int ALU3dGridFaceInfo< type, Comm >::duneTwist(const int faceIdx, const int aluTwist) const
  {
    typedef ElementTopologyMapping<type> ElementTopo;
    typedef FaceTopologyMapping<type> FaceTopo;

    const int mappedZero =
      FaceTopo::twist(ElementTopo::dune2aluFaceVertex( faceIdx, 0), aluTwist);

    const int twist =
      (ElementTopo::faceOrientation( faceIdx ) * sign(aluTwist) < 0 ?
       mappedZero : -mappedZero-1);

    // see topology.* files
    return FaceTopo :: aluTwistMap( twist );
  }

  template< ALU3dGridElementType type, class Comm >
  inline int ALU3dGridFaceInfo< type, Comm >::outerTwist() const
  {
    // don't check ghost boundaries here
    //assert( (outerBoundary_) ?
    //          (outerTwist_ == boundaryFace().twist(0)) :
    //          (! ghostBoundary_) ?
    //          (outerTwist_ == outerEntity().twist(outerALUFaceIndex())) : true
    //      );
    return outerTwist_;
  }

  template< ALU3dGridElementType type, class Comm >
  inline int ALU3dGridFaceInfo< type, Comm >::innerALUFaceIndex() const {
    return innerFaceNumber_;
  }

  template< ALU3dGridElementType type, class Comm >
  inline int ALU3dGridFaceInfo< type, Comm >::outerALUFaceIndex() const {
    return outerFaceNumber_;
  }

  template< ALU3dGridElementType type, class Comm >
  typename ALU3dGridFaceInfo< type, Comm >::ConformanceState
  inline ALU3dGridFaceInfo< type, Comm >::conformanceState() const
  {
    assert( conformanceState_ != UNDEFINED );
    return conformanceState_;
  }

  // calculate conformance state
  template< ALU3dGridElementType type, class Comm >
  typename ALU3dGridFaceInfo< type, Comm >::ConformanceState
  inline ALU3dGridFaceInfo< type, Comm >::getConformanceState(const int innerLevel) const
  {
    ConformanceState result = CONFORMING;

    // A boundary is always unrefined
    int levelDifference = 0 ;
    if ( isElementLike() )
      levelDifference = innerLevel - outerEntity().level();
    else
      levelDifference = innerLevel - boundaryFace().level();

    if (levelDifference < 0) {
      result = REFINED_OUTER;
    }
    else if (levelDifference > 0) {
      result = REFINED_INNER;
    }

    return result;
  }

  template< ALU3dGridElementType type, class Comm >
  inline ALU3dGridGeometricFaceInfoBase< type, Comm >::
  ALU3dGridGeometricFaceInfoBase(const ConnectorType& connector) :
    connector_(connector),
    coordsSelfLocal_(-1.0),
    coordsNeighborLocal_(-1.0),
    generatedGlobal_(false),
    generatedLocal_(false)
  {}

  template< ALU3dGridElementType type, class Comm >
  inline void
  ALU3dGridGeometricFaceInfoBase< type, Comm >::
  resetFaceGeom()
  {
    generatedGlobal_ = false;
    generatedLocal_  = false;
  }

  template< ALU3dGridElementType type, class Comm >
  inline ALU3dGridGeometricFaceInfoBase< type, Comm >::
  ALU3dGridGeometricFaceInfoBase ( const ALU3dGridGeometricFaceInfoBase &orig )
    : connector_(orig.connector_),
      coordsSelfLocal_(orig.coordsSelfLocal_),
      coordsNeighborLocal_(orig.coordsNeighborLocal_),
      generatedGlobal_(orig.generatedGlobal_),
      generatedLocal_(orig.generatedLocal_)
  {}

  template< ALU3dGridElementType type, class Comm >
  inline const typename ALU3dGridGeometricFaceInfoBase< type, Comm >::CoordinateType&
  ALU3dGridGeometricFaceInfoBase< type, Comm >::intersectionSelfLocal() const {
    generateLocalGeometries();
    assert(generatedLocal_);
    return coordsSelfLocal_;
  }

  template< ALU3dGridElementType type, class Comm >
  inline const typename ALU3dGridGeometricFaceInfoBase< type, Comm >::CoordinateType&
  ALU3dGridGeometricFaceInfoBase< type, Comm >::intersectionNeighborLocal() const {
    assert(!connector_.outerBoundary());
    generateLocalGeometries();
    assert(generatedLocal_);
    return coordsNeighborLocal_;
  }


  //sepcialisation for tetra and hexa
  template< class Comm >
  inline ALU3dGridGeometricFaceInfoTetra< Comm >::
  ALU3dGridGeometricFaceInfoTetra(const ConnectorType& connector)
    : Base( connector ), normalUp2Date_( false )
  {}

  template< class Comm >
  inline void ALU3dGridGeometricFaceInfoTetra< Comm >::
  resetFaceGeom()
  {
    Base::resetFaceGeom();
    normalUp2Date_ = false;
  }

  template< class Comm >
  inline ALU3dGridGeometricFaceInfoTetra< Comm >::
  ALU3dGridGeometricFaceInfoTetra(const ALU3dGridGeometricFaceInfoTetra& orig)
    : Base( orig ), normalUp2Date_( orig.normalUp2Date_ )
  {}

  template< class Comm >
  template <class GeometryImp>
  inline void
  ALU3dGridGeometricFaceInfoTetra< Comm >::
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

  template< class Comm >
  inline FieldVector<alu3d_ctype, 3> &
  ALU3dGridGeometricFaceInfoTetra< Comm >::
  outerNormal(const FieldVector<alu3d_ctype, 2>& local) const
  {
    // if geomInfo was not reseted then normal is still correct
    if(!normalUp2Date_)
    {
      // calculate the normal
      const GEOFaceType & face = this->connector_.face();
      const alu3d_ctype (&_p0)[3] = face.myvertex(0)->Point();
      const alu3d_ctype (&_p1)[3] = face.myvertex(1)->Point();
      const alu3d_ctype (&_p2)[3] = face.myvertex(2)->Point();

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
  template< class Comm >
  inline ALU3dGridGeometricFaceInfoHexa< Comm >::
  ALU3dGridGeometricFaceInfoHexa(const ConnectorType& connector)
    : Base( connector )
      , mappingGlobal_()
      , mappingGlobalUp2Date_(false)
  {}

  template< class Comm >
  inline void ALU3dGridGeometricFaceInfoHexa< Comm >::
  resetFaceGeom()
  {
    Base::resetFaceGeom();
    mappingGlobalUp2Date_ = false;
  }

  template< class Comm >
  inline ALU3dGridGeometricFaceInfoHexa< Comm >::
  ALU3dGridGeometricFaceInfoHexa(const ALU3dGridGeometricFaceInfoHexa& orig)
    : Base( orig )
      , mappingGlobal_(orig.mappingGlobal_)
      , mappingGlobalUp2Date_(orig.mappingGlobalUp2Date_)
  {}

  template< class Comm >
  template <class GeometryImp>
  inline void
  ALU3dGridGeometricFaceInfoHexa< Comm >::
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

  template< class Comm >
  inline FieldVector<alu3d_ctype, 3> &
  ALU3dGridGeometricFaceInfoHexa< Comm >::
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

  template< ALU3dGridElementType type, class Comm >
  inline void ALU3dGridGeometricFaceInfoBase< type, Comm >::
  generateLocalGeometries() const
  {
    if (!generatedLocal_) {
      // Get the coordinates of the face in the reference element of the
      // adjoining inner and outer elements and initialise the respective
      // geometries
      switch (connector_.conformanceState())
      {
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

  template< ALU3dGridElementType type, class Comm >
  inline int ALU3dGridGeometricFaceInfoBase< type, Comm >::
  globalVertexIndex(const int duneFaceIndex,
                    const int aluFaceTwist,
                    const int duneFaceVertexIndex) const
  {
    const int localALUIndex =
      FaceTopo::dune2aluVertex(duneFaceVertexIndex,
                               aluFaceTwist);

    // get local ALU vertex number on the element's face
    const int localDuneIndex = ElementTopo::
                               alu2duneFaceVertex(ElementTopo::dune2aluFace(duneFaceIndex),
                                                  localALUIndex);

    return getReferenceElement().subEntity(duneFaceIndex, 1, localDuneIndex, 3);
  }


  template< ALU3dGridElementType type, class Comm >
  inline void ALU3dGridGeometricFaceInfoBase< type, Comm >::
  referenceElementCoordinatesRefined(SideIdentifier side,
                                     CoordinateType& result) const
  {
    // this is a dune face index
    const int faceIndex =
      (side == INNER ?
       ElementTopo::alu2duneFace(connector_.innerALUFaceIndex()) :
       ElementTopo::alu2duneFace(connector_.outerALUFaceIndex()));
    const int faceTwist =
      (side == INNER ?
       connector_.innerTwist() :
       connector_.outerTwist());

    const ReferenceElementType& refElem = getReferenceElement();

    for (int i = 0; i < numVerticesPerFace; ++i)
    {
      int duneVertexIndex = globalVertexIndex(faceIndex, faceTwist, i);
      result[i] = refElem.position(duneVertexIndex, 3);
    }
  }
} //end namespace Dune
#endif
