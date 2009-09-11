// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/grid/genericgeometry/conversion.hh>

#include "grid.hh"
#include "mappings.hh"

namespace Dune {
  // --Geometry

  //- Tetra specialization
  template<int mydim, int cdim>
  inline ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> >::
  ALU3dGridGeometry()
    : coord_(0.0)
      , Jinv_ ()
      , volume_(-1.0)
      , A_()
      , AT_x_()
      , localCoord_()
      , globalCoord_()
      , tmpV_()
      , tmpU_()
      , myGeomType_(GeometryType::simplex,mydim)
      , builtinverse_ (false) , builtA_ (false)
      , builtDetDF_ (false)
  {}

  //   B U I L T G E O M   - - -
  template<int mydim, int cdim>
  inline void ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> >::
  calcElMatrix () const
  {
    if(!builtA_)
    {
      // creat Matrix A (=Df)               INDIZES: row/col
      // Mapping: R^dim -> R^3,  F(x) = A x + p_0
      // columns:    p_1 - p_0  |  p_2 - p_0  |  p_3 - p_0

      for (int i=0; i<mydim; ++i) {
        for(int j=0; j<cdim; ++j)
          A_[j][i] = coord_[i+1][j] - coord_[0][j];
      }
      builtA_ = true;
    }
  }

  //dim = dimworld = 3
  template<int mydim, int cdim>
  inline void ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> > ::
  buildJacobianInverseTransposed() const
  {
    if(!builtinverse_)
    {
      calcElMatrix();

      // DetDf = integrationElement, invert transposed Matrix
      detDF_ = std::abs( FMatrixHelp::invertMatrix_retTransposed(A_,Jinv_) );
      builtinverse_ = builtDetDF_ = true;
    }
  }

  template<> //dim = 2 , dimworld = 3
  inline void ALU3dGridGeometry<2,3, const ALU3dGrid<3,3,tetra> > :: buildJacobianInverseTransposed() const
  {
    if(!builtinverse_)
    {
      // calc A and stores it in A_
      calcElMatrix();

      // calc ret = A^T*A
      FMatrixHelp::multTransposedMatrix(A_,AT_A_);

      // calc Jinv_ = A (A^T*A)^-1
      //FMatrixHelp::invertMatrix_retTransposed(AT_A_,Jinv_);
      FieldMatrix< alu3d_ctype, 2, 2 > inv_AT_A;
      FMatrixHelp :: invertMatrix( AT_A_, inv_AT_A );
      FMatrixHelp :: multMatrix( A_, inv_AT_A, Jinv_ );

      enum { dim = 3 };

      tmpV_ = coord_[1] - coord_[0];
      tmpU_ = coord_[2] - coord_[1];

      // calculate scaled outer normal
      for(int i=0; i<dim; i++)
      {
        globalCoord_[i] = (  tmpU_[(i+1)%dim] * tmpV_[(i+2)%dim]
                             - tmpU_[(i+2)%dim] * tmpV_[(i+1)%dim] );
      }

      detDF_ = std::abs ( globalCoord_.two_norm() );

      //std::cout << det << " calced| detDF " << detDF_ << "\n";

      builtinverse_ = builtDetDF_ = true;
    }
  }

  template<> //dim = 1 , dimworld = 3
  inline void ALU3dGridGeometry<1,3, const ALU3dGrid<3,3,tetra> > :: buildJacobianInverseTransposed() const
  {
    if(!builtinverse_)
    {
      // calc A and stores it in A_
      calcElMatrix();

      // calc ret = A^T*A
      FMatrixHelp::multTransposedMatrix(A_,AT_A_);

      // calc Jinv_ = A (A^T*A)^-1
      //FMatrixHelp::invertMatrix_retTransposed(AT_A_,Jinv_);
      FieldMatrix< alu3d_ctype, 1, 1 > inv_AT_A;
      FMatrixHelp :: invertMatrix( AT_A_, inv_AT_A );
      FMatrixHelp :: multMatrix( A_, inv_AT_A, Jinv_ );

      // create vectors of face
      globalCoord_ = coord_[1] - coord_[0];
      detDF_ = std::abs ( globalCoord_.two_norm() );
      builtinverse_ = builtDetDF_ = true;
    }
  }

  template<> //dim = 1 , dimworld = 3
  inline void ALU3dGridGeometry<0,3, const ALU3dGrid<3,3,tetra> > :: buildJacobianInverseTransposed() const
  {
    if(!builtinverse_)
    {
      detDF_ = 1.0;
      builtinverse_ = builtDetDF_ = true;
    }
  }

  // built Geometry
  template <int mydim, int cdim>
  template <class GeometryType>
  inline bool
  ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, tetra> >::
  buildGeomInFather(const GeometryType &fatherGeom , const GeometryType & myGeom)
  {
    // reset flags, because mappings need to be calculated again
    builtinverse_ = builtA_ = false;

    // compute the local coordinates in father refelem
    for(int i=0; i < myGeom.corners() ; ++i)
    {
      coord_[i] = fatherGeom.local( myGeom.corner( i ) );
      // to avoid rounding errors
      for(int j=0; j<cdim; ++j)
        if ( coord_[i][j] < 1e-16) coord_[i][j] = 0.0;
    }

    calculateDeterminant();
    return true;
  }

  template<int mydim, int cdim>
  template<class coord_t>
  inline void ALU3dGridGeometry<mydim,cdim, const ALU3dGrid<3, 3, tetra> > ::
  copyCoordVec(const coord_t& point,
               FieldVector<alu3d_ctype,cdim> & coord ) const
  {
    assert( cdim == 3 );
    coord[0] = point[0];
    coord[1] = point[1];
    coord[2] = point[2];
  }

  template <>
  inline bool ALU3dGridGeometry<3,3, const ALU3dGrid<3,3,tetra> > ::
  buildGeom(const IMPLElementType & item)
  {
    enum { dim = 3 };
    enum { dimworld = 3};

    builtinverse_ = builtA_ = false;

#ifndef NDEBUG
    builtDetDF_ = true ;

    // if this assertion is thrown, use ElementTopo::dune2aluVertex instead
    // of number when calling myvertex
    for(int i=0; i<4; ++i)
      assert( ElementTopo::dune2aluVertex(i) == i );

#endif
    // get volume and calc elDet
    volume_ = item.volume();
    detDF_  = 6.0 * volume_;

    copyCoordVec(item.myvertex(0)->Point(), coord_[0]);
    copyCoordVec(item.myvertex(1)->Point(), coord_[1]);
    copyCoordVec(item.myvertex(2)->Point(), coord_[2]);
    copyCoordVec(item.myvertex(3)->Point(), coord_[3]);
    return true;
  }

  template <>
  inline bool ALU3dGridGeometry<2,3, const ALU3dGrid<3,3,tetra> > ::
  buildGeom(const HFaceType & item, int twist, int)
  {
    enum { dim = 2 };
    enum { dimworld = 3};

    builtinverse_ = builtA_ = builtDetDF_ = false;

    for (int i=0; i<(dim+1); ++i)
    {
      // Transform Dune index to ALU index and apply twist
      const int rotatedALUIndex = FaceTopo::twist(
        FaceTopo::dune2aluVertex(i), twist);

      // copy coordinates
      copyCoordVec(static_cast<const GEOFaceType
                               &>(item).myvertex(rotatedALUIndex)->Point(), coord_[i]);
    }

    return true;
  }

  template <>
  template <class coord_t>
  inline bool
  ALU3dGridGeometry<2,3, const ALU3dGrid<3, 3, tetra> > ::
  buildGeom(const coord_t& p0,
            const coord_t& p1,
            const coord_t& p2)
  {

    // copy coordinates
    copyCoordVec( p0, coord_[0] );
    copyCoordVec( p1, coord_[1] );
    copyCoordVec( p2, coord_[2] );

    // reset flags
    builtinverse_ = builtA_ = builtDetDF_ = false;
    return true;
  }

  template <>
  inline bool
  ALU3dGridGeometry<2,3, const ALU3dGrid<3, 3, tetra> > ::
  buildGeom(const FaceCoordinatesType& coords)
  {
    return buildGeom( coords[0], coords[1], coords[2] );
  }

  template <> // for edges
  inline bool ALU3dGridGeometry<1,3, const ALU3dGrid<3,3,tetra> > ::
  buildGeom(const HEdgeType & item, int twist, int)
  {
    enum { dim = 1 };
    enum { dimworld = 3};

    builtinverse_ = builtA_ = builtDetDF_ = false;

    for (int i=0; i<(dim+1); ++i)
    {
      // copy coordinates
      copyCoordVec(static_cast<const GEOEdgeType &> (item).myvertex((i+twist)%2)->Point(),
                   coord_[i]);
    }

    buildJacobianInverseTransposed();
    return true;
  }

  template <> // for Vertices ,i.e. Points (note that twist is a dummy parameter here, needed for consistency)
  inline bool ALU3dGridGeometry<0,3, const ALU3dGrid<3,3,tetra> > ::
  buildGeom(const VertexType & item, int, int)
  {
    enum { dim = 0 };
    enum { dimworld = 3};

    builtinverse_ = builtA_ = builtDetDF_ = false;

    const double (&p)[3] = static_cast<const GEOVertexType &> (item).Point();
    for (int j=0; j<dimworld; ++j) coord_[0][j] = p[j];

    buildJacobianInverseTransposed();
    return true;
  }


  /* Comment in for adaptation to new GeometryType */
  template< int mydim, int cdim >
  inline GeometryType
  ALU3dGridGeometry< mydim, cdim, const ALU3dGrid< 3, 3, tetra > > :: type () const
  {
    return myGeomType_;
  }

  template<int mydim, int cdim>
  inline int ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> > ::corners () const
  {
    return corners_;
  }

  template<int mydim, int cdim>
  inline const FieldVector<alu3d_ctype, cdim>&
  ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> > :: operator[] (int i) const
  {
    assert((i>=0) && (i < mydim+1));
    return coord_[i];
  }

  template<int mydim, int cdim>
  inline FieldVector<alu3d_ctype, cdim>
  ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> > :: corner (int i) const
  {
    typedef GenericGeometry::MapNumberingProvider< mydim > Numbering;
    const unsigned int tid = GenericGeometry::topologyId( type() );
    const int j = Numbering::template generic2dune< mydim >( tid, i );
    return this->operator[](j);
  }

  //   G L O B A L   - - -

  // dim = 1,2,3 dimworld = 3
  template<int mydim, int cdim>
  inline FieldVector<alu3d_ctype, cdim>
  ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> >::
  global(const FieldVector<alu3d_ctype, mydim>& local) const
  {
    calcElMatrix();

    globalCoord_ = coord_[0];
    // multiply with transposed because AT is also transposed
    A_.umv(local,globalCoord_);
    return globalCoord_;
  }

  template<> // dim = dimworld = 3
  inline FieldVector<alu3d_ctype, 3>
  ALU3dGridGeometry<3,3,const ALU3dGrid<3,3,tetra> > ::
  local(const FieldVector<alu3d_ctype, 3>& global) const
  {
    if (!builtinverse_) buildJacobianInverseTransposed();

    globalCoord_  = global;
    globalCoord_ -= coord_[0];

    localCoord_ = 0.0;
    // multiply with transposed because Jinv_ is already transposed
    Jinv_.umtv(globalCoord_,localCoord_);
    return localCoord_;
  }

  template<int mydim, int cdim> // dim = dimworld = 3
  inline FieldVector<alu3d_ctype, mydim>
  ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3,3,tetra> > ::
  local(const FieldVector<alu3d_ctype, cdim>& global) const
  {
    if(!builtinverse_)
      buildJacobianInverseTransposed();

    globalCoord_  = global;
    globalCoord_ -= coord_[0];

#if 0
    AT_x_ = 0.0;
    A_.umtv(globalCoord_,AT_x_);

    localCoord_ = 0.0;
    //multipy with transposed, because Jinv_ ist transposed
    Jinv_.umtv(AT_x_,localCoord_);
#endif
    FMatrixHelp :: multAssignTransposed( Jinv_, globalCoord_, localCoord_ );
    return localCoord_;
  }


#if 0
  template<int mydim, int cdim>
  inline bool ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> > ::
  checkInside(const FieldVector<alu3d_ctype, mydim>& local) const
  {
    alu3d_ctype sum = 0.0;

    for(int i=0; i<mydim; i++)
    {
      sum += local[i];
      if(local[i] < 0.0)
      {
        if(std::abs(local[i]) > ALUnumericEpsilon)
        {
          return false;
        }
      }
    }

    if( sum > 1.0 )
    {
      if(sum > (1.0 + ALUnumericEpsilon))
        return false;
    }
    return true;
  }
#endif


  // for elements this is already calculated
  template<>
  inline alu3d_ctype
  ALU3dGridGeometry<3,3,const ALU3dGrid<3, 3, tetra> > ::
  integrationElement (const FieldVector<alu3d_ctype, 3>& local) const
  {
    assert( builtDetDF_ );
    return detDF_;
  }

  template<int mydim, int cdim>
  inline void
  ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> > ::
  calculateDeterminant() const
  {
    calcElMatrix();

    // only can do this for mydim == dim, otherwise not defined
    if(mydim == cdim)
    {
      detDF_ = A_.determinant();
    }
    else
    {
      // buildJacobianInverseTransposed() also calculates determinant
      buildJacobianInverseTransposed();
    }
    assert(detDF_ > 0.0);

    enum { factor = Factorial<mydim>::factorial };
    volume_ = detDF_ / factor ;

    builtDetDF_ = true;
  }

  template<int mydim, int cdim>
  inline alu3d_ctype
  ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> > ::integrationElement (const FieldVector<alu3d_ctype, mydim>& local) const
  {
    if(builtDetDF_)
      return detDF_;

    calculateDeterminant() ;
    return detDF_;
  }

  template<>
  inline alu3d_ctype
  ALU3dGridGeometry<3,3,const ALU3dGrid<3, 3, tetra> > ::volume () const
  {
    assert( builtDetDF_ );
    return volume_;
  }

  template<int mydim, int cdim>
  inline alu3d_ctype
  ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> > ::volume () const
  {
    if( builtDetDF_ )
      return volume_;

    calculateDeterminant() ;
    return volume_;
  }


  //  J A C O B I A N _ I N V E R S E  - - -

  template<> // dim = dimworld = 3
  inline const FieldMatrix<alu3d_ctype,3,3> &
  ALU3dGridGeometry<3,3, const ALU3dGrid<3,3,tetra> >::
  jacobianInverseTransposed (const FieldVector<alu3d_ctype, 3>& local) const
  {
    if (!builtinverse_) buildJacobianInverseTransposed();
    return Jinv_;
  }

  template<> // dim = dimworld = 3
  inline const FieldMatrix< alu3d_ctype, 3, 2 > &
  ALU3dGridGeometry<2,3, const ALU3dGrid<3,3,tetra> >::
  jacobianInverseTransposed (const FieldVector<alu3d_ctype, 2>& local) const
  {
    if (!builtinverse_) buildJacobianInverseTransposed();
    return Jinv_;
  }

  // print the ElementInformation
  template<int mydim, int cdim>
  inline void ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> >::print (std::ostream& ss) const
  {
    ss << "ALU3dGridGeometry<" << mydim << "," << cdim << ", tetra> = {\n";
    for(int i=0; i<corners(); i++)
    {
      ss << " corner " << i << " ";
      ss << "{" << ((*this)[i]) << "}"; ss << std::endl;
    }
    ss << "} \n";
  }

  template<int mydim, int cdim>
  inline bool ALU3dGridGeometry<mydim,cdim,const ALU3dGrid<3, 3, tetra> > ::
  affine() const
  {
    return true;
  }

  /////////////////////////////////////////////////////////////////////////
  //
  //- Hexahedron specialization
  //
  /////////////////////////////////////////////////////////////////////////
  template <int mydim, int cdim>
  inline ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> >::
  ALU3dGridGeometry()
    : geoImpl_(),
      volume_(1.0)
  {}

  template <>
  inline ALU3dGridGeometry<3, 3, const ALU3dGrid<3, 3, hexa> >::
  ALU3dGridGeometry()
    : geoImpl_(),
      volume_(-1.0)
  {}

  template <>
  inline ALU3dGridGeometry<2, 3, const ALU3dGrid<3, 3, hexa> >::
  ALU3dGridGeometry()
    : geoImpl_(),
      volume_(-1.0)
  {}

  template <int mydim, int cdim>
  ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> >::
  ~ALU3dGridGeometry() {}

  template< int mydim, int cdim >
  inline GeometryType
  ALU3dGridGeometry< mydim, cdim, const ALU3dGrid< 3, 3, hexa > > :: type () const
  {
    return GeometryType( GeometryType :: cube, mydim );
  }

  template <int mydim, int cdim>
  inline int
  ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> >::corners() const {
    return corners_;
  }

  template <int mydim, int cdim>
  inline const FieldVector<alu3d_ctype, cdim>&
  ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> >::
  operator[] (int i) const
  {
    return geoImpl_[i];
  }

  template <int mydim, int cdim>
  inline FieldVector<alu3d_ctype, cdim>
  ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> >::
  corner (int i) const
  {
    typedef GenericGeometry::MapNumberingProvider< mydim > Numbering;
    const unsigned int tid = GenericGeometry::topologyId( type() );
    const int j = Numbering::template generic2dune< mydim >( tid, i );
    return this->operator[](j);
  }


  template <>
  inline FieldVector<alu3d_ctype, 3>
  ALU3dGridGeometry<3, 3, const ALU3dGrid<3, 3, hexa> >::
  global (const FieldVector<alu3d_ctype, 3>& local) const
  {
    FieldVector<alu3d_ctype, 3> tmp2_;
    geoImpl_.mapping().map2world(local, tmp2_);
    return tmp2_;
  }

  template <>
  inline FieldVector<alu3d_ctype, 3>
  ALU3dGridGeometry<2, 3, const ALU3dGrid<3, 3, hexa> >::
  global (const FieldVector<alu3d_ctype, 2>& local) const
  {
    FieldVector<alu3d_ctype, 3> tmp2_;
    geoImpl_.mapping().map2world(local, tmp2_);
    return tmp2_;
  }

  template <>
  inline FieldVector<alu3d_ctype, 3>
  ALU3dGridGeometry<3, 3, const ALU3dGrid<3, 3, hexa> >::
  local (const FieldVector<alu3d_ctype, 3>& global) const
  {
    FieldVector<alu3d_ctype, 3> tmp2_;
    geoImpl_.mapping().world2map(global, tmp2_);
    return tmp2_;
  }

  template <>
  inline FieldVector<alu3d_ctype, 2>
  ALU3dGridGeometry<2, 3, const ALU3dGrid<3, 3, hexa> >::
  local (const FieldVector<alu3d_ctype, 3>& global) const
  {
    FieldVector<alu3d_ctype, 2> tmp1_;
    geoImpl_.mapping().world2map(global, tmp1_);
    return tmp1_;
  }


#if 0
  template <int mydim, int cdim>
  inline bool
  ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> >::
  checkInside(const FieldVector<alu3d_ctype, mydim>& local) const
  {
    bool result = true;
    for (int i = 0; i < mydim; i++ ) {
      result &= (local[i] >= - ALUnumericEpsilon &&
                 local[i] <= 1.0 +  ALUnumericEpsilon);
    }
    return result;
  }
#endif


  template<>
  inline alu3d_ctype
  ALU3dGridGeometry<3, 3, const ALU3dGrid<3, 3, hexa> >::
  integrationElement (const FieldVector<alu3d_ctype, 3>& local) const
  {
    return geoImpl_.mapping().det(local);
  }

  template<>
  inline alu3d_ctype
  ALU3dGridGeometry<2, 3, const ALU3dGrid<3, 3, hexa> >::
  integrationElement (const FieldVector<alu3d_ctype, 2>& local) const
  {
    FieldVector<alu3d_ctype, 3> tmp2_;
    geoImpl_.mapping().normal(local, tmp2_);
    return tmp2_.two_norm();
  }

  template<>
  inline alu3d_ctype
  ALU3dGridGeometry<1, 3, const ALU3dGrid<3, 3, hexa> >::
  integrationElement (const FieldVector<alu3d_ctype, 1>& local) const
  {
    return (this->operator[] (0) - this->operator[] (1)).two_norm();
  }

  template<>
  inline alu3d_ctype
  ALU3dGridGeometry<3, 3, const ALU3dGrid<3, 3, hexa> >::
  volume() const
  {
    return volume_;
  }

  template<int mydim, int cdim>
  inline alu3d_ctype
  ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> >::
  volume () const
  {
    return integrationElement(FieldVector<alu3d_ctype, mydim> (0.5));
  }

  template <int mydim, int cdim>
  inline bool
  ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> >::
  affine() const
  {
    return geoImpl_.mapping().affine();
  }

  template <>
  inline const FieldMatrix<alu3d_ctype, 3, 3>&
  ALU3dGridGeometry<3, 3, const ALU3dGrid<3, 3, hexa> >::
  jacobianInverseTransposed (const FieldVector<alu3d_ctype, 3>& local) const
  {
    return geoImpl_.mapping().jacobianInverseTransposed(local);
  }

  template <>
  inline const FieldMatrix<alu3d_ctype, 3, 2>&
  ALU3dGridGeometry<2, 3, const ALU3dGrid<3, 3, hexa> >::
  jacobianInverseTransposed (const FieldVector<alu3d_ctype, 2>& local) const
  {
    return geoImpl_.mapping().jacobianInverseTransposed(local);
  }

  template <int mydim, int cdim>
  inline void
  ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> >::
  print (std::ostream& ss) const {
    ss << "ALU3dGridGeometry<" << mydim << "," << cdim << ", hexa> = {\n";
    for(int i=0; i<corners(); i++)
    {
      ss << " corner " << i << " ";
      ss << "{" << ((*this)[i]) << "}"; ss << std::endl;
    }
    ss << "} \n";
  }

  // built Geometry
  template <int mydim, int cdim>
  template <class GeometryType>
  inline bool
  ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> >::
  buildGeomInFather(const GeometryType &fatherGeom , const GeometryType & myGeom)
  {
    // update geo impl
    geoImpl_.updateInFather( fatherGeom, myGeom );

    // my volume is a part of 1
    volume_ = myGeom.volume() / fatherGeom.volume();

    return true;
  }

  //--hexaBuildGeom
  template <>
  inline bool
  ALU3dGridGeometry<3, 3, const ALU3dGrid<3, 3, hexa> >::
  buildGeom(const IMPLElementType& item)
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

    // get volume of element
    volume_ = item.volume();

    return true;
  }

  template <>
  inline bool
  ALU3dGridGeometry<2,3, const ALU3dGrid<3, 3, hexa> > ::
  buildGeom(const HFaceType & item, int twist, int duneFace )
  {
    // get geo face
    const GEOFaceType& face = static_cast<const GEOFaceType&> (item);

    //assert( duneFace >= 0 && duneFace < 6 );
    if(duneFace < 0 ) duneFace = 0;

    // for all vertices of this face get rotatedIndex
    int rotatedALUIndex[4];
    for (int i = 0; i < 4; ++i)
    {
      // Transform Dune index to ALU index and apply twist
      int localALUIndex = ElementTopo::dune2aluFaceVertex(duneFace,i);
      rotatedALUIndex[i] = FaceTopo::twist(localALUIndex, twist);
    }

    // update geometry implementation
    geoImpl_.update( face.myvertex(rotatedALUIndex[0])->Point(),
                     face.myvertex(rotatedALUIndex[1])->Point(),
                     face.myvertex(rotatedALUIndex[2])->Point(),
                     face.myvertex(rotatedALUIndex[3])->Point() );

    return true;
  }

  template <>
  template <class coord_t>
  inline bool
  ALU3dGridGeometry<2,3, const ALU3dGrid<3, 3, hexa> > ::
  buildGeom(const coord_t& p0,
            const coord_t& p1,
            const coord_t& p2,
            const coord_t& p3)
  {
    // update geometry implementation
    geoImpl_.update( p0, p1, p2, p3 );
    return true;
  }


  template <>
  inline bool
  ALU3dGridGeometry<2,3, const ALU3dGrid<3, 3, hexa> > ::
  buildGeom(const FaceCoordinatesType& coords)
  {
    return buildGeom( coords[0], coords[1], coords[2], coords[3] );
  }

  template <> // for edges
  inline bool
  ALU3dGridGeometry<1,3, const ALU3dGrid<3, 3, hexa> >::
  buildGeom(const HEdgeType & item, int twist, int)
  {
    const GEOEdgeType & edge = static_cast<const GEOEdgeType &> (item);
    // update geometry implementation
    geoImpl_.update( edge.myvertex((twist)  %2)->Point(),
                     edge.myvertex((1+twist)%2)->Point() );
    return true;
  }

  template <> // for Vertices ,i.e. Points
  inline bool
  ALU3dGridGeometry<0,3, const ALU3dGrid<3,3,hexa> >::
  buildGeom(const VertexType & item, int twist, int)
  {
    // update geometry implementation
    geoImpl_.update( static_cast<const GEOVertexType &> (item).Point() );
    return true;
  }


} // end namespace Dune
