// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDGEOMETRYIMP_CC
#define DUNE_ALU2DGRIDGEOMETRYIMP_CC

#include <dune/grid/genericgeometry/conversion.hh>

namespace Dune {

  // implementation of methods

  //**********************************************************************
  //
  // --ALU2dGridGeometry
  // --Geometry
  //**********************************************************************


  template <int mydim, int cdim, class GridImp>
  inline ALU2dGridGeometry<mydim, cdim, GridImp> :: ALU2dGridGeometry()
    : coord_(0.0)
      , refCoord_(0.0)
      , builtElMat_(false)
      , builtinverse_(false)
#ifndef NDEBUG
      , calcedDet_(false)
#endif
      , up2Date_(false)
  {
    setReferenceCoordinates();
  }

  template <int mydim, int cdim, class GridImp>
  inline ALU2dGridGeometry<mydim,cdim,GridImp>::
  ALU2dGridGeometry(const int child, const int orientation)
    : coord_(0.0)
      , refCoord_(0.0)
      , builtElMat_(false)
      , builtinverse_(false)
#ifndef NDEBUG
      , calcedDet_(false)
#endif
      , up2Date_(false)
  {
    setReferenceCoordinates();

    // make empty element
    buildGeomInFather(child,orientation);
  }

  template <int mydim, int cdim, class GridImp>
  inline void ALU2dGridGeometry<mydim, cdim, GridImp> :: setReferenceCoordinates()
  {
    // zero all entries
    refCoord_ = 0.0;

    // point 1
    refCoord_[1][0] = 1;
    // point 2
    refCoord_[2][1] = 1;

    // length of faces
    refCoord_[0][2] = M_SQRT2;
    refCoord_[1][2] = 1.0;
    refCoord_[2][2] = 1.0;
  }

  //! print the GeometryInformation
  template <int mydim, int cdim, class GridImp>
  inline void ALU2dGridGeometry <mydim, cdim, GridImp> :: print (std::ostream& ss) const {

    ss << "ALU2dGridGeometry<" << mydim << "," << cdim << "> = {\n";
    for(int i=0; i<corners(); i++)
    {
      ss << " corner " << i << " ";
      ss << "{" << ((*this)[i]) << "}"; ss << std::endl;
    }
    ss << "} \n";
  }

  //! return the element type identifier
  //! line , triangle or tetrahedron, depends on dim
  template <int mydim, int cdim, class GridImp>
  inline const GeometryType ALU2dGridGeometry< mydim, cdim, GridImp >::type () const
  {
    return GeometryType( GeometryType::simplex, mydim );
  }

  //! return the number of corners of this element. Corners are numbered 0...mydim+1
  template <int mydim, int cdim, class GridImp>
  inline int ALU2dGridGeometry<mydim, cdim, GridImp> :: corners () const {
    return mydim+1;
  }

  ///////////////////////////////////////////////////////////////////////

  //! access to coordinates of corners. Index is the number of the corner
  template <int mydim, int cdim, class GridImp>
  inline const FieldVector<alu2d_ctype, cdim>& ALU2dGridGeometry<mydim, cdim, GridImp> :: operator[] (int i) const {
    assert((i >= 0) && (i < corners()));
    return coord_[i];
  }

  //! access to coordinates of corners. Index is the number of the corner
  template <int mydim, int cdim, class GridImp>
  inline FieldVector<alu2d_ctype, cdim> ALU2dGridGeometry<mydim, cdim, GridImp> :: corner(int i) const {
    typedef GenericGeometry::MapNumberingProvider< mydim > Numbering;
    const unsigned int tid = GenericGeometry::topologyId( type() );
    const int j = Numbering::template generic2dune< mydim >( tid, i );
    return this->operator[](j);
  }


  ///////////////////////////////////////////////////////////////////////

  template <int mydim, int cdim, class GridImp>
  inline FieldVector<alu2d_ctype, cdim>& ALU2dGridGeometry<mydim,cdim,GridImp>::
  getCoordVec (int i)
  {
    assert( i >= 0 );
    assert( i < mydim+1 );

    // if global, or jacobianInverse is called then
    // matrix has to be calculated again , because coord might have changed
    builtinverse_ = false;
    builtElMat_   = false;

#ifndef NDEBUG
    calcedDet_    = false;
#endif

    return coord_[i];
  }

  template <class GridImp, int mydim, int cdim>
  struct CalcElementMatrix
  {
    enum { matdim = (mydim > 0) ? mydim : 1 };
    static bool calcElMatrix(const FieldMatrix<alu2d_ctype,mydim+1,cdim> & coord,
                             FieldMatrix<alu2d_ctype,matdim,matdim> & elMat)
    {
      std::string text;
      text += "ALU2dGridGeometry<";
      char fake[128];
      sprintf(fake,"%d",mydim);
      text += fake; text += ",";
      sprintf(fake,"%d",cdim); text += ">::calcElMatrix: No default implementation!";
      //DUNE_THROW(AlbertaError, text);
      return false;
    }
  };


  template <class GridImp>
  struct CalcElementMatrix<GridImp,1,2>
  {
    enum { mydim  = 1 };
    enum { cdim   = 2 };
    static bool calcElMatrix(const FieldMatrix<alu2d_ctype,mydim+1,cdim> & coord,
                             FieldMatrix<alu2d_ctype,cdim,mydim> & elMat)
    {
      //      column 0
      // A = ( P1 - P0 )
      for (int i=0; i<cdim; ++i)
      {
        elMat[i][0] = coord[1][i] - coord[0][i];
      }
      return true;
    }
  };

  template <class GridImp>
  struct CalcElementMatrix<GridImp,2,2>
  {
    enum { mydim  = 2 };
    enum { cdim   = 2 };
    enum { matdim = 2 };
    static bool calcElMatrix(const FieldMatrix<alu2d_ctype,mydim+1,cdim> & coord,
                             FieldMatrix<alu2d_ctype,matdim,matdim> & elMat)
    {
      //       column 0 , column 1
      // A = ( P1 - P0  , P2 - P0 )
      for (int i=0; i<cdim; ++i)
      {
        elMat[i][0] = coord[1][i] - coord[0][i];
        elMat[i][1] = coord[2][i] - coord[0][i];
      }
      return true;
    }
  };

  template <int mydim, int cdim, class GridImp>
  inline void ALU2dGridGeometry<mydim,cdim,GridImp>::calcElMatrix () const
  {
    if(!builtElMat_)
    {
      // build mapping from reference element to actual element
      builtElMat_ = CalcElementMatrix<GridImp,mydim,cdim>::calcElMatrix(coord_,elMat_);
    }
  }

  //! maps a local coordinate within reference element to
  //! global coordinate in element
  template <int mydim, int cdim, class GridImp>
  inline FieldVector<alu2d_ctype, cdim> ALU2dGridGeometry<mydim, cdim, GridImp> ::
  global (const FieldVector<alu2d_ctype, mydim>& local) const
  {
    calcElMatrix();

    globalCoord_ = coord_[0];
    elMat_.umv(local,globalCoord_);
    return globalCoord_;
  }

  //! maps a local coordinate within reference element to
  //! global coordinate in element
  template <>
  inline FieldVector<alu2d_ctype, 2> ALU2dGridGeometry<0, 2, const ALU2dGrid<2,2> > ::
  global (const FieldVector<alu2d_ctype, 0>& local) const
  {
    return coord_[0];
  }

  //! maps a global coordinate within the element to a
  //! local coordinate in its reference element
  template <int mydim, int cdim, class GridImp>
  inline FieldVector<alu2d_ctype,  mydim> ALU2dGridGeometry <mydim, cdim, GridImp> ::
  local (const FieldVector<alu2d_ctype, cdim>& global) const
  {
    if (!builtinverse_)
      buildJacobianInverseTransposed();

    globalCoord_  = global;
    globalCoord_ -= coord_[0];

    //AT_x_ = FMatrixHelp::multTransposed(elMat_,globalCoord_);
    //localCoord_ = FMatrixHelp::mult(Jinv_,AT_x_);
    FMatrixHelp :: multAssignTransposed( Jinv_, globalCoord_, localCoord_ );
    return localCoord_;
  }

  template <>
  inline FieldVector<alu2d_ctype, 2> ALU2dGridGeometry<2,2,const ALU2dGrid<2,2> >::
  local(const FieldVector<alu2d_ctype, 2>& global) const
  {
    if(!builtinverse_)
      buildJacobianInverseTransposed();

    globalCoord_  = global;
    globalCoord_ -= coord_[0];
    FMatrixHelp::multAssignTransposed(Jinv_,globalCoord_,localCoord_);

    return localCoord_;
  }

  template <>
  inline FieldVector<alu2d_ctype, 0> ALU2dGridGeometry<0,2,const ALU2dGrid<2,2> >::
  local(const FieldVector<alu2d_ctype, 2>& global) const
  {
    return FieldVector<alu2d_ctype, 0> (1);
  }

  //! determinant of one Geometry, here line
  template <>
  inline alu2d_ctype ALU2dGridGeometry<1,2,const ALU2dGrid<2,2> >::elDeterminant () const
  {
    // volume is length of edge
    tmpZ_ = coord_[0] - coord_[1];
    return std::abs ( tmpZ_.two_norm() );
  }

  //! determinant of one Geometry, here triangle
  template <>
  inline alu2d_ctype ALU2dGridGeometry<2,2,const ALU2dGrid<2,2> >::elDeterminant () const
  {
    calcElMatrix();
    return std::abs ( elMat_.determinant () );
  }

  //! volume of one Geometry, here point
  template <>
  inline alu2d_ctype ALU2dGridGeometry<0,2,const ALU2dGrid<2,2> >::elDeterminant () const
  {
    return 1.0;
  }

  template <>
  inline void ALU2dGridGeometry<1,2,const ALU2dGrid<2,2> >::
  buildJacobianInverseTransposed() const
  {
    // calc A and stores it in elMat_
    calcElMatrix();
    assert( builtElMat_ == true );

    // calc ret = A^T*A
    FMatrixHelp::multTransposedMatrix(elMat_,elMatT_elMat_);

    // calc Jinv_ = A (A^T*A)^-1
    FieldMatrix< alu2d_ctype, matdim, matdim > inv_elMatT_elMat;
    FMatrixHelp::invertMatrix(elMatT_elMat_, inv_elMatT_elMat);
    FMatrixHelp :: multMatrix( elMat_, inv_elMatT_elMat, Jinv_ );
    //std::abs( FMatrixHelp::invertMatrix(elMatT_elMat_,Jinv_) );
    builtinverse_ = true;
  }

  // this method is for (dim==dimworld) = 2
  template <int mydim, int cdim, class GridImp>
  inline void ALU2dGridGeometry<mydim,cdim,GridImp>::
  buildJacobianInverseTransposed() const
  {
    //******************************************************
    //
    //  the mapping is:
    //
    //  F(T) = D where T is the reference element
    //  and D the actual element
    //
    //  F(x) = A * x + b    with   A := ( P_0 , P_1 )
    //
    //  A consist of the column vectors P_0 and P_1 and
    //  is calculated by the method calcElMatrix
    //
    //******************************************************

    // calc A and stores it in elMat_
    calcElMatrix();

    // Jinv = A^-1^T
    assert( builtElMat_ == true );
    // here the transposed jacobian inverse is calculated
    elDet_ = std::abs( FMatrixHelp::invertMatrix_retTransposed(elMat_,Jinv_) );

    assert(elDet_ > 1.0E-25 );

#ifndef NDEBUG
    calcedDet_ = true;
#endif

    builtinverse_ = true;
    return;
  }

  template <int mydim, int cdim, class GridImp>
  inline alu2d_ctype ALU2dGridGeometry<mydim,cdim,GridImp>::
  integrationElement (const FieldVector<alu2d_ctype, mydim>& local) const
  {
    assert( calcedDet_ );
    return elDet_;
  }

  template <>
  inline alu2d_ctype ALU2dGridGeometry<2,2, const ALU2dGrid<2,2> >:: volume () const
  {
    assert( calcedDet_ );
    return 0.5*elDet_;
  }

  template <>
  inline alu2d_ctype ALU2dGridGeometry<1,2, const ALU2dGrid<2,2> >:: volume () const
  {
    assert( calcedDet_ );
    return elDet_;
  }


  template <int mydim, int cdim, class GridImp>
  inline alu2d_ctype ALU2dGridGeometry<mydim,cdim,GridImp>:: volume () const
  {
    assert( mydim == 0 );
    return 1.;
  }

  template <int mydim, int cdim, class GridImp>
  inline const FieldMatrix<alu2d_ctype,cdim,mydim>& ALU2dGridGeometry<mydim,cdim,GridImp>::
  jacobianInverseTransposed (const FieldVector<alu2d_ctype, mydim>& local) const
  {
    if(builtinverse_)
      return Jinv_;

    // builds the jacobian inverse and calculates the volume
    buildJacobianInverseTransposed();
    return Jinv_;
  }

  //! returns true if the point in local coordinates is inside reference element
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry <mydim, cdim, GridImp> :: checkInside(const FieldVector<alu2d_ctype, mydim>& local) const {
    alu2d_ctype sum = 0.0;

    for(int i=0; i<mydim; i++)
    {
      sum += local[i];
      if(local[i] < 0.0)
      {
        if(std::abs(local[i]) > 1e-8)
        {
          return false;
        }
      }
    }

    if( sum > 1.0 )
    {
      if(sum > (1.0 + 1e-8))
        return false;
    }

    return true;
  }

  //! built Geometry
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  builtGeom(const HElementType & item, int face)
  {
    builtinverse_ = false;
    builtElMat_   = false;

    assert( &item );

    // defined in geometry.hh
    elDet_ = CopyCoordinates<GeometryImp,mydim>::copy(item,coord_,face);
#ifndef NDEBUG
    calcedDet_ = true;
#endif

    // geom is up2date
    up2Date_ = true;

    // geometry built
    return true;
  }


  //! built Geometry
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  builtGeom(const ALU2DSPACE Vertex & item, int )
  {
    builtinverse_ = false;
    builtElMat_   = false;

    // set coordinates to zero
    coord_ = 0.0;

    if(&item!=0)
    {
      coord_[0][0] = item.coord()[0];
      coord_[0][1] = item.coord()[1];

      elDet_     = 1.0; // inant();
#ifndef NDEBUG
      calcedDet_ = true;
#endif
      // geom is up2date
      up2Date_ = true;

      // geometry built
      return true;
    }

    // default values
    elDet_     = 0.0;
#ifndef NDEBUG
    calcedDet_ = false;
#endif

    // geom is not up2date
    up2Date_ = false;

    // geometry not built
    return false;
  }

  // built Geometry
  template <int mydim, int cdim, class GridImp>
  template <class GeometryType, class LocalGeometryType >
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  builtLocalGeom(const GeometryType &geo, const LocalGeometryType & localGeom)
  {
    builtinverse_ = false;
    builtElMat_   = false;

    // just map the point of the global intersection to the local
    // coordinates , this is the default procedure
    // for simplices this is not so bad
    for( int i = 0; i <= mydim; ++i )
      coord_[ i ] = geo.local( localGeom.corner( i ) );

    elDet_  = elDeterminant();
#ifndef NDEBUG
    calcedDet_ = true;
#endif

    // geom is up2date
    up2Date_ = true;

    // geometry built
    return true;
  }

  // built Geometry (faceNumber is in generic numbering)
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  builtLocalGeom(const int faceNumber, const int twist)
  {
    builtinverse_ = false;
    builtElMat_   = false;

    assert( twist == 0 || twist == 1 );
    assert( mydim == 1 );

    // just map the point of the global intersection to the local
    // coordinates , this is the default procedure
    // for simplices this is not so bad

    CopyCoordinates<GeometryImp,mydim>::
    copyData( &refCoord_[(2 - faceNumber + (twist%2) + 1)%3][0] , coord_[0]);

    CopyCoordinates<GeometryImp,mydim>::
    copyData( &refCoord_[(2 - faceNumber + ((twist+1)%2) + 1)%3][0] , coord_[1]);

    // get length of faces
    elDet_ = refCoord_[2 - faceNumber][2];
    assert( std::abs(elDet_ - elDeterminant()) < 1e-10 );

#ifndef NDEBUG
    calcedDet_ = true;
#endif

    // geom is up2date
    up2Date_ = true;

    // geometry built
    return true;
  }

  // built Geometry
  template <int mydim, int cdim, class GridImp >
  inline bool ALU2dGridGeometry<mydim, cdim,GridImp>::
  buildGeomInFather(const Geometry & fatherGeom , const Geometry & myGeom)
  {
    // reset flags, because mappings need to be calculated again
    builtinverse_ = builtElMat_ = false;

    // compute the local coordinates in father refelem
    for(int i=0; i < myGeom.corners() ; ++i)
    {
      coord_[i] = fatherGeom.local( myGeom.corner( i ) );
      // to avoid rounding errors
      for(int j=0; j<cdim; ++j)
        if ( coord_[i][j] < 1e-16) coord_[i][j] = 0.0;
    }

    buildJacobianInverseTransposed();

    // geom is up2date
    up2Date_ = true;

    return true;
  }

} //end namespace Dune

#endif
