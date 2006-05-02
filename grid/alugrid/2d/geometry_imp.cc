// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDGEOMETRYIMP_CC
#define DUNE_ALU2DGRIDGEOMETRYIMP_CC

namespace Dune {

  // implementation of methods

  //**********************************************************************
  //
  // --ALU2dGridGeometry
  // --Geometry
  //**********************************************************************

  // default, do nothing
  template <int mydim, int cdim, class GridImp>
  inline int ALU2dGridGeometry<mydim,cdim,GridImp>::mapVertices (int i) const
  {
    // there is a specialisation for each combination of mydim and coorddim
    //return ALBERTA AlbertHelp :: MapVertices<mydim,cdim>::mapVertices(i,face_,edge_,vertex_);
    return 0;
  }


  template <int mydim, int cdim, class GridImp>
  inline ALU2dGridGeometry<mydim, cdim, GridImp> :: ALU2dGridGeometry() :
    myGeomType_(GeometryType::simplex,mydim)
  {
    // make empty element
    initGeom();
  }

  template <int mydim, int cdim, class GridImp>
  inline ALU2dGridGeometry<mydim,cdim,GridImp>::
  ALU2dGridGeometry(const int child, const int orientation) : myGeomType_(GeometryType::simplex,mydim)
  {
    // make empty element
    buildGeomInFather(child,orientation);
  }

  /*
     template <int mydim, int cdim, class GridImp>
     inline ALBERTA EL_INFO * AlbertaGridGeometry<mydim,cdim,GridImp>::
     makeEmptyElInfo()
     {
     ALBERTA EL_INFO * elInfo = &statElInfo[mydim];

     elInfo->mesh = 0;
     elInfo->el = 0;
     elInfo->parent = 0;
     elInfo->macro_el = 0;
     elInfo->level = 0;
     #if DIM > 2
     elInfo->orientation = 0;
     elInfo->el_type = 0;
     #endif

     for(int i =0; i<mydim+1; i++)
     {
      for(int j =0; j< cdim; j++)
      {
        elInfo->coord[i][j] = 0.0;
        elInfo->opp_coord[i][j] = 0.0;
      }
      elInfo->bound[i] = 0;
     }
     return elInfo;
     }
   */

  template <int mydim, int cdim, class GridImp>
  inline void ALU2dGridGeometry<mydim,cdim,GridImp>::
  initGeom()
  {
    //elInfo_ = 0;
    face_ = 0;
    edge_ = 0;
    vertex_ = 0;
    builtinverse_ = false;
    builtElMat_   = false;
    calcedDet_    = false;
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
  inline const GeometryType & ALU2dGridGeometry<mydim, cdim, GridImp> :: type () const {
    return myGeomType_;
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
    calcedDet_    = false;

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
  inline FieldVector<alu2d_ctype, cdim> ALU2dGridGeometry<mydim, cdim, GridImp> :: global (const FieldVector<alu2d_ctype, mydim>& local) const {
    calcElMatrix();

    globalCoord_ = coord_[0];
    elMat_.umv(local,globalCoord_);
    return globalCoord_;
  }

  //! maps a global coordinate within the element to a
  //! local coordinate in its reference element
  template <int mydim, int cdim, class GridImp>
  inline FieldVector<alu2d_ctype,  mydim> ALU2dGridGeometry <mydim, cdim, GridImp> :: local (const FieldVector<alu2d_ctype, cdim>& global) const {
    if (!builtinverse_)
      buildJacobianInverseTransposed();

    globalCoord_  = global;
    globalCoord_ -= coord_[0];

    AT_x_ = FMatrixHelp::multTransposed(elMat_,globalCoord_);
    localCoord_ = FMatrixHelp::mult(Jinv_,AT_x_);
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

    // calc Jinv_ = (A^T*A)^-1
    std::abs( FMatrixHelp::invertMatrix(elMatT_elMat_,Jinv_) );
    builtinverse_ = true;
    return;
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

    calcedDet_ = true;
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

  template <int mydim, int cdim, class GridImp>
  inline alu2d_ctype ALU2dGridGeometry<mydim,cdim,GridImp>:: volume () const {
    if (cdim == 0) {
      assert( calcedDet_ );
      return 0.5*elDet_;
    }
    else if(cdim == 1) {
      tmpZ_[0] = coord_[0][0] - coord_[1][0];
      tmpZ_[1] = coord_[0][1] - coord_[1][1];
      return tmpZ_.two_norm();
    }
    else
      return 0.;
  }


  template <int mydim, int cdim, class GridImp>
  inline const FieldMatrix<alu2d_ctype,mydim,mydim>& ALU2dGridGeometry<mydim,cdim,GridImp>::
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
        if(std::abs(local[i]) > 1e-15)
        {
          return false;
        }
      }
    }

    if( sum > 1.0 )
    {
      if(sum > (1.0 + 1e-15))
        return false;
    }

    return true;
  }

  /*
     //! built Geometry
     template <int mydim, int cdim, class GridImp>
     inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
     builtGeom(ALU2dGridGeometry::ElementType * item, int face,
            int edge, int vertex)
     {
     //elInfo_ = elInfo;
     face_ = face;
     edge_ = edge;
     vertex_ = vertex;
     builtinverse_ = false;
     builtElMat_   = false;

     if(item)
     {
      // copy coordinates
      for(int i=0; i<mydim+1; ++i)
      {
        // copy coordinates
        for(int j=0; j<cdim; ++j) coord_[i][j] = item->vertex(i)->coord()[j];
      }

      elDet_     = elDeterminant();
      calcedDet_ = true;
      // geometry built
      return true;
     }
     else
     {
      elDet_     = 0.0;
      calcedDet_ = false;
     }
     // geometry not built
     return false;
     }
   */

  //! built Geometry
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  builtGeom(const HElementType & item, int face,
            int edge, int vertex)
  {
    //elInfo_ = elInfo;
    face_ = face;
    edge_ = edge;
    vertex_ = vertex;
    builtinverse_ = false;
    builtElMat_   = false;

    if(&item!=0)
    {
      // copy coordinates
      for(int i=0; i<mydim+1; ++i)
      {
        // copy coordinates
        if (mydim == 1)
          for(int j=0; j<cdim; ++j) coord_[i][j] = item.vertex((i+edge+1)%3)->coord()[j];
        else
          for(int j=0; j<cdim; ++j) coord_[i][j] = item.vertex(i)->coord()[j];
      }

      elDet_     = elDeterminant();
      calcedDet_ = true;
      // geometry built
      return true;
    }
    else
    {
      elDet_     = 0.0;
      calcedDet_ = false;
    }
    // geometry not built
    return false;
  }


  //! built Geometry
  template <int mydim, int cdim, class GridImp>
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  builtGeom(const ALU2DSPACE Vertex & item, int face,
            int edge, int vertex)
  {
    //elInfo_ = elInfo;
    face_ = face;
    edge_ = edge;
    vertex_ = vertex;
    builtinverse_ = false;
    builtElMat_   = false;

    coord_[0][0] = 0;
    coord_[0][1] = 0;
    coord_[1][0] = 0;
    coord_[1][1] = 0;
    coord_[2][0] = 0;
    coord_[2][1] = 0;

    if(&item!=0) {
      coord_[0][0] = item.coord()[0];
      coord_[0][1] = item.coord()[1];

      elDet_     = elDeterminant();
      calcedDet_ = true;
      // geometry built
      return true;
    }
    else
    {
      elDet_     = 0.0;
      calcedDet_ = false;
    }
    // geometry not built
    return false;
  }

  // built Geometry
  template <int mydim, int cdim, class GridImp>
  template <class GeometryType, class LocalGeometryType >
  inline bool ALU2dGridGeometry<mydim,cdim,GridImp>::
  builtLocalGeom(const GeometryType &geo, const LocalGeometryType & localGeom,
                 HElementType * item,int face)
  {
    //elInfo_ = elInfo;
    face_ = face;
    edge_   = 0;
    vertex_ = 0;
    builtinverse_ = false;
    builtElMat_   = false;

    if(item)
    {
      // just map the point of the global intersection to the local
      // coordinates , this is the default procedure
      // for simplices this is not so bad
      for(int i=0; i<mydim+1; i++)
      {
        coord_[i] = geo.local( localGeom[i] );
      }

      elDet_     = elDeterminant();
      calcedDet_ = true;

      // geometry built
      return true;
    }
    else
    {
      elDet_     = 0.0;
      calcedDet_ = false;
    }

    // geometry not built
    return false;
  }

  /*
     // built Geometry
     template <int mydim, int cdim, class GridImp>
     inline void ALU2dGridGeometry<mydim,cdim,GridImp>::
     buildGeomInFather(const int child, const int orientation )
     {
     initGeom();
     // reset coordinate vectors
     coord_ = 0.0;

     assert( (child == 0) || (child == 1) );
     if(mydim == 2)
     {
   */
  /*
     //////////////////////////////////////////////
     //
     //               (0,1)
     //                /|\
     //               /0|1\
     //              /  |  \
     //             /   |   \
     //            /    |    \
     //           /     |     \
     //          /      |      \
     //         / ch 0  | ch 1  \
     //        /1      2|2      0\
     //        -------------------
     //    (0,0)     (0.5,0)    (1,0)
     //
     //
     ///////////////////////////////////////////
   */
  /*
      if( child == 0 )
      {
        coord_[0][1] = 1.0; // (0,1)
        coord_[1]    = 0.0; // (0,0)
        coord_[2][0] = 0.5; // (0.5,0)
      }
      if( child == 1 )
      {
        coord_[0][0] = 1.0; // (1,0)
        coord_[1][1] = 1.0; // (0,1)
        coord_[2][0] = 0.5; // (0.5,0)
      }
     }

     // its a child of the reference element ==> det = 0.5
     elDet_ = elDeterminant();
     calcedDet_ = true;

     if( elDet_ > 0.0 ) return;
     DUNE_THROW(NotImplemented,"wrong dimension given!");
     }
   */

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
      coord_[i] = fatherGeom.local( myGeom[i] );
      // to avoid rounding errors
      for(int j=0; j<cdim; ++j)
        if ( coord_[i][j] < 1e-16) coord_[i][j] = 0.0;
    }

    buildJacobianInverseTransposed();
    return true;
  }

  // altes Zeug!
  /*
     //! generate the geometry for out of given ALU2dGridElement
     template <int mydim, int cdim, class GridImp>
     inline bool ALU2dGridGeometry<mydim, cdim, GridImp> :: buildGeom(const HElementType & item) {
     enum { dim = 2 };
     enum { dimworld = 2};


     builtinverse_ = builtA_ =  false;
     detDF_ = 2.*item.area();
     builtDetDF_ = true;

     for(int i=0; i<3; ++i)
      for(int j=0; j<2;++j)
        coord_[i][j]=item.vertex(i)->coord()[j];
     return true;
     }

     //! generate the geometry for out of given ALU2dGridElement
     template <int mydim, int cdim, class GridImp>
     inline bool ALU2dGridGeometry<mydim, cdim, GridImp> :: buildGeomInFather(const GeometryImp &fatherGeom , const GeometryImp & myGeom)
     {
     // reset flags, because mappings need to be calculated again
     builtinverse_ = builtA_ = builtDetDF_ = false;

     // compute the local coordinates in father refelem
     for(int i=0; i < myGeom.corners() ; i++)
      for(int j=0; j<2;++j)
        coord_[i] = fatherGeom.local( myGeom[i] );
     return true;
     }

     //! generate transposed Jacobian Inverse and calculate integration_element
     template <int mydim, int cdim, class GridImp>
     inline void ALU2dGridGeometry <mydim, cdim, GridImp> :: buildJacobianInverseTransposed() const {
     if(!builtinverse_)  {
      calcElMatrix();

      // DetDf = integrationElement
      detDF_ = std::abs( FMatrixHelp::invertMatrix(AT_,Jinv_) );
      // transpose Jinv_
      for (int i = 0; i < matdim; ++i) {
        for (int j = i+1; j < matdim; ++j) {
          alu2d_ctype tmp = Jinv_[i][j];
          Jinv_[i][j] = Jinv_[j][i];
          Jinv_[j][i] = tmp;
        }
      }
      //builtinverse_ = builtDetDF_ = true;
      builtinverse_  = true;
     }
     }

     //! calculates the element matrix for calculation of the jacobian inverse
     template <int mydim, int cdim, class GridImp>
     inline void ALU2dGridGeometry <mydim, cdim, GridImp> :: calcElMatrix () const {
     if(!builtA_)
     {
      // create Matrix A (=Df)               INDIZES: row/col
      // Mapping: R^dim -> R^2,  F(x) = A x + p_0
      // columns:    p_1 - p_0  |  p_2 - p_0

      for (int i=0; i<mydim; i++) {
        AT_[i] = coord_[i+1] - coord_[0];
      }
      builtA_ = true;
     }
     }

     //! A(l) , see grid.hh
     template <int mydim, int cdim, class GridImp>
     inline alu2d_ctype ALU2dGridGeometry <mydim, cdim, GridImp> :: integrationElement (const FieldVector<alu2d_ctype, mydim>& local) const {
     // if (!builtDetDF_)
     //   buildJacobianInverseTransposed();
     return detDF_;
     }

     //! can only be called for dim=dimworld!
     template <int mydim, int cdim, class GridImp>
     inline const FieldMatrix<alu2d_ctype,mydim,mydim>& ALU2dGridGeometry <mydim, cdim, GridImp> ::
     jacobianInverseTransposed (const FieldVector<alu2d_ctype, cdim>& local) const {
     if(!builtinverse_)
      buildJacobianInverseTransposed();
     return Jinv_;
     }

   */

} //end namespace Dune

#endif
