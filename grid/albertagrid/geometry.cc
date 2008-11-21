// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GEOMETRY_CC
#define DUNE_ALBERTA_GEOMETRY_CC

#include <dune/grid/albertagrid/geometry.hh>

namespace Dune
{

  static ALBERTA EL_INFO statElInfo[DIM+1];

  // AlbertaGridGeometry
  // -------------------

  template <int mydim, int cdim, class GridImp>
  inline int AlbertaGridGeometry<mydim,cdim,GridImp>::mapVertices (int i) const
  {
    // there is a specialisation for each combination of mydim and coorddim
    return ALBERTA AlbertHelp :: MapVertices<mydim,cdim>::mapVertices(i,face_,edge_,vertex_);
  }

  template <int mydim, int cdim, class GridImp>
  inline AlbertaGridGeometry<mydim,cdim,GridImp>::
  AlbertaGridGeometry()
    : detFactor_ ( (albertCtype) 1.0/ Factorial<mydim> :: factorial )
  {
    // make empty element
    initGeom();
  }

  template <int mydim, int cdim, class GridImp>
  inline AlbertaGridGeometry<mydim,cdim,GridImp>::
  AlbertaGridGeometry(const int child, const int orientation)
    : detFactor_ ( (albertCtype) 1.0/ Factorial<mydim> :: factorial )
  {
    // make empty element
    buildGeomInFather(child,orientation);
  }

  template <int mydim, int cdim, class GridImp>
  inline void AlbertaGridGeometry<mydim,cdim,GridImp>::
  initGeom()
  {
    face_ = 0;
    edge_ = 0;
    vertex_ = 0;
    builtinverse_ = false;
    builtElMat_   = false;
    calcedDet_    = false;
  }

  // print the GeometryInformation
  template <int mydim, int cdim, class GridImp>
  inline void AlbertaGridGeometry<mydim,cdim,GridImp>::print (std::ostream& ss) const
  {
    ss << "AlbertaGridGeometry<" << mydim << "," << cdim << "> = { \n";
    for(int i=0; i<corners(); i++)
    {
      ss << " corner " << i << " = ";
      ss << ((*this)[i]); ss << "\n";
    }
    ss << "} \n";
  }

  template< int mydim, int cdim, class GridImp >
  inline GeometryType
  AlbertaGridGeometry< mydim, cdim, GridImp > :: type () const
  {
    return GeometryType( GeometryType :: simplex, mydim );
  }

  template <int mydim, int cdim, class GridImp>
  inline int AlbertaGridGeometry<mydim,cdim,GridImp>::corners() const
  {
    return (mydim+1);
  }

  ///////////////////////////////////////////////////////////////////////
  template <int mydim, int cdim, class GridImp>
  inline const FieldVector<albertCtype, cdim>& AlbertaGridGeometry<mydim,cdim,GridImp>::
  operator [](int i) const
  {
    return coord_[i];
  }

  ///////////////////////////////////////////////////////////////////////
  template <int mydim, int cdim, class GridImp>
  inline FieldVector<albertCtype, cdim>& AlbertaGridGeometry<mydim,cdim,GridImp>::
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
  struct AlbertaCalcElementMatrix
  {
    enum { matdim = (mydim > 0) ? mydim : 1 };
    static bool calcElMatrix(const FieldMatrix<albertCtype,mydim+1,cdim> & coord,
                             FieldMatrix<albertCtype,matdim,matdim> & elMat)
    {
      std::string text;
      text += "AlbertaGridGeometry<";
      char fake[128];
      sprintf(fake,"%d",mydim);
      text += fake; text += ",";
      sprintf(fake,"%d",cdim); text += ">::calcElMatrix: No default implementation!";
      DUNE_THROW(AlbertaError, text);
      return false;
    }
  };

  template <class GridImp>
  struct AlbertaCalcElementMatrix<GridImp,1,2>
  {
    enum { mydim  = 1 };
    enum { cdim   = 2 };
    static bool calcElMatrix(const FieldMatrix<albertCtype,mydim+1,cdim> & coord,
                             FieldMatrix<albertCtype,cdim,mydim> & elMat)
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
  struct AlbertaCalcElementMatrix<GridImp,2,2>
  {
    enum { mydim  = 2 };
    enum { cdim   = 2 };
    enum { matdim = 2 };
    static bool calcElMatrix(const FieldMatrix<albertCtype,mydim+1,cdim> & coord,
                             FieldMatrix<albertCtype,matdim,matdim> & elMat)
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

  template <class GridImp>
  struct AlbertaCalcElementMatrix<GridImp,2,3>
  {
    enum { mydim  = 2 };
    enum { cdim   = 3 };
    static bool calcElMatrix(const FieldMatrix<albertCtype,mydim+1,cdim> & coord,
                             FieldMatrix<albertCtype,cdim,mydim> & elMat)
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

  template <class GridImp>
  struct AlbertaCalcElementMatrix<GridImp,3,3>
  {
    enum { mydim  = 3 };
    enum { cdim   = 3 };
    enum { matdim = 3 };
    static bool calcElMatrix(const FieldMatrix<albertCtype,mydim+1,cdim> & coord,
                             FieldMatrix<albertCtype,matdim,matdim> & elMat)
    {
      const FieldVector<albertCtype, cdim> & coord0 = coord[0];
      for(int i=0 ; i<cdim; ++i)
      {
        elMat[i][0] = coord[1][i] - coord0[i];
        elMat[i][1] = coord[2][i] - coord0[i];
        elMat[i][2] = coord[3][i] - coord0[i];
      }
      return true;
    }
  };

  template <int mydim, int cdim, class GridImp>
  inline void AlbertaGridGeometry<mydim,cdim,GridImp>::calcElMatrix () const
  {
    if(!builtElMat_)
    {
      // build mapping from reference element to actual element
      builtElMat_ = AlbertaCalcElementMatrix<GridImp,mydim,cdim>::calcElMatrix(coord_,elMat_);
    }
  }

  template< int mydim, int cdim, class GridImp >
  inline FieldVector< albertCtype, cdim >
  AlbertaGridGeometry< mydim, cdim, GridImp>
  :: global ( const FieldVector< albertCtype, mydim > &local ) const
  {
    calcElMatrix();

    FieldVector< albertCtype, cdim > y = coord_[ 0 ];
    elMat_.umv( local, y );
    return y;
  }

  //local implementation for mydim < cdim
  template< int mydim, int cdim, class GridImp >
  inline FieldVector< albertCtype, mydim >
  AlbertaGridGeometry<mydim,cdim,GridImp>
  :: local ( const FieldVector< albertCtype, cdim > &global ) const
  {
    if( !builtinverse_ )
      buildJacobianInverseTransposed();

    FieldVector< albertCtype, cdim > y = global;
    y -= coord_[ 0 ];

    FieldVector< albertCtype, mydim > x;
    FMatrixHelp::multAssignTransposed( Jinv_, y, x );
    return x;
  }

  // determinant of one Geometry, here line
  template <>
  inline albertCtype AlbertaGridGeometry<1,2,const AlbertaGrid<2,2> >::elDeterminant () const
  {
    // volume is length of edge
    tmpZ_ = coord_[0] - coord_[1];
    return std::abs ( tmpZ_.two_norm() );
  }

  // determinant of one Geometry, here line
  template <>
  inline albertCtype AlbertaGridGeometry<1,3,const AlbertaGrid<3,3> >::elDeterminant () const
  {
    // volume is length of edge
    tmpZ_ = coord_[0] - coord_[1];
    return std::abs ( tmpZ_.two_norm() );
  }

  // determinant of one Geometry, here triangle
  template <>
  inline albertCtype AlbertaGridGeometry<2,2,const AlbertaGrid<2,2> >::elDeterminant () const
  {
    calcElMatrix();
    return std::abs ( elMat_.determinant () );
  }

  // determinant of one Geometry, here triangle in 3d
  template <>
  inline albertCtype AlbertaGridGeometry<2,3,const AlbertaGrid<3,3> >::elDeterminant () const
  {
    enum { dim = 3 };

    // create vectors of face
    tmpV_ = coord_[1] - coord_[0];
    tmpU_ = coord_[2] - coord_[1];

    // calculate scaled outer normal
    for(int i=0; i<dim; i++)
    {
      tmpZ_[i] = (  tmpU_[(i+1)%dim] * tmpV_[(i+2)%dim]
                    - tmpU_[(i+2)%dim] * tmpV_[(i+1)%dim] );
    }

    // tmpZ is the same as 2.0 * the outer normal
    return std::abs( tmpZ_.two_norm() );
  }

  // volume of one Geometry, here therahedron
  template <>
  inline albertCtype AlbertaGridGeometry<3,3,const AlbertaGrid<3,3> >::elDeterminant () const
  {
    calcElMatrix();
    return std::abs(elMat_.determinant ());
  }

  // volume of one Geometry, here point
  template <>
  inline albertCtype AlbertaGridGeometry<0,2,const AlbertaGrid<2,2> >::elDeterminant () const
  {
    return 1.0;
  }
  // volume of one Geometry, here point
  template <>
  inline albertCtype AlbertaGridGeometry<0,3,const AlbertaGrid<3,3> >::elDeterminant () const
  {
    return 1.0;
  }

  // this method is for (mydim < cdim)
  template< int mydim, int cdim, class GridImp >
  inline void AlbertaGridGeometry< mydim, cdim, GridImp >
  :: buildJacobianInverseTransposed () const
  {
    // calc A and stores it in elMat_
    calcElMatrix();
    assert( builtElMat_ == true );

    // calc ret = A^T*A
    FMatrixHelp::multTransposedMatrix(elMat_,elMatT_elMat_);

    // calc Jinv_ = A (A^T*A)^-1
    FieldMatrix< albertCtype, mydim, mydim > inv_elMatT_elMat;
    FMatrixHelp :: invertMatrix( elMatT_elMat_, inv_elMatT_elMat );
    FMatrixHelp :: multMatrix( elMat_, inv_elMatT_elMat, Jinv_ );

    builtinverse_ = true;
  }

  // this method is for (mydim = cdim = 2)
  template<>
  inline void AlbertaGridGeometry< 2, 2, const AlbertaGrid< 2, 2 > >
  :: buildJacobianInverseTransposed () const
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
  }

  // this method is for (mydim = cdim = 2)
  template<>
  inline void AlbertaGridGeometry< 3, 3, const AlbertaGrid< 3, 3 > >
  :: buildJacobianInverseTransposed () const
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
  }

  template <int mydim, int cdim, class GridImp>
  inline albertCtype AlbertaGridGeometry<mydim,cdim,GridImp>::
  volume () const
  {
    assert( calcedDet_ );
    return detFactor_ * elDet_;
  }

  template <int mydim, int cdim, class GridImp>
  inline albertCtype AlbertaGridGeometry<mydim,cdim,GridImp>::
  integrationElement (const FieldVector<albertCtype, mydim>& local) const
  {
    assert( calcedDet_ );
    return elDet_;
  }

  template <int mydim, int cdim, class GridImp>
  inline const FieldMatrix<albertCtype,cdim,mydim>& AlbertaGridGeometry<mydim,cdim,GridImp>::
  jacobianInverseTransposed (const FieldVector<albertCtype, mydim>& local) const
  {
    if(builtinverse_)
      return Jinv_;

    // builds the jacobian inverse and calculates the volume
    buildJacobianInverseTransposed();
    return Jinv_;
  }

  template <int mydim, int cdim, class GridImp>
  inline bool AlbertaGridGeometry<mydim,cdim,GridImp>::
  checkInside(const FieldVector<albertCtype, mydim> &local) const
  {
    albertCtype sum = 0.0;

    for(int i=0; i<mydim; i++)
    {
      sum += local[i];
      if(local[i] < 0.0)
      {
        if(std::abs(local[i]) > 1e-13)
        {
          return false;
        }
      }
    }

    if( sum > 1.0 )
    {
      if(sum > (1.0 + 1e-13))
        return false;
    }

    return true;
  }

  // built Geometry
  template <int mydim, int cdim, class GridImp>
  inline bool AlbertaGridGeometry<mydim,cdim,GridImp>::
  builtGeom(const GridImp & grid, ALBERTA EL_INFO *elInfo, int face,
            int edge, int vertex)
  {
    face_ = face;
    edge_ = edge;
    vertex_ = vertex;
    builtinverse_ = false;
    builtElMat_   = false;

    if( elInfo != NULL )
    {
      // copy coordinates
      for(int i=0; i<mydim+1; ++i)
      {
        const ALBERTA REAL_D &elcoord = grid.getCoord( elInfo, mapVertices( i ) );
        // copy coordinates
        for(int j=0; j<cdim; ++j) coord_[i][j] = elcoord[j];
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


  // specialization yields speed up, because vertex_ .. is not copied
  template <>
  inline bool AlbertaGridGeometry<2,2,const AlbertaGrid<2,2> >::
  builtGeom(const GridType & grid, ALBERTA EL_INFO *elInfo, int face,
            int edge, int vertex)
  {
    typedef AlbertaGrid<2,2> GridImp;
    enum { dim = 2 };
    enum { dimworld = 2 };

    builtinverse_ = false;
    builtElMat_   = false;

    if( elInfo != NULL )
    {
      // copy coordinates
      for(int i=0; i<dim+1; ++i)
      {
        const ALBERTA REAL_D &elcoord = grid.getCoord( elInfo, mapVertices( i ) );
        for(int j=0; j<dimworld; ++j)
          coord_[i][j] = elcoord[j];
      }

      const ALBERTA EL *el = elInfo->el;
      assert( el );
      // if leaf element, get determinant from leaf data
      if( IS_LEAF_EL(el) )
      {
        typedef GridImp :: LeafDataType::Data LeafDataType;
        LeafDataType * ldata = (LeafDataType *) el->child[1];
        assert( ldata );
        elDet_ = ldata->determinant;
        assert( std::abs( elDet_ ) > 0.0 );
        calcedDet_ = true;
      }
      else
      {
        elDet_     = elDeterminant();
        calcedDet_ = true;
      }

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

  template <>
  inline bool AlbertaGridGeometry<3,3,const AlbertaGrid<3,3> >::
  builtGeom(const GridType & grid, ALBERTA EL_INFO *elInfo, int face,
            int edge, int vertex)
  {
    typedef AlbertaGrid<3,3> GridImp;
    enum { dim = 3 };
    enum { dimworld = 3 };

    builtinverse_ = false;
    builtElMat_   = false;

    if( elInfo != NULL )
    {
      // copy coordinates
      for(int i=0; i<dim+1; ++i)
      {
        const ALBERTA REAL_D &elcoord = grid.getCoord( elInfo, mapVertices( i ) );
        for(int j=0; j<dimworld; ++j)
          coord_[i][j] = elcoord[j];
      }

      const ALBERTA EL *el = elInfo->el;
      assert( el );
      // if leaf element, get determinant from leaf data
      if( IS_LEAF_EL(el) )
      {
        typedef GridImp :: LeafDataType::Data LeafDataType;
        LeafDataType * ldata = (LeafDataType *) el->child[1];
        assert( ldata );

        elDet_ = ldata->determinant;
        assert( std::abs( elDet_ ) > 0.0 );
        calcedDet_ = true;
      }
      else
      {
        elDet_     = elDeterminant();
        calcedDet_ = true;
      }

      //assert(elDeterminant() > 0.0);
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
  inline bool AlbertaGridGeometry<mydim,cdim,GridImp>::
  builtLocalGeom(const GeometryType &geo, const LocalGeometryType & localGeom,
                 ALBERTA EL_INFO *elInfo,int face)
  {
    face_ = face;
    edge_   = 0;
    vertex_ = 0;
    builtinverse_ = false;
    builtElMat_   = false;

    if( elInfo != NULL )
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

  // built Geometry
  template <int mydim, int cdim, class GridImp>
  inline void AlbertaGridGeometry<mydim,cdim,GridImp>::
  buildGeomInFather(const int child, const int orientation )
  {
    initGeom();
    // reset coordinate vectors
    coord_ = 0.0;

    assert( (child == 0) || (child == 1) );
    if(mydim == 2)
    {
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

    if(mydim == 3)
    {
      if( child == 0 )
      {
        coord_[0] = 0.0; // (0,0,0)
        coord_[1][1] = 1.0; // (0,1,0)
        coord_[2][2] = 1.0; // (0,0,1)
        coord_[3][0] = 0.5; // (0.5,0,0)
      }
      if( child == 1 )
      {
        coord_[0][0] = 1.0; // (1,0,0)

        // depending on orientation the coood 1 and 2 are swaped
        // see Alberta Docu page 5
        if(orientation > 0)
        {
          coord_[1][1] = 1.0; // (0,1,0)
          coord_[2][2] = 1.0; // (0,0,1)
        }
        else
        {
          coord_[1][2] = 1.0; // (0,1,0)
          coord_[2][1] = 1.0; // (0,0,1)
        }
        coord_[3][0] = 0.5; // (0.5,0,0)
      }
    }

    // its a child of the reference element ==> det = 0.5
    elDet_ = elDeterminant();
    calcedDet_ = true;

    if( elDet_ > 0.0 ) return;
    DUNE_THROW(NotImplemented,"wrong dimension given!");
  }

}

#endif
