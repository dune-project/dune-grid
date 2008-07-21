// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
//****************************************************************
//
// --UGGridGeometry
//
//****************************************************************

#include <algorithm>

///////////////////////////////////////////////////////////////////////
//
// General implementation of UGGridGeometry mydim-> coorddim
//
///////////////////////////////////////////////////////////////////////


template< int mydim, int coorddim, class GridImp>
inline Dune::GeometryType Dune::UGGridGeometry<mydim,coorddim,GridImp>::type() const
{
  switch (mydim)
  {
  case 0 : return GeometryType(0);
  case 1 : return GeometryType(1);
  case 2 :

    switch (UG_NS<coorddim>::Tag(target_)) {
    case UG::D2::TRIANGLE :
      return GeometryType(GeometryType::simplex,2);
    case UG::D2::QUADRILATERAL :
      return GeometryType(GeometryType::cube,2);
    default :
      DUNE_THROW(GridError, "UGGridGeometry::type():  ERROR:  Unknown type "
                 << UG_NS<coorddim>::Tag(target_) << " found!");
    }

  case 3 :
    switch (UG_NS<coorddim>::Tag(target_)) {

    case UG::D3::TETRAHEDRON :
      return GeometryType(GeometryType::simplex,3);
    case UG::D3::PYRAMID :
      return GeometryType(GeometryType::pyramid,3);
    case UG::D3::PRISM :
      return GeometryType(GeometryType::prism,3);
    case UG::D3::HEXAHEDRON :
      return GeometryType(GeometryType::cube,3);
    default :
      DUNE_THROW(GridError, "UGGridGeometry::type():  ERROR:  Unknown type "
                 << UG_NS<coorddim>::Tag(target_) << " found!");

    }
  }

}

template<int mydim, int coorddim, class GridImp>
inline const Dune::FieldVector<typename GridImp::ctype, coorddim>& Dune::UGGridGeometry<mydim,coorddim,GridImp>::
operator [](int i) const
{
  // This geometry is a vertex
  if (mydim==0)
    return reinterpret_cast<FieldVector<typename GridImp::ctype, coorddim>&>(((typename UG_NS<coorddim>::Node*)target_)->myvertex->iv.x);

  // ////////////////////////////////
  //  This geometry is an element
  // ////////////////////////////////
  assert(mydim==coorddim);

  i = UGGridRenumberer<mydim>::verticesDUNEtoUG(i,type());

  if (mode_==element_mode)
    return reinterpret_cast<FieldVector<typename GridImp::ctype, coorddim>&>(UG_NS<coorddim>::Corner(((typename UG_NS<coorddim>::Element*)target_),i)->myvertex->iv.x);

  return coord_[i];
}

template< int mydim, int coorddim, class GridImp>
inline Dune::FieldVector<typename GridImp::ctype, coorddim> Dune::UGGridGeometry<mydim,coorddim,GridImp>::
global(const FieldVector<UGCtype, mydim>& local) const
{
  FieldVector<UGCtype, coorddim> globalCoord;

  if (mode_==element_mode) {

    UGCtype* cornerCoords[corners()];
    UG_NS<coorddim>::Corner_Coordinates(target_, cornerCoords);

    // Actually do the computation
    UG_NS<coorddim>::Local_To_Global(corners(), cornerCoords, local, globalCoord);

  } else {
    UG_NS<coorddim>::Local_To_Global(corners(), cornerpointers_, local, globalCoord);
  }

  return globalCoord;
}


// Maps a global coordinate within the element to a
// local coordinate in its reference element
template< int mydim, int coorddim, class GridImp>
inline Dune::FieldVector<typename GridImp::ctype, mydim> Dune::UGGridGeometry<mydim,coorddim, GridImp>::
local (const Dune::FieldVector<typename GridImp::ctype, coorddim>& global) const
{
  FieldVector<UGCtype, mydim> result;

  if (mode_==element_mode)
  {
    // coorddim*coorddim is an upper bound for the number of vertices
    UGCtype* cornerCoords[coorddim*coorddim];
    UG_NS<coorddim>::Corner_Coordinates(target_, cornerCoords);

    // Actually do the computation
    /** \todo Why is this const_cast necessary? */
    UG_NS<coorddim>::GlobalToLocal(corners(), const_cast<const double**>(cornerCoords), &global[0], &result[0]);
  }
  else
  {
    // Actually do the computation
    /** \todo Why is this const_cast necessary? */
    UG_NS<coorddim>::GlobalToLocal(corners(), const_cast<const double**>(cornerpointers_), &global[0], &result[0]);
  }

  return result;
}

template<int mydim, int coorddim, class GridImp>
inline bool Dune::UGGridGeometry<mydim,coorddim,GridImp>::
checkInside(const Dune::FieldVector<UGCtype, mydim> &loc) const
{
  switch (mydim) {

  case 0 :  // vertex
    return false;
  case 1 :  // line
    return 0 <= loc[0] && loc[0] <= 1;
  case 2 :

    switch (UG_NS<coorddim>::Tag(target_)) {
    case UG::D2::TRIANGLE :
      return 0 <= loc[0] && 0 <= loc[1] && (loc[0]+loc[1])<=1;
    case UG::D2::QUADRILATERAL :
      return 0 <= loc[0] && loc[0] <= 1
             && 0 <= loc[1] && loc[1] <= 1;
    default :
      DUNE_THROW(GridError, "UGGridGeometry::checkInside():  ERROR:  Unknown type "
                 << UG_NS<coorddim>::Tag(target_) << " found!");
    }

  case 3 :
    switch (UG_NS<coorddim>::Tag(target_)) {

    case UG::D3::TETRAHEDRON :
      return 0 <= loc[0] && 0 <= loc[1] && 0 <= loc[2]
             && (loc[0]+loc[1]+loc[2]) <= 1;
    case UG::D3::PYRAMID :
      return 0 <= loc[0] && 0 <= loc[1] && 0 <= loc[2]
             && (loc[0]+loc[2]) <= 1
             && (loc[1]+loc[2]) <= 1;
    case UG::D3::PRISM :
      return 0 <= loc[0] && 0 <= loc[1]
             && (loc[0]+loc[1])<=1
             && 0 <= loc[2] && loc[2] <= 1;
    case UG::D3::HEXAHEDRON :
      return 0 <= loc[0] && loc[0] <= 1
             && 0 <= loc[1] && loc[1] <= 1
             && 0 <= loc[2] && loc[2] <= 1;
    default :
      DUNE_THROW(GridError, "UGGridGeometry::checkInside():  ERROR:  Unknown type "
                 << UG_NS<coorddim>::Tag(target_) << " found!");

    }
  }

}


template< int mydim, int coorddim, class GridImp>
inline typename GridImp::ctype Dune::UGGridGeometry<mydim,coorddim,GridImp>::
integrationElement (const Dune::FieldVector<typename GridImp::ctype, mydim>& local) const
{
  /** \todo No need to recompute the determinant every time on a simplex */
  return std::abs(1/jacobianInverseTransposed(local).determinant());
}


template< int mydim, int coorddim, class GridImp>
inline const Dune::FieldMatrix<typename GridImp::ctype, mydim,mydim>& Dune::UGGridGeometry<mydim,coorddim, GridImp>::
jacobianInverseTransposed (const Dune::FieldVector<typename GridImp::ctype, mydim>& local) const
{
  if (jacobianInverseIsUpToDate_)
    return jac_inverse_;

  if (mode_==element_mode) {

    // compile array of pointers to corner coordinates
    UGCtype* cornerCoords[corners()];
    UG_NS<coorddim>::Corner_Coordinates(target_, cornerCoords);

    // compute the transformation onto the reference element (or vice versa?)
    UG_NS<coorddim>::Transformation(corners(), cornerCoords, local, jac_inverse_);

  } else
  {
    // compute the transformation onto the reference element (or vice versa?)
    UG_NS<coorddim>::Transformation(corners(), cornerpointers_, local, jac_inverse_);
  }

  if (type().isSimplex())
    jacobianInverseIsUpToDate_ = true;

  return jac_inverse_;
}


///////////////////////////////////////////////////////////////////////
//
// The specializations 1->2, 2->3
// (only methods that are not yet defined in header file)
//
///////////////////////////////////////////////////////////////////////


// Specialization for dim==1, coorddim==2.  This is necessary
// because we specialized the whole class
template <class GridImp>
inline Dune::FieldVector<typename GridImp::ctype, 2> Dune::UGGridGeometry<1,2,GridImp>::
global(const FieldVector<typename GridImp::ctype, 1>& local) const
{
  FieldVector<UGCtype, 2> globalCoord;

  globalCoord[0] = local[0]*coord_[1][0] + (1-local[0])*coord_[0][0];
  globalCoord[1] = local[0]*coord_[1][1] + (1-local[0])*coord_[0][1];

  return globalCoord;
}

// Specialization for dim==2, coorddim==3.  This is necessary
// because we specialized the whole class
template <class GridImp>
inline Dune::FieldVector<typename GridImp::ctype, 3> Dune::UGGridGeometry<2,3,GridImp>::
global(const FieldVector<typename GridImp::ctype, 2>& local) const
{

  FieldVector<UGCtype, 3> result;

  if (elementType_.isSimplex()) {

    for (int i=0; i<3; i++)
      result[i] = (1.0-local[0]-local[1])*coord_[0][i]
                  + local[0]*coord_[1][i]
                  + local[1]*coord_[2][i];

  } else {

    // quadrilateral
    for (int i=0; i<3; i++)
      result[i] = (1.0-local[0])*(1.0-local[1])*coord_[0][i]
                  + local[0]*(1.0-local[1])*coord_[1][i]
                  + local[0]*local[1]*coord_[2][i]
                  + (1.0-local[0])*local[1]*coord_[3][i];

  }

  return result;

}


template <class GridImp>
inline Dune::FieldVector<typename GridImp::ctype, 2> Dune::UGGridGeometry<2,3,GridImp>::
local(const FieldVector<typename GridImp::ctype, 3>& global) const
{
  // The UG method UG_GlobalToLocalBnd pretends to do what this method does,
  // but in reality it is buggy!
  FieldVector<UGCtype,2> result;

  FieldMatrix<UGCtype,2,2> M;
  FieldVector<UGCtype,2> partialDiff;

  FieldVector<UGCtype,3> diff = global;
  diff -= coord_[0];

  if (elementType_.isTriangle()) {

    /* the simplex case:
       Let $l \in R^2$ be the local coordinates and $g \in R^3$ the global ones.
       Let $p_0, p_1, p_2$ be the three triangle corners.
       Then $g = p_0 + D l$, where $D$ is a $3\times 2$-matrix with the first
       column being $p_1 - p_0$ and the second one being $p_2 - p_0$.  In other
       words $g - p_0 = D l$ or $l = D^+ (g - p_0)$, where $D^+$ is the
       Moore-Penrose pseudo inverse $D^+ = (D^T D)^{-1} D^T$.  Using this
       we get that $l$ is the solution of the linear system
       $(D^T D)l = D^T (g-p_0)$ and that is what is implemented here.
     */
    FieldMatrix<UGCtype,2,3> matrixTransp;
    for (int i=0; i<3; i++) {
      matrixTransp[0][i] = coord_[1][i] - coord_[0][i];
      matrixTransp[1][i] = coord_[2][i] - coord_[0][i];
    }

    FieldMatrix<UGCtype,2,2> matrixTranspTimesMatrix(0);
    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
        for (int k=0; k<3; k++)
          matrixTranspTimesMatrix[i][j] += matrixTransp[i][k]*matrixTransp[j][k];

    FieldVector<UGCtype,2> scaledDiff(0);
    matrixTransp.umv(diff,scaledDiff);
    matrixTranspTimesMatrix.solve(result, scaledDiff);

    return result;

  } else {
    assert(elementType_.isQuadrilateral());

    //  Newton - Iteration zum Invertieren der Abbildung f.
    double err;
    FieldVector<UGCtype,3> map(0);
    int count = 0 ;

    FieldMatrix<double,3,3> Df;                // set in method map2worldlinear
    Dune::FieldVector<double,3> normal_;       // in method map2worldnormal
    Dune::FieldMatrix<double,4,3> _b;
    Dune::FieldMatrix<double,3,3> _n;

    // buildmapping
    for (int i=0; i<3; i++) {

      _b [0][i] = coord_[0] [i] ;
      _b [1][i] = coord_[1] [i] - coord_[0] [i] ;
      _b [2][i] = coord_[3] [i] - coord_[0] [i] ;
      _b [3][i] = coord_[2] [i] - coord_[3] [i] - _b [1][i] ;

    }

    _n [0][0] = _b [1][1] * _b [2][2] - _b [1][2] * _b [2][1] ;
    _n [0][1] = _b [1][2] * _b [2][0] - _b [1][0] * _b [2][2] ;
    _n [0][2] = _b [1][0] * _b [2][1] - _b [1][1] * _b [2][0] ;
    _n [1][0] = _b [1][1] * _b [3][2] - _b [1][2] * _b [3][1] ;
    _n [1][1] = _b [1][2] * _b [3][0] - _b [1][0] * _b [3][2] ;
    _n [1][2] = _b [1][0] * _b [3][1] - _b [1][1] * _b [3][0] ;
    _n [2][0] = _b [3][1] * _b [2][2] - _b [3][2] * _b [2][1] ;
    _n [2][1] = _b [3][2] * _b [2][0] - _b [3][0] * _b [2][2] ;
    _n [2][2] = _b [3][0] * _b [2][1] - _b [3][1] * _b [2][0] ;

    do {
      FieldVector<UGCtype,3> upd ;
      //map2worldnormal (map_[0],map_[1],map_[2], upd, normal_,_b, _n) ;
      for (int i=0; i<3; i++)
        normal_ [i] = -(_n [0][i] + _n [1][i] * map[0] + _n [2][i] * map[1]);

      for (int i=0; i<3; i++)
        upd[i] = _b [0][i] + map[0] * _b [1][i] + map[1] * _b [2][i] + map[0]*map[1] * _b [3][i] + map[2]*normal_[0];

      //inverse (map_, normal_, _b, _n, Df,Dfi) ;
      for (int i=0; i<3; i++) {
        Df[i][0] = _b [1][i] + map[1] * _b [3][i]+ map[2]*_n[1][i] ;
        Df[i][1] = _b [2][i] + map[0] * _b [3][i]+ map[2]*_n[2][i] ;
        Df[i][2] = normal_[i];
      }

      upd -= global;

      FieldVector<UGCtype,3> c(0);
      Df.solve(c,upd);

      map -= c;

      err = c.two_norm();

      if (count++ > 100)
        DUNE_THROW(NotImplemented, "Newton solver did not converge!");

    } while (err > 1e-5) ;

    result[0]=map[0];
    result[1]=map[1];

  }

  return result;
}

template <class GridImp>
inline typename GridImp::ctype Dune::UGGridGeometry<1,2,GridImp>::
integrationElement (const Dune::FieldVector<typename GridImp::ctype, 1>& local) const
{
  // We could call UG_NS<2>::SurfaceElement, but this is faster, and not more complicated
  FieldVector<UGCtype, 2> diff = coord_[0];
  diff -= coord_[1];
  return diff.two_norm();
}

template <class GridImp>
const Dune::FieldMatrix<typename GridImp::ctype,2,1>& Dune::UGGridGeometry<1,2,GridImp>::
jacobianInverseTransposed(const Dune::FieldVector<typename GridImp::ctype, 1>& local) const
{
  DUNE_THROW(NotImplemented, "Not implemented yet!");
  return jacobianInverseTransposed_;
}

template <class GridImp>
inline typename GridImp::ctype Dune::UGGridGeometry<2,3,GridImp>::
integrationElement (const Dune::FieldVector<typename GridImp::ctype, 2>& local) const
{
  // Check whether the memory layout of double[4][3] and std::array<FieldVector<double,3>,4>
  // match.  This may not be the case because there are people who like to play
  // with the alignment settings of FieldVector.  The conditional is static and
  // does not cost run-time.
  typedef typename GridImp::ctype ctype;
  if (sizeof(FieldVector<ctype,3>)==3*sizeof(ctype)) {

    // Memory layout matches.  Simply cast
    return UG_NS<3>::SurfaceElement(corners(), (const ctype(*)[3])&coord_,&local[0]);

  } else {

    // Memory layout does not match.  We need to copy into a temporary object.
    ctype tmp[4][3];
    for (int i=0; i<4; i++)
      for (int j=0; j<3; j++)
        tmp[i][j] = coord_[i][j];

    return UG_NS<3>::SurfaceElement(corners(), tmp, &local[0]);

  }
}

template <class GridImp>
const Dune::FieldMatrix<typename GridImp::ctype,3,2>& Dune::UGGridGeometry<2,3,GridImp>::
jacobianInverseTransposed(const Dune::FieldVector<typename GridImp::ctype, 2>& local) const
{
  DUNE_THROW(NotImplemented, "Not implemented yet!");
  return jacobianInverseTransposed_;
}
