// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
template<class GridImp>
inline const FieldVector<typename GridImp::ctype, GridImp::dimensionworld>&
UGGridLevelIntersectionIterator <GridImp>::outerNormal (const FieldVector<UGCtype, GridImp::dimension-1>& local) const
{
  // //////////////////////////////////////////////////////
  //   Implementation for 3D
  // //////////////////////////////////////////////////////

  if (dim == 3) {

    if (UG_NS<dim>::Corners_Of_Side(center_, neighborCount_) == 3) {

      // A triangular intersection.  The normals are constant
      const UGCtype* aPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, 0))->myvertex->iv.x;
      const UGCtype* bPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, 1))->myvertex->iv.x;
      const UGCtype* cPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, 2))->myvertex->iv.x;

      FieldVector<UGCtype, 3> ba, ca;

      for (int i=0; i<3; i++) {
        ba[i] = bPos[i] - aPos[i];
        ca[i] = cPos[i] - aPos[i];
      }

      outerNormal_[0] = ba[1]*ca[2] - ba[2]*ca[1];
      outerNormal_[1] = ba[2]*ca[0] - ba[0]*ca[2];
      outerNormal_[2] = ba[0]*ca[1] - ba[1]*ca[0];

    } else {

      // A quadrilateral: compute the normal in each corner and do bilinear interpolation
      // The cornerNormals array uses UG corner numbering
      FieldVector<UGCtype,3> cornerNormals[4];
      for (int i=0; i<4; i++) {

        // Compute the normal on the i-th corner
        const UGCtype* aPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_,neighborCount_,i))->myvertex->iv.x;
        const UGCtype* bPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_,neighborCount_,(i+1)%4))->myvertex->iv.x;
        const UGCtype* cPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_,neighborCount_,(i+3)%4))->myvertex->iv.x;

        FieldVector<UGCtype, 3> ba, ca;

        for (int j=0; j<3; j++) {
          ba[j] = bPos[j] - aPos[j];
          ca[j] = cPos[j] - aPos[j];
        }

        cornerNormals[i][0] = ba[1]*ca[2] - ba[2]*ca[1];
        cornerNormals[i][1] = ba[2]*ca[0] - ba[0]*ca[2];
        cornerNormals[i][2] = ba[0]*ca[1] - ba[1]*ca[0];
      }

      // Bilinear interpolation
      for (int i=0; i<3; i++)
        outerNormal_[i] = (1-local[0])*(1-local[1])*cornerNormals[0][i]
                          + local[0]     * (1-local[1]) * cornerNormals[1][i]
                          + local[0]     * local[1]     * cornerNormals[2][i]
                          + (1-local[0]) * local[1]     * cornerNormals[3][i];

    }

  } else {     // if (dim==3) ... else

    // //////////////////////////////////////////////////////
    //   Implementation for 2D
    // //////////////////////////////////////////////////////

    // Get the vertices of this side.
    const UGCtype* aPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, 0))->myvertex->iv.x;
    const UGCtype* bPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, 1))->myvertex->iv.x;

    // compute normal
    outerNormal_[0] = bPos[1] - aPos[1];
    outerNormal_[1] = aPos[0] - bPos[0];

  }

  return outerNormal_;
}

template< class GridImp>
inline const typename UGGridLevelIntersectionIterator<GridImp>::LocalGeometry&
UGGridLevelIntersectionIterator<GridImp>::
intersectionSelfLocal() const
{
  int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_, neighborCount_);

  selfLocal_.setNumberOfCorners(numCornersOfSide);

  for (int i=0; i<numCornersOfSide; i++)
  {
    // get number of corner in UG's numbering system
    int cornerIdx = UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, i);

    // we need a temporary to be filled
    FieldVector<UGCtype, dim> tmp;

    // get the corners local coordinates
    UG_NS<dim>::getCornerLocal(center_,cornerIdx,tmp);

    // and poke them into the Geometry
    selfLocal_.setCoords(i,tmp);
  }

  return selfLocal_;
}

template< class GridImp>
inline const typename UGGridLevelIntersectionIterator<GridImp>::Geometry&
UGGridLevelIntersectionIterator<GridImp>::
intersectionGlobal() const
{
  int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_, neighborCount_);

  neighGlob_.setNumberOfCorners(numCornersOfSide);

  for (int i=0; i<numCornersOfSide; i++) {

    int cornerIdx = UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, i);
    typename UG_NS<dim>::Node* node = UG_NS<dim>::Corner(center_, cornerIdx);

    neighGlob_.setCoords(i, node->myvertex->iv.x);

  }

  return neighGlob_;
}

template< class GridImp>
inline const typename UGGridLevelIntersectionIterator<GridImp>::LocalGeometry&
UGGridLevelIntersectionIterator<GridImp>::
intersectionNeighborLocal() const
{
  typename UG_NS<dim>::Element *other;

  // if we have a neighbor on this level, then return it
  if (UG_NS<dim>::NbElem(center_, neighborCount_)!=NULL)
    other = UG_NS<dim>::NbElem(center_, neighborCount_);
  else
    DUNE_THROW(GridError,"no neighbor found");

  // ///////////////////////////////////////
  // go on and get the local coordinates
  // ///////////////////////////////////////
  int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_,neighborCount_);
  neighLocal_.setNumberOfCorners(numCornersOfSide);

  for (int i=0; i<numCornersOfSide; i++) {

    // get the node in this element
    int localCornerNumber = UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, i);
    const typename UG_NS<dim>::Node* node = UG_NS<dim>::Corner(center_,localCornerNumber);

    // get this node's local index in the neighbor element
    int j;
    for (j=0; j<UG_NS<dim>::Corners_Of_Elem(other); j++)
      if (UG_NS<dim>::Corner(other, j) == node)
        break;

    assert(j<UG_NS<dim>::Corners_Of_Elem(other));

    // get the local coordinate there
    FieldVector<UGCtype, dim> tmp;
    UG_NS<dim>::getCornerLocal(other,j,tmp);

    // and poke them into the Geometry
    neighLocal_.setCoords(i,tmp);
  }

  return neighLocal_;
}

template< class GridImp>
inline int UGGridLevelIntersectionIterator<GridImp>::
numberInNeighbor () const
{
  const typename UG_NS<dim>::Element *other;

  // Look for a neighbor on this level
  if ((other = UG_NS<dim>::NbElem(center_, neighborCount_)) == NULL)
    DUNE_THROW(GridError,"There is no neighbor element!");

  // Find the corresponding side in the neighbor element
  const int nSides = UG_NS<dim>::Sides_Of_Elem(other);
  int i;
  for (i=0; i<nSides; i++)
    if (UG_NS<dim>::NbElem(other,i) == center_)
      break;

  // now we have to renumber the side i
  return UGGridRenumberer<dim>::facesUGtoDUNE(i, nSides);
}


// /////////////////////////////////////////////////////////////////////////////
//   Implementations for the class UGGridLeafIntersectionIterator
// /////////////////////////////////////////////////////////////////////////////

/** \todo Needs to be checked for the nonconforming case */
template<class GridImp>
inline const FieldVector<typename GridImp::ctype, GridImp::dimensionworld>&
UGGridLeafIntersectionIterator <GridImp>::outerNormal (const FieldVector<UGCtype, GridImp::dimension-1>& local) const
{
  // //////////////////////////////////////////////////////
  //   Implementation for 3D
  // //////////////////////////////////////////////////////

  if (dim == 3) {

    if (UG_NS<dim>::Corners_Of_Side(center_, neighborCount_) == 3) {

      // A triangular intersection.  The normals are constant
      const UGCtype* aPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, 0))->myvertex->iv.x;
      const UGCtype* bPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, 1))->myvertex->iv.x;
      const UGCtype* cPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, 2))->myvertex->iv.x;

      FieldVector<UGCtype, 3> ba, ca;

      for (int i=0; i<3; i++) {
        ba[i] = bPos[i] - aPos[i];
        ca[i] = cPos[i] - aPos[i];
      }

      outerNormal_[0] = ba[1]*ca[2] - ba[2]*ca[1];
      outerNormal_[1] = ba[2]*ca[0] - ba[0]*ca[2];
      outerNormal_[2] = ba[0]*ca[1] - ba[1]*ca[0];

    } else {

      // A quadrilateral: compute the normal in each corner and do bilinear interpolation
      // The cornerNormals array uses UG corner numbering
      FieldVector<UGCtype,3> cornerNormals[4];
      for (int i=0; i<4; i++) {

        // Compute the normal on the i-th corner
        const UGCtype* aPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_,neighborCount_,i))->myvertex->iv.x;
        const UGCtype* bPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_,neighborCount_,(i+1)%4))->myvertex->iv.x;
        const UGCtype* cPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_,neighborCount_,(i+3)%4))->myvertex->iv.x;

        FieldVector<UGCtype, 3> ba, ca;

        for (int j=0; j<3; j++) {
          ba[j] = bPos[j] - aPos[j];
          ca[j] = cPos[j] - aPos[j];
        }

        cornerNormals[i][0] = ba[1]*ca[2] - ba[2]*ca[1];
        cornerNormals[i][1] = ba[2]*ca[0] - ba[0]*ca[2];
        cornerNormals[i][2] = ba[0]*ca[1] - ba[1]*ca[0];
      }

      // Bilinear interpolation
      for (int i=0; i<3; i++)
        outerNormal_[i] = (1-local[0])*(1-local[1])*cornerNormals[0][i]
                          + local[0]     * (1-local[1]) * cornerNormals[1][i]
                          + local[0]     * local[1]     * cornerNormals[2][i]
                          + (1-local[0]) * local[1]     * cornerNormals[3][i];

    }

  } else {     // if (dim==3) ... else

    // //////////////////////////////////////////////////////
    //   Implementation for 2D
    // //////////////////////////////////////////////////////

    // Get the vertices of this side.
    const UGCtype* aPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, 0))->myvertex->iv.x;
    const UGCtype* bPos = UG_NS<dim>::Corner(center_,UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, 1))->myvertex->iv.x;

    // compute normal
    outerNormal_[0] = bPos[1] - aPos[1];
    outerNormal_[1] = aPos[0] - bPos[0];

  }

  return outerNormal_;
}

/** \todo Needs to be checked for the nonconforming case */
template< class GridImp>
inline const typename UGGridLeafIntersectionIterator<GridImp>::LocalGeometry&
UGGridLeafIntersectionIterator<GridImp>::
intersectionSelfLocal() const
{
  if (leafSubFaces_.size() == 1) {

    // //////////////////////////////////////////////////////
    //   The easy case: a conforming intersection
    // //////////////////////////////////////////////////////

    int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_, neighborCount_);

    selfLocal_.setNumberOfCorners(numCornersOfSide);

    for (int i=0; i<numCornersOfSide; i++)
    {
      // get number of corner in UG's numbering system
      int cornerIdx = UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, i);

      // we need a temporary to be filled
      FieldVector<UGCtype, dim> tmp;

      // get the corners local coordinates
      UG_NS<dim>::getCornerLocal(center_,cornerIdx,tmp);

      // and poke them into the Geometry
      selfLocal_.setCoords(i,tmp);
    }

  } else {
    DUNE_THROW(NotImplemented, "no intersectionSelfLocal() for nonconforming UGGrids");
  }

  return selfLocal_;
}

/** \todo Needs to be checked for the nonconforming case */
template< class GridImp>
inline const typename UGGridLeafIntersectionIterator<GridImp>::Geometry&
UGGridLeafIntersectionIterator<GridImp>::
intersectionGlobal() const
{
  if (leafSubFaces_.size() == 1) {

    // //////////////////////////////////////////////////////
    //   The easy case: a conforming intersection
    // //////////////////////////////////////////////////////

    int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_, neighborCount_);

    neighGlob_.setNumberOfCorners(numCornersOfSide);

    for (int i=0; i<numCornersOfSide; i++) {

      int cornerIdx = UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, i);
      typename UG_NS<dim>::Node* node = UG_NS<dim>::Corner(center_, cornerIdx);

      neighGlob_.setCoords(i, node->myvertex->iv.x);

    }

  } else {
    DUNE_THROW(NotImplemented, "no intersectionGlobal() for nonconforming UGGrids");
  }

  return neighGlob_;
}

/** \todo Needs to be checked for the nonconforming case */
template< class GridImp>
inline const typename UGGridLeafIntersectionIterator<GridImp>::LocalGeometry&
UGGridLeafIntersectionIterator<GridImp>::
intersectionNeighborLocal() const
{
  if (leafSubFaces_.size() == 1) {

    // //////////////////////////////////////////////////////
    //   The easy case: a conforming intersection
    // //////////////////////////////////////////////////////

    const typename UG_NS<dim>::Element *other = leafSubFaces_[subNeighborCount_].first;

    // ///////////////////////////////////////
    // go on and get the local coordinates
    // ///////////////////////////////////////
    int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_,neighborCount_);
    neighLocal_.setNumberOfCorners(numCornersOfSide);

    for (int i=0; i<numCornersOfSide; i++) {

      // get the node in this element
      int localCornerNumber = UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, i);
      const typename UG_NS<dim>::Node* node = UG_NS<dim>::Corner(center_,localCornerNumber);

      // get this node's local index in the neighbor element
      int j;
      for (j=0; j<UG_NS<dim>::Corners_Of_Elem(other); j++)
        // Compare vertices because the nodes may be on different levels, but the nodes are the same
        if (UG_NS<dim>::Corner(other, j)->myvertex == node->myvertex)
          break;

      assert(j<UG_NS<dim>::Corners_Of_Elem(other));

      // get the local coordinate there
      FieldVector<UGCtype, dim> tmp;
      UG_NS<dim>::getCornerLocal(other,j,tmp);

      // and poke them into the Geometry
      neighLocal_.setCoords(i,tmp);
    }

  } else {

    DUNE_THROW(NotImplemented, "no intersectionNeighborLocal for nonconforming UGGrids");

  }

  return neighLocal_;
}

template< class GridImp>
inline int UGGridLeafIntersectionIterator<GridImp>::
numberInNeighbor () const
{
  if (leafSubFaces_[subNeighborCount_].first == NULL)
    DUNE_THROW(GridError,"There is no neighbor!");

  const int nSides = UG_NS<dim>::Sides_Of_Elem(leafSubFaces_[subNeighborCount_].first);

  assert(leafSubFaces_[subNeighborCount_].second < nSides);

  // Renumber to DUNE numbering
  return UGGridRenumberer<dim>::facesUGtoDUNE(leafSubFaces_[subNeighborCount_].second, nSides);
}

template< class GridImp>
inline void UGGridLeafIntersectionIterator<GridImp>::constructLeafSubfaces() {

  // Do nothing if level neighbor doesn't exit
  typename UG_NS<dim>::Element* levelNeighbor = UG_NS<dim>::NbElem(center_, neighborCount_);

  if (levelNeighbor != NULL && UG_NS<dim>::isLeaf(levelNeighbor)) {
    leafSubFaces_.resize(1);
    leafSubFaces_[0] = Face(levelNeighbor, numberInNeighbor(center_, levelNeighbor));
    return;
  }

  // Go down
  if (levelNeighbor == NULL) {

    leafSubFaces_.resize(1);

    // I am a leaf and the neighbor does not exist: go down
    typename UG_NS<dim>::Element* father = UG_NS<GridImp::dimensionworld>::EFather(center_);
    while (father != NULL) {

      /** \bug Does not work for nonconforming grids, because there the father is not a copy */
      if (!UG_NS<dim>::hasCopy(father)) {
        break;         // father must be a copy
      }
      if (UG_NS<dim>::NbElem(father, neighborCount_)!=NULL)       // check existence of neighbor
        if (UG_NS<dim>::isLeaf(UG_NS<dim>::NbElem(father, neighborCount_))) {         // check leafness
          leafSubFaces_[0] = Face(UG_NS<dim>::NbElem(father, neighborCount_),
                                  numberInNeighbor(father, UG_NS<dim>::NbElem(father, neighborCount_)));
          return;
        }

      father = UG_NS<dim>::EFather(father);

    }

    // Nothing found
    leafSubFaces_[0] = Face(NULL, 0);
    return;
  }


  // ///////////////
  //   init list
  // ///////////////

  SLList<Face> list;
  int levelNeighborSide = numberInNeighbor(center_, levelNeighbor);

  int Sons_of_Side = 0;
  typename UG_NS<dim>::Element* SonList[UG_NS<dim>::MAX_SONS];
  int SonSides[UG_NS<dim>::MAX_SONS];

  int rv = Get_Sons_of_ElementSide(levelNeighbor,
                                   levelNeighborSide,
                                   &Sons_of_Side,
                                   SonList,        // the output elements
                                   SonSides,       // Output element side numbers
                                   true,          // Element sons are not precomputed
                                   false,          // ioflag: Obsolete debugging flag
                                   true);

  if (rv!=0)
    DUNE_THROW(GridError, "Get_Sons_of_ElementSide returned with error value " << rv);

  for (int i=0; i<Sons_of_Side; i++)
    list.push_back(Face(SonList[i],SonSides[i]));

  // //////////////////////////////////////////////////
  //   Traverse and collect all children of the side
  // //////////////////////////////////////////////////

  typename SLList<Face>::iterator f = list.begin();
  for (; f!=list.end(); ++f) {

    typename UG_NS<dim>::Element* theElement = f->first;

    int Sons_of_Side = 0;
    typename UG_NS<dim>::Element* SonList[UG_NS<dim>::MAX_SONS];
    int SonSides[UG_NS<dim>::MAX_SONS];

    if (!UG_NS<dim>::isLeaf(theElement)) {

      Get_Sons_of_ElementSide(theElement,
                              f->second,        // Input element side number
                              &Sons_of_Side,       // Number of topological sons of the element side
                              SonList,            // Output elements
                              SonSides,           // Output element side numbers
                              true,
                              false,
                              true);

      for (int i=0; i<Sons_of_Side; i++)
        list.push_back(Face(SonList[i],SonSides[i]));

    }

  }

  // //////////////////////////////
  //   Extract result from stack
  // //////////////////////////////
  leafSubFaces_.resize(0);

  for (f = list.begin(); f!=list.end(); ++f) {

    // Set element
    if (UG_NS<dim>::isLeaf(f->first))
      leafSubFaces_.push_back(*f);

  }

}
