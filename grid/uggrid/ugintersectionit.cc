// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
template<class GridImp>
inline bool UGGridLevelIntersectionIterator< GridImp >::neighbor() const
{
  return UG_NS<dim>::NbElem(center_, neighborCount_) != NULL;
}

template<class GridImp>
inline bool
UGGridLevelIntersectionIterator<GridImp>::boundary() const
{
  return UG_NS<dim>::Side_On_Bnd(center_, neighborCount_);
}

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
  typename UG_NS<dim>::Element *other,*self;

  // if we have a neighbor on this level, then return it
  if (UG_NS<dim>::NbElem(center_, neighborCount_)!=NULL)
  {
    other = UG_NS<dim>::NbElem(center_, neighborCount_);
    self = center_;
  }
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
numberInSelf ()  const
{
  return UGGridRenumberer<dim>::facesUGtoDUNE(neighborCount_, UG_NS<dimworld>::Sides_Of_Elem(center_));
}

template< class GridImp>
inline int UGGridLevelIntersectionIterator<GridImp>::
numberInNeighbor () const
{
  typename UG_NS<dim>::Element *other, *self;

  // if we have a neighbor on this level, then return it
  if (UG_NS<dim>::NbElem(center_, neighborCount_)!=NULL)
  {
    other = UG_NS<dim>::NbElem(center_, neighborCount_);
    self = center_;
  }
  else
  {
    // now go down the stack of copies to find a lower level leaf neighbor
    typename UG_NS<dim>::Element* father_ = UG_NS<dim>::EFather(center_);
    while (father_!=0)
    {
      if (!UG_NS<dim>::hasCopy(father_))
        DUNE_THROW(GridError,"no neighbor found");
      if (UG_NS<dim>::NbElem(father_, neighborCount_)!=NULL)             // check existence of neighbor
        if (UG_NS<dim>::isLeaf(UG_NS<dim>::NbElem(father_, neighborCount_)))
        {
          other = UG_NS<dim>::NbElem(father_, neighborCount_);
          self = father_;
          break;
        }
      // try father
      father_ = UG_NS<dim>::EFather(father_);
    }
    if (father_==0)
      DUNE_THROW(GridError,"no neighbor found");
  }

  // we have other and self
  const int nSides = UG_NS<dim>::Sides_Of_Elem(other);
  int i;
  for (i=0; i<nSides; i++)
    if (UG_NS<dim>::NbElem(other,i) == self)
      break;

  // now we have to renumber the side i
  return UGGridRenumberer<dim>::facesUGtoDUNE(i, nSides);
}


// /////////////////////////////////////////////////////////////////////////////
//   Implementations for the class UGGridLeafIntersectionIterator
// /////////////////////////////////////////////////////////////////////////////

// returns a neighbor that is a leaf or nothing (neighbor might be on the same level)
// works only on leaf elements!
template<class GridImp>
inline typename UG_NS<GridImp::dimension>::Element* UGGridLeafIntersectionIterator< GridImp >::getLeafNeighbor () const
{
  // if the level neighbor exists and is a leaf then return it
  typename UG_NS<dim>::Element* p = UG_NS<dim>::NbElem(center_, neighborCount_);
  if (p!=NULL)
    if (UG_NS<dim>::isLeaf(p))
      return p;

  // now I must be a leaf to proceed
  if (!UG_NS<dim>::isLeaf(center_))
    return NULL;

  // up or down ?
  if (p==NULL)
  {
    // I am a leaf and the neighbor does not exist: go down
    typename UG_NS<dim>::Element* father_ = UG_NS<GridImp::dimensionworld>::EFather(center_);
    while (father_!=0)
    {
      if (!UG_NS<dim>::hasCopy(father_)) break;             // father must be a copy
      if (UG_NS<dim>::NbElem(father_, neighborCount_)!=NULL)             // check existence of neighbor
        if (UG_NS<dim>::isLeaf(UG_NS<dim>::NbElem(father_, neighborCount_)))                 // check leafness
          return UG_NS<dim>::NbElem(father_, neighborCount_);
      father_ = UG_NS<dim>::EFather(father_);
    }
  }
  else
  {
    // I am a leaf and the neighbor exists and the neighbor is not a leaf: go up
    while (p!=0)
    {
      if (!UG_NS<dim>::hasCopy(p)) break;             // element must be copy refined
      typename UG_NS<dim>::Element *sons[32];
      UG_NS<dim>::GetSons(p,sons);
      p = sons[0];
      if (UG_NS<dim>::isLeaf(p))
        return p;
    }
  }

  // nothing found, return 0 (might be a processor boundary)
  return NULL;
}

template<class GridImp>
inline bool UGGridLeafIntersectionIterator< GridImp >::neighbor() const
{
  return getLeafNeighbor() != NULL;
}

template<class GridImp>
inline bool
UGGridLeafIntersectionIterator<GridImp>::boundary() const
{
  return UG_NS<dim>::Side_On_Bnd(center_, neighborCount_);
}

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

template< class GridImp>
inline const typename UGGridLeafIntersectionIterator<GridImp>::LocalGeometry&
UGGridLeafIntersectionIterator<GridImp>::
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
inline const typename UGGridLeafIntersectionIterator<GridImp>::Geometry&
UGGridLeafIntersectionIterator<GridImp>::
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
inline const typename UGGridLeafIntersectionIterator<GridImp>::LocalGeometry&
UGGridLeafIntersectionIterator<GridImp>::
intersectionNeighborLocal() const
{
  typename UG_NS<dim>::Element *other,*self;

  // if we have a neighbor on this level, then return it
  if (UG_NS<dim>::NbElem(center_, neighborCount_)!=NULL)
  {
    other = UG_NS<dim>::NbElem(center_, neighborCount_);
    self = center_;
  }
  else
  {
    // now go down the stack of copies to find a lower level leaf neighbor
    typename UG_NS<dim>::Element* father_ = UG_NS<dim>::EFather(center_);
    while (father_!=0)
    {
      if (!UG_NS<dim>::hasCopy(father_))
        DUNE_THROW(GridError,"no neighbor found");
      if (UG_NS<dim>::NbElem(father_, neighborCount_)!=NULL)             // check existence of neighbor
        if (UG_NS<dim>::isLeaf(UG_NS<dim>::NbElem(father_, neighborCount_)))
        {
          other = UG_NS<dim>::NbElem(father_, neighborCount_);
          self = father_;
          break;
        }
      // try father
      father_ = UG_NS<dim>::EFather(father_);
    }
    if (father_==0)
      DUNE_THROW(GridError,"no neighbor found");
  }

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

  return neighLocal_;
}

template< class GridImp>
inline int UGGridLeafIntersectionIterator<GridImp>::
numberInSelf ()  const
{
  return UGGridRenumberer<dim>::facesUGtoDUNE(neighborCount_, UG_NS<dimworld>::Sides_Of_Elem(center_));
}

template< class GridImp>
inline int UGGridLeafIntersectionIterator<GridImp>::
numberInNeighbor () const
{
  typename UG_NS<dim>::Element *other, *self;

  // if we have a neighbor on this level, then return it
  if (UG_NS<dim>::NbElem(center_, neighborCount_)!=NULL)
  {
    other = UG_NS<dim>::NbElem(center_, neighborCount_);
    self = center_;
  }
  else
  {
    // now go down the stack of copies to find a lower level leaf neighbor
    typename UG_NS<dim>::Element* father_ = UG_NS<dim>::EFather(center_);
    while (father_!=0)
    {
      if (!UG_NS<dim>::hasCopy(father_))
        DUNE_THROW(GridError,"no neighbor found");
      if (UG_NS<dim>::NbElem(father_, neighborCount_)!=NULL)             // check existence of neighbor
        if (UG_NS<dim>::isLeaf(UG_NS<dim>::NbElem(father_, neighborCount_)))
        {
          other = UG_NS<dim>::NbElem(father_, neighborCount_);
          self = father_;
          break;
        }
      // try father
      father_ = UG_NS<dim>::EFather(father_);
    }
    if (father_==0)
      DUNE_THROW(GridError,"no neighbor found");
  }

  // we have other and self
  const int nSides = UG_NS<dim>::Sides_Of_Elem(other);
  int i;
  for (i=0; i<nSides; i++)
    if (UG_NS<dim>::NbElem(other,i) == self)
      break;

  // now we have to renumber the side i
  return UGGridRenumberer<dim>::facesUGtoDUNE(i, nSides);
}
