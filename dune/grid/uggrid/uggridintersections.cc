// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridintersections.hh>

#include <set>


template<class GridImp>
const typename Dune::UGGridLevelIntersection<GridImp>::WorldVector&
Dune::UGGridLevelIntersection<GridImp>::outerNormal
  (const typename Dune::UGGridLevelIntersection<GridImp>::FaceVector& local) const
{
  // Implementation uses the GenericGeometry facilities.  For this we need to set up a
  // generic geometry of the _inside_ element, which then has a method computing the normal.
  std::vector<FieldVector<typename GridImp::ctype, dim> > corners(inside()->template count<dim>());

  for (int i=0; i<corners.size(); i++)
    corners[i] = inside()->geometry().corner(i);

  GenericGeometry::BasicGeometry<dim, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,dim,dim> > insideGeometry(inside()->type(), corners);

  // Actually compute the normal.
  // Note: The local coordinates that have to be provided are with respect to the inside() element
  outerNormal_ = insideGeometry.normal(indexInInside(),
                                       geometryInInside().global(local));

  return outerNormal_;
}

template< class GridImp>
const typename Dune::UGGridLevelIntersection<GridImp>::LocalGeometry&
Dune::UGGridLevelIntersection<GridImp>::geometryInInside () const
{
  int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_, neighborCount_);
  std::vector<FieldVector<UGCtype,dim> > coordinates(numCornersOfSide);
  GeometryType intersectionGeometryType( (numCornersOfSide==4) ? GeometryType::cube : GeometryType::simplex ,dim-1);

  for (int i=0; i<numCornersOfSide; i++) {

    // get number of corner in UG's numbering system
    int ugIdx     = UGGridRenumberer<dim-1>::verticesDUNEtoUG(i, intersectionGeometryType);
    int cornerIdx = UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, ugIdx);

    // get the corners local coordinates
    UG_NS<dim>::getCornerLocal(center_,cornerIdx,coordinates[i]);

  }

  GridImp::getRealImplementation(selfLocal_).setup(intersectionGeometryType, coordinates);

  return selfLocal_;
}

template< class GridImp>
const typename Dune::UGGridLevelIntersection<GridImp>::Geometry&
Dune::UGGridLevelIntersection<GridImp>::geometry () const
{
  int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_, neighborCount_);
  std::vector<FieldVector<UGCtype,dim> > coordinates(numCornersOfSide);
  GeometryType intersectionGeometryType( (numCornersOfSide==4) ? GeometryType::cube : GeometryType::simplex ,dim-1);

  for (int i=0; i<numCornersOfSide; i++) {

    int ugIdx     = UGGridRenumberer<dim-1>::verticesDUNEtoUG(i, intersectionGeometryType);
    int cornerIdx = UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, ugIdx);
    typename UG_NS<dim>::Node* node = UG_NS<dim>::Corner(center_, cornerIdx);

    for (int j=0; j<dim; j++)
      coordinates[i][j] = node->myvertex->iv.x[j];

  }

  GridImp::getRealImplementation(neighGlob_).setup(intersectionGeometryType, coordinates);

  return neighGlob_;
}

template< class GridImp>
const typename Dune::UGGridLevelIntersection<GridImp>::LocalGeometry&
Dune::UGGridLevelIntersection<GridImp>::geometryInOutside () const
{
  typename UG_NS<dim>::Element *other;

  // if we have a neighbor on this level, then return it
  other = UG_NS<dim>::NbElem(center_, neighborCount_);
  if (!other)
    DUNE_THROW(GridError,"no neighbor found");

  // ///////////////////////////////////////
  // go on and get the local coordinates
  // ///////////////////////////////////////
  int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_,neighborCount_);
  std::vector<FieldVector<UGCtype,dim> > coordinates(numCornersOfSide);
  GeometryType intersectionGeometryType( (numCornersOfSide==4) ? GeometryType::cube : GeometryType::simplex ,dim-1);

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
    UG_NS<dim>::getCornerLocal(other,j,coordinates[UGGridRenumberer<dim-1>::verticesUGtoDUNE(i, intersectionGeometryType)]);

  }

  GridImp::getRealImplementation(neighLocal_).setup(intersectionGeometryType, coordinates);

  return neighLocal_;
}

template< class GridImp>
int Dune::UGGridLevelIntersection<GridImp>::indexInOutside () const
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
  const int number = UGGridRenumberer<dim>::facesUGtoDUNE(i, nSides);

  typedef GenericGeometry::MapNumberingProvider< dim > Numbering;
  const unsigned int tid = GenericGeometry::topologyId( outside()->type() );
  return Numbering::template dune2generic< 1 >( tid, number );
}


// /////////////////////////////////////////////////////////////////////////////
//   Implementations for the class UGGridLeafIntersection
// /////////////////////////////////////////////////////////////////////////////

/** \bug Doesn't work properly for nonplanar nonconforming quadrilateral faces,
   because the local coordinates are interpreted as being with respect to the element
   face.  Instead, they should be interpreted with respect to the intersection.
   If the face is flat this doesn't matter.
 */
template<class GridImp>
const typename Dune::UGGridLeafIntersection<GridImp>::WorldVector&
Dune::UGGridLeafIntersection<GridImp>::outerNormal
  (const typename Dune::UGGridLeafIntersection<GridImp>::FaceVector& local) const
{
  // Implementation uses the GenericGeometry facilities.  For this we need to set up a
  // generic geometry of the _inside_ element, which then has a method computing the normal.
  std::vector<FieldVector<typename GridImp::ctype, dim> > corners(inside()->template count<dim>());

  for (int i=0; i<corners.size(); i++)
    corners[i] = inside()->geometry().corner(i);

  GenericGeometry::BasicGeometry<dim, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,dim,dim> > insideGeometry(inside()->type(), corners);

  // Actually compute the normal.
  // Note: The local coordinates that have to be provided are with respect to the inside() element
  outerNormal_ = insideGeometry.normal(indexInInside(),
                                       geometryInInside().global(local));

  return outerNormal_;
}

template< class GridImp>
const typename Dune::UGGridLeafIntersection<GridImp>::LocalGeometry&
Dune::UGGridLeafIntersection< GridImp >::geometryInInside () const
{
  if (leafSubFaces_[0].first == NULL         // boundary intersection
      // or if this face is the intersection
      || UG_NS<dim>::myLevel(leafSubFaces_[subNeighborCount_].first) <= UG_NS<dim>::myLevel(center_)
      || (UG_NS<dim>::myLevel(leafSubFaces_[subNeighborCount_].first) > UG_NS<dim>::myLevel(center_)
          && leafSubFaces_.size()==1)
      ) {

    // //////////////////////////////////////////////////////
    //   The easy case: a conforming intersection
    // //////////////////////////////////////////////////////

    int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_, neighborCount_);
    std::vector<FieldVector<UGCtype,dim> > coordinates(numCornersOfSide);
    GeometryType intersectionGeometryType( (numCornersOfSide==4) ? GeometryType::cube : GeometryType::simplex ,dim-1);

    for (int i=0; i<numCornersOfSide; i++)
    {
      // get number of corner in UG's numbering system
      int cornerIdx = UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, i);

      // get the corners local coordinates
      UG_NS<dim>::getCornerLocal(center_,cornerIdx,coordinates[UGGridRenumberer<dim-1>::verticesUGtoDUNE(i, intersectionGeometryType)]);

    }

    GridImp::getRealImplementation(selfLocal_).setup(intersectionGeometryType, coordinates);

  } else {

    const typename UG_NS<dim>::Element* other = leafSubFaces_[subNeighborCount_].first;
    int otherSide                             = leafSubFaces_[subNeighborCount_].second;

    int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(other, otherSide);
    std::vector<FieldVector<UGCtype,dim> > coordinates(numCornersOfSide);
    GeometryType intersectionGeometryType( (numCornersOfSide==4) ? GeometryType::cube : GeometryType::simplex ,dim-1);

    for (int i=0; i<numCornersOfSide; i++) {

      // Get world coordinate of other element's vertex
      const UGCtype* worldPos = UG_NS<dim>::Corner(other,UG_NS<dim>::Corner_Of_Side(other,otherSide,i))->myvertex->iv.x;

      // Get the local coordinate with respect to this element
      // coorddim*coorddim is an upper bound for the number of vertices
      UGCtype* cornerCoords[dim*dim];
      UG_NS<dim>::Corner_Coordinates(center_, cornerCoords);

      // Actually do the computation
      /** \todo Why is this const_cast necessary? */
      UG_NS<dim>::GlobalToLocal(UG_NS<dim>::Corners_Of_Elem(center_),
                                const_cast<const double**>(cornerCoords), worldPos,
                                &coordinates[UGGridRenumberer<dim-1>::verticesUGtoDUNE(i, intersectionGeometryType)][0]);

    }

    GridImp::getRealImplementation(selfLocal_).setup(intersectionGeometryType, coordinates);
  }

  return selfLocal_;
}

template< class GridImp>
const typename Dune::UGGridLeafIntersection<GridImp>::Geometry&
Dune::UGGridLeafIntersection< GridImp >::geometry () const
{
  if (leafSubFaces_[0].first == NULL         // boundary intersection
      // or if this face is the intersection
      || UG_NS<dim>::myLevel(leafSubFaces_[subNeighborCount_].first) <= UG_NS<dim>::myLevel(center_)
      || (UG_NS<dim>::myLevel(leafSubFaces_[subNeighborCount_].first) > UG_NS<dim>::myLevel(center_)
          && leafSubFaces_.size()==1)
      ) {

    // //////////////////////////////////////////////////////
    //   The easy case: a conforming intersection
    // //////////////////////////////////////////////////////

    int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_, neighborCount_);
    std::vector<FieldVector<UGCtype,dim> > coordinates(numCornersOfSide);
    GeometryType intersectionGeometryType( (numCornersOfSide==4) ? GeometryType::cube : GeometryType::simplex ,dim-1);

    for (int i=0; i<numCornersOfSide; i++) {

      int cornerIdx = UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, i);
      const typename UG_NS<dim>::Node* node = UG_NS<dim>::Corner(center_, cornerIdx);

      for (int j=0; j<dim; j++)
        coordinates[UGGridRenumberer<dim-1>::verticesUGtoDUNE(i, intersectionGeometryType)][j] = node->myvertex->iv.x[j];

    }

    GridImp::getRealImplementation(neighGlob_).setup(intersectionGeometryType, coordinates);

  } else {

    const typename UG_NS<dim>::Element* other = leafSubFaces_[subNeighborCount_].first;
    int otherSide                             = leafSubFaces_[subNeighborCount_].second;

    int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(other, otherSide);
    std::vector<FieldVector<UGCtype,dim> > coordinates(numCornersOfSide);
    GeometryType intersectionGeometryType( (numCornersOfSide==4) ? GeometryType::cube : GeometryType::simplex ,dim-1);

    for (int i=0; i<numCornersOfSide; i++) {

      // get number of corner in UG's numbering system
      int cornerIdx = UG_NS<dim>::Corner_Of_Side(other, otherSide, i);

      // Get world coordinate of other element's vertex
      const UGCtype* worldPos = UG_NS<dim>::Corner(other,cornerIdx)->myvertex->iv.x;

      // and poke them into the Geometry
      for (int j=0; j<dim; j++)
        coordinates[UGGridRenumberer<dim-1>::verticesUGtoDUNE(i, intersectionGeometryType)][j] = worldPos[j];

    }

    GridImp::getRealImplementation(neighGlob_).setup(intersectionGeometryType, coordinates);

  }

  return neighGlob_;
}

/** \todo Needs to be checked for the nonconforming case */
template< class GridImp>
const typename Dune::UGGridLeafIntersection<GridImp>::LocalGeometry&
Dune::UGGridLeafIntersection< GridImp >::geometryInOutside () const
{
  if (leafSubFaces_[0].first == NULL)
    DUNE_THROW(GridError, "There is no neighbor!");

  if (  // if this face is the intersection
    UG_NS<dim>::myLevel(leafSubFaces_[subNeighborCount_].first) <= UG_NS<dim>::myLevel(center_)
    || (UG_NS<dim>::myLevel(leafSubFaces_[subNeighborCount_].first) > UG_NS<dim>::myLevel(center_)
        && leafSubFaces_.size()==1)
    ) {

    const typename UG_NS<dim>::Element* other = leafSubFaces_[subNeighborCount_].first;

    int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(center_, neighborCount_);
    std::vector<FieldVector<UGCtype,dim> > coordinates(numCornersOfSide);
    GeometryType intersectionGeometryType( (numCornersOfSide==4) ? GeometryType::cube : GeometryType::simplex ,dim-1);

    for (int i=0; i<numCornersOfSide; i++) {

      // get number of corner in UG's numbering system
      int cornerIdx = UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, i);

      // Get world coordinate of this element's vertex
      const UGCtype* worldPos = UG_NS<dim>::Corner(center_,cornerIdx)->myvertex->iv.x;

      // Get the local coordinate with respect to the other element
      // coorddim*coorddim is an upper bound for the number of vertices
      UGCtype* cornerCoords[dim*dim];
      UG_NS<dim>::Corner_Coordinates(other, cornerCoords);

      // Actually do the computation
      /** \todo Why is this const_cast necessary? */
      UG_NS<dim>::GlobalToLocal(UG_NS<dim>::Corners_Of_Elem(other),
                                const_cast<const double**>(cornerCoords), worldPos,
                                &coordinates[UGGridRenumberer<dim-1>::verticesUGtoDUNE(i, intersectionGeometryType)][0]);

    }

    GridImp::getRealImplementation(neighLocal_).setup(intersectionGeometryType, coordinates);

  } else {

    const typename UG_NS<dim>::Element* other = leafSubFaces_[subNeighborCount_].first;
    int otherSide                             = leafSubFaces_[subNeighborCount_].second;

    int numCornersOfSide = UG_NS<dim>::Corners_Of_Side(other, otherSide);
    std::vector<FieldVector<UGCtype,dim> > coordinates(numCornersOfSide);
    GeometryType intersectionGeometryType( (numCornersOfSide==4) ? GeometryType::cube : GeometryType::simplex ,dim-1);

    for (int i=0; i<numCornersOfSide; i++) {

      // get the local coordinate of j-th corner
      int v = UG_NS<dim>::Corner_Of_Side(other,otherSide,i);
      UG_NS<dim>::getCornerLocal(other, v, coordinates[UGGridRenumberer<dim-1>::verticesUGtoDUNE(i, intersectionGeometryType)]);

    }

    GridImp::getRealImplementation(neighLocal_).setup(intersectionGeometryType, coordinates);

  }

  return neighLocal_;
}

template< class GridImp>
int Dune::UGGridLeafIntersection<GridImp>::indexInOutside () const
{
  if (leafSubFaces_[subNeighborCount_].first == NULL)
    DUNE_THROW(GridError,"There is no neighbor!");

  const int nSides = UG_NS<dim>::Sides_Of_Elem(leafSubFaces_[subNeighborCount_].first);

  assert(leafSubFaces_[subNeighborCount_].second < nSides);

  // Renumber to DUNE numbering
  const int number = UGGridRenumberer<dim>::facesUGtoDUNE(leafSubFaces_[subNeighborCount_].second, nSides);

  typedef GenericGeometry::MapNumberingProvider< dim > Numbering;
  const unsigned int tid = GenericGeometry::topologyId( outside()->type() );
  return Numbering::template dune2generic< 1 >( tid, number );
}

template< class GridImp>
void Dune::UGGridLeafIntersection<GridImp>::constructLeafSubfaces() {

  // Do nothing if level neighbor doesn't exit
  typename UG_NS<dim>::Element* levelNeighbor = UG_NS<dim>::NbElem(center_, neighborCount_);

  // If the level neighbor exists and is leaf, then there is only a single leaf intersection
  if (levelNeighbor != NULL && UG_NS<dim>::isLeaf(levelNeighbor)) {
    leafSubFaces_.resize(1);
    leafSubFaces_[0] = Face(levelNeighbor, numberInNeighbor(center_, levelNeighbor));
    return;
  }

  // If the level neighbor does not exist, then leaf intersections exist only with neighbors
  // on lower levels, if they exist at all.  Therefore we descend in the hierarchy towards
  // the coarsest grid until we have found a level neighbor.
  if (levelNeighbor == NULL) {

    leafSubFaces_.resize(1);

    // I am a leaf and the neighbor does not exist: go down
    const typename UG_NS<dim>::Element* me = center_;
    const typename UG_NS<dim>::Element* father = UG_NS<GridImp::dimensionworld>::EFather(center_);

    int side = neighborCount_;
    int fatherSide = -1;

    while (father != NULL) {

      // //////////////////////////////////////////////////////////////
      //   Find the topological father face
      // //////////////////////////////////////////////////////////////
      if (dim==2) {

        // Get the two nodes
        const typename UG_NS<dim>::Node* n0 = UG_NS<dim>::Corner(me,UG_NS<dim>::Corner_Of_Side(me, side, 0));
        const typename UG_NS<dim>::Node* n1 = UG_NS<dim>::Corner(me,UG_NS<dim>::Corner_Of_Side(me, side, 1));

        const typename UG_NS<dim>::Node* fatherN0, *fatherN1;
        assert(!(UG::D2::ReadCW(n0, UG::D2::NTYPE_CE) == UG::D2::MID_NODE
                 && UG::D2::ReadCW(n1, UG::D2::NTYPE_CE) == UG::D2::MID_NODE));

        if (UG::D2::ReadCW(n1, UG::D2::NTYPE_CE) == UG::D2::MID_NODE) {

          // n0 exists on the coarser level, but n1 doesn't
          const typename UG_NS<dim>::Edge* fatherEdge = (const typename UG_NS<dim>::Edge*)n1->father;
          fatherN0 = fatherEdge->links[0].nbnode;
          fatherN1 = fatherEdge->links[1].nbnode;

        } else if (UG::D2::ReadCW(n0, UG::D2::NTYPE_CE) == UG::D2::MID_NODE) {

          // n1 exists on the coarser level, but n0 doesn't
          const typename UG_NS<dim>::Edge* fatherEdge = (const typename UG_NS<dim>::Edge*)n0->father;
          fatherN0 = fatherEdge->links[0].nbnode;
          fatherN1 = fatherEdge->links[1].nbnode;

        } else {

          // This edge is a copy
          fatherN0 = (const typename UG_NS<dim>::Node*)n0->father;
          fatherN1 = (const typename UG_NS<dim>::Node*)n1->father;

        }

        // Find the corresponding side on the father element
        int i;
        for (i=0; i<UG_NS<dim>::Sides_Of_Elem(father); i++) {
          if ( (fatherN0 == UG_NS<dim>::Corner(father,UG_NS<dim>::Corner_Of_Side(father, i, 0))
                && fatherN1 == UG_NS<dim>::Corner(father,UG_NS<dim>::Corner_Of_Side(father, i, 1)))
               || (fatherN0 == UG_NS<dim>::Corner(father,UG_NS<dim>::Corner_Of_Side(father, i, 1))
                   && fatherN1 == UG_NS<dim>::Corner(father,UG_NS<dim>::Corner_Of_Side(father, i, 0))))
            break;
        }

        assert (i<UG_NS<dim>::Sides_Of_Elem(father));
        fatherSide = i;

      } else {        //  dim==3

        // Get the nodes
        int nNodes = UG_NS<dim>::Corners_Of_Side(me,side);
        const typename UG_NS<dim>::Node* n[nNodes];
        for (int i=0; i<nNodes; i++)
          n[i] = UG_NS<dim>::Corner(me,UG_NS<dim>::Corner_Of_Side(me, side, i));

        std::set<const typename UG_NS<dim>::Node*> fatherNodes;          // No more than four father nodes

        for (int i=0; i<nNodes; i++) {

          switch (UG::D3::ReadCW(n[i], UG::D3::NTYPE_CE)) {

          case UG::D3::CORNER_NODE :
            fatherNodes.insert((const typename UG_NS<dim>::Node*)n[i]->father);
            break;
          case UG::D3::MID_NODE :
            fatherNodes.insert( ((const typename UG_NS<dim>::Edge*)n[i]->father)->links[0].nbnode );
            fatherNodes.insert( ((const typename UG_NS<dim>::Edge*)n[i]->father)->links[1].nbnode );
            break;
          default :
            break;
            // Do nothing
          }

        }

        /* When explicitly using anisotropic refinement rules without green closure there
           may be the case that a quad face is split into two quads and a triangle.
           In that case the triangle has two CORNER_NODEs and one MID_NODES.  The current
           code does not know how to handle this situation. */
        if (fatherNodes.size() < 3)
          DUNE_THROW(NotImplemented, "Anisotropic nonconforming grids are not fully implemented!");

        // Find the corresponding side on the father element
        int i;
        for (i=0; i<UG_NS<dim>::Sides_Of_Elem(father); i++) {
          unsigned int found = 0;
          typename std::set<const typename UG_NS<dim>::Node*>::iterator fNIt = fatherNodes.begin();
          for (; fNIt != fatherNodes.end(); ++fNIt)
            for (int k=0; k<UG_NS<dim>::Corners_Of_Side(father,i); k++)
              if (*fNIt == UG_NS<dim>::Corner(father,UG_NS<dim>::Corner_Of_Side(father, i, k))) {
                found++;
                break;
              }

          if (found==fatherNodes.size()) {
            fatherSide = i;
            break;
          }

        }
      }

      // Do we have a neighbor on across this side?
      typename UG_NS<dim>::Element* otherElement = UG_NS<dim>::NbElem(father, fatherSide);
      if (otherElement) {
        const int nSides = UG_NS<dim>::Sides_Of_Elem(otherElement);
        for (int i=0; i<nSides; i++)
          if (UG_NS<dim>::NbElem(otherElement,i) == father) {
            leafSubFaces_[0] = Face(otherElement, i);
            return;
          }

      }

      // Go further down
      me     = father;
      father = UG_NS<dim>::EFather(father);
      side   = fatherSide;

    }

    // Nothing found
    leafSubFaces_[0] = Face( (typename UG_NS<dim>::Element*)NULL, 0);
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

// Explicit template instantiations to compile the stuff in this file
template class Dune::UGGridLevelIntersection<const Dune::UGGrid<2> >;
template class Dune::UGGridLevelIntersection<const Dune::UGGrid<3> >;

template class Dune::UGGridLeafIntersection<const Dune::UGGrid<2> >;
template class Dune::UGGridLeafIntersection<const Dune::UGGrid<3> >;
