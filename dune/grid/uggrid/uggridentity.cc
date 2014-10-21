// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridentity.hh>
#include <dune/grid/uggrid/uggridrenumberer.hh>

//*************************************************************************
//
//  --UGGridEntity
//  --Entity
//
//*************************************************************************


//
//  codim > 0
//
//*********************************************************************


//*****************************************************************8
// count
template <int codim, int dim, class GridImp>
template <int cc>
int Dune::UGGridEntity<codim,dim,GridImp>::count () const
{
  DUNE_THROW(GridError, "UGGridEntity<" << codim << ", " << dim
                                        << ">::count() not implemented yet!");
  return -1;
}

template< int codim, int dim, class GridImp>
Dune::GeometryType Dune::UGGridEntity<codim,dim,GridImp>::type() const
{
  static_assert(codim!=0, "The code for codim==0 is in the corresponding class specialization");

  // Vertex
  if (dim-codim == 0)
    return GeometryType(0);

  // Edge
  if (dim-codim == 1)
    return GeometryType(1);

  // Face in 3d
  switch (UG_NS<dim>::Tag(target_)) {
  case UG::D2::TRIANGLE :
    return GeometryType(GeometryType::simplex,2);
  case UG::D2::QUADRILATERAL :
    return GeometryType(GeometryType::cube,2);
  default :
    DUNE_THROW(GridError, "UGGridGeometry::type():  ERROR:  Unknown type "
               << UG_NS<dim>::Tag(target_) << " found!");
  }

}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//     Specializations for codim == 0                                     //
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template< int dim, class GridImp>
bool Dune::UGGridEntity < 0, dim ,GridImp >::isNew () const
{
  return UG_NS<dim>::ReadCW(target_, UG_NS<dim>::NEWEL_CE);
}

template< int dim, class GridImp>
bool Dune::UGGridEntity < 0, dim ,GridImp >::mightVanish () const
{
  if ((!isRegular()) || (UG_NS<dim>::ReadCW(target_, UG_NS<dim>::COARSEN_CE)))
    return true;

  // maybe it isn't marked, but a sibling is
  if(hasFather())
  {
    typename UG_NS<dim>::Element* father_ = UG_NS<dim>::EFather(target_);
    typename UG_NS<dim>::Element* sons_[UG_NS<dim>::MAX_SONS];
    UG_NS<dim>::GetSons(father_,sons_);
    for (int i = 0; i < UG_NS<dim>::MAX_SONS; ++i)
    {
      if (sons_[i] == NULL) // last son visited
        break;
      if (!UG_NS<dim>::isRegular(sons_[i]) || UG_NS<dim>::ReadCW(sons_[i], UG_NS<dim>::COARSEN_CE))
        return true;
    }
  }

  // neither the element nor a sibling is marked
  return false;
}

//*****************************************************************8
// count
template <int dim, class GridImp>
template <int cc>
int Dune::UGGridEntity<0,dim,GridImp>::count() const
{
  if (dim==3) {

    switch (cc) {
    case 0 :
      return 1;
    case 1 :
      return UG_NS<dim>::Sides_Of_Elem(target_);
    case 2 :
      return UG_NS<dim>::Edges_Of_Elem(target_);
    case 3 :
      return UG_NS<dim>::Corners_Of_Elem(target_);
    }

  } else {

    switch (cc) {
    case 0 :
      return 1;
    case 1 :
      return UG_NS<dim>::Edges_Of_Elem(target_);
    case 2 :
      return UG_NS<dim>::Corners_Of_Elem(target_);
    default :
      // do nothing for cc = 3
      break;
    }

  }
  DUNE_THROW(GridError, "You can't call UGGridEntity<0,dim>::count<codim> "
             << "with dim==" << dim << " and codim==" << cc << "!");
}

template <int dim, class GridImp>
template <int cc>
typename GridImp::template Codim<cc>::Entity
Dune::UGGridEntity<0,dim,GridImp>::subEntity ( int i ) const
{
  assert(i>=0 && i<count<cc>());

  typename UG_NS<dim>::template Entity<cc>::T* subEntity;

  if (cc==dim)
    subEntity = reinterpret_cast<typename UG_NS<dim>::template Entity<cc>::T*>( UG_NS<dim>::Corner(target_,UGGridRenumberer<dim>::verticesDUNEtoUG(i, this->type())) );
  else if (cc==dim - 1)
    subEntity = reinterpret_cast<typename UG_NS<dim>::template Entity<cc>::T*>( UG_NS<dim>::ElementEdge(target_,UGGridRenumberer<dim>::edgesDUNEtoUG(i, this->type())) );
  else if (cc==0)
    subEntity = reinterpret_cast<typename UG_NS<dim>::template Entity<cc>::T*>( target_ );
  else if (dim==3 and cc==1)
    subEntity = reinterpret_cast<typename UG_NS<dim>::template Entity<cc>::T*>( UG_NS<dim>::SideVector(target_,UGGridRenumberer<dim>::facesDUNEtoUG(i, this->type())) );
  else
    DUNE_THROW(GridError, "subEntity for nonexisting codimension requested!" );

  return typename GridImp::template Codim<cc>::Entity(Dune::UGGridEntity<cc,dim,GridImp>((typename UG_NS<dim>::template Entity<cc>::T*)subEntity, gridImp_));
}

template<int dim, class GridImp>
Dune::GeometryType Dune::UGGridEntity<0,dim,GridImp>::type() const
{
  if (dim==2) {

    switch (UG_NS<dim>::Tag(target_)) {
    case UG::D2::TRIANGLE :
      return GeometryType(GeometryType::simplex,2);
    case UG::D2::QUADRILATERAL :
      return GeometryType(GeometryType::cube,2);
    default :
      DUNE_THROW(GridError, "UGGridGeometry::type():  ERROR:  Unknown type "
                 << UG_NS<dim>::Tag(target_) << " found!");
    }

  } else {    // dim==3

    switch (UG_NS<dim>::Tag(target_)) {

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
                 << UG_NS<dim>::Tag(target_) << " found!");

    }

  }

}


template<int dim, class GridImp>
void Dune::UGGridEntity < 0, dim ,GridImp >::
setToTarget(typename UG_NS<dim>::Element* target, const GridImp* gridImp)
{
  target_ = target;
  geo_.setToTarget(target);
  gridImp_ = gridImp;
}

template<int dim, class GridImp>
Dune::UGGridHierarchicIterator<GridImp>
Dune::UGGridEntity < 0, dim ,GridImp >::hbegin(int maxlevel) const
{
  UGGridHierarchicIterator<GridImp> it(maxlevel,gridImp_);

  if (level()<maxlevel) {

    typename UG_NS<dim>::Element* sonList[UG_NS<dim>::MAX_SONS];
    UG_NS<dim>::GetSons(target_, sonList);

    // Load sons of old target onto the iterator stack
    for (int i=0; i<UG_NS<dim>::nSons(target_); i++)
      it.elementStack_.push(sonList[i]);

    it.virtualEntity_.setToTarget( (it.elementStack_.empty())
                                   ? NULL
                                   : it.elementStack_.top(),
                                   gridImp_
                                   );

  } else {
    it.virtualEntity_.setToTarget(nullptr,nullptr);
  }

  return it;
}


template< int dim, class GridImp>
Dune::UGGridHierarchicIterator<GridImp>
Dune::UGGridEntity < 0, dim ,GridImp >::hend(int maxlevel) const
{
  return UGGridHierarchicIterator<GridImp>(maxlevel,gridImp_);
}

template<int dim, class GridImp>
typename Dune::UGGridEntity<0,dim,GridImp>::LocalGeometry
Dune::UGGridEntity < 0, dim, GridImp>::geometryInFather () const
{
  // we need to have a father element
  typename UG_NS<dim>::Element* fatherElement = UG_NS<dim>::EFather(target_);
  if (!fatherElement)
    DUNE_THROW(GridError, "Called geometryInFather() for an entity which doesn't have a father!");

  // The task is to find out the positions of the vertices of this element
  // in the local coordinate system of the father.

  // Get the 'context' of the father element.  In UG-speak, the context is
  // the set of all nodes of an element's sons.  They appear in a fixed
  // order, therefore we can infer the local positions in the father element.
  const int MAX_CORNERS_OF_ELEM = 8;  // this is too much in 2d, but UG is that way
  const int MAX_NEW_CORNERS_DIM = (dim==2) ? 5 : 19;
  const typename UG_NS<dim>::Node* context[MAX_CORNERS_OF_ELEM + MAX_NEW_CORNERS_DIM];
  UG_NS<dim>::GetNodeContext(fatherElement, context);

  // loop through all corner nodes
  std::vector<FieldVector<typename GridImp::ctype,dim> > cornerCoordinates(UG_NS<dim>::Corners_Of_Elem(target_));

  for (int i=0; i<UG_NS<dim>::Corners_Of_Elem(target_); i++) {

    // get corner node pointer
    const typename UG_NS<dim>::Node* fnode = UG_NS<dim>::Corner(target_,i);

    // Find out where in the father's context this node is
    int idx = -1;
    for (int j=0; j<MAX_CORNERS_OF_ELEM + MAX_NEW_CORNERS_DIM; j++)
      if (context[j] == fnode) {
        idx = j;
        break;
      }

    // Node has not been found.  There must be a programming error somewhere
    assert(idx!=-1);

    int duneIdx = UGGridRenumberer<dim>::verticesUGtoDUNE(i,type());

    if (dim==2) {
      switch (UG_NS<dim>::Tag(fatherElement)) {

      case UG::D2::TRIANGLE : {

        assert(idx<6);
        const double coords[6][2] = {
          // The corners
          {0,0}, {1,0}, {0,1},
          // The edge midpoints
          {0.5,0}, {0.5,0.5}, {0,0.5}
        };

        for (int j=0; j<dim; j++)
          cornerCoordinates[duneIdx][j] = coords[idx][j];
        break;
      }
      case UG::D2::QUADRILATERAL : {

        assert(idx<9);
        const double coords[9][2] = {
          // The corners
          {0,0}, {1,0}, {1,1}, {0,1},
          // The edge midpoints
          {0.5,0}, {1,0.5}, {0.5,1}, {0,0.5},
          // The element midpoint
          {0.5,0.5}
        };

        for (int j=0; j<dim; j++)
          cornerCoordinates[duneIdx][j] = coords[idx][j];
        break;
      }

      }

    } else {
      switch (UG_NS<dim>::Tag(fatherElement)) {
      case UG::D3::TETRAHEDRON : {

        // If this assert fails a refinement rule has been appeared which inserts
        // side midpoints.  These have to be added then.
        assert(idx!=10 && idx!=11 && idx!=12 && idx!=13 && idx<15);
        const double coords[15][3] = {
          // The corners
          {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1},
          // The edge midpoints
          {0.5,0,0}, {0.5,0.5,0}, {0,0.5,0},
          {0,0,0.5}, {0.5,0,0.5}, {0,0.5,0.5},
          // Four side midpoints.  I don't know whether they can actually appear
          {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
          // Element midpoints.  Appears in some closure refinement rules
          {0.25,0.25,0.25}
        };

        for (int j=0; j<dim; j++)
          cornerCoordinates[duneIdx][j] = coords[idx][j];
        break;
      }
      case UG::D3::PYRAMID : {

        // If this assert fails a refinement rule has been appeared which inserts
        // side midpoints.  These have to be added then.
        assert( idx<14 || idx==23);
        const double coords[24][3] = {
          // The corners
          {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1},
          // The edge midpoints
          {0.5,0,0}, {1,0.5,0}, {0.5,1,0}, {0,0.5,0},
          {0,0,0.5}, {0.5,0,0.5}, {0.5,0.5,0.5}, {0,0.5,0.5},
          // The bottom face midpoint
          {0.5,0.5,0},
          // Side midpoints of the four triangular faces
          {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
          // Padding due to suboptimal implementation in UG
          {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
          // Element midpoint
          {0.4,0.4,0.2}
        };

        for (int j=0; j<dim; j++)
          cornerCoordinates[duneIdx][j] = coords[idx][j];
        break;
      }
      case UG::D3::PRISM : {
        // If this assert fails a refinement rule has been appeared which inserts
        // side midpoints.  These have to be added then.
        assert(idx!=15 && !(idx>=19 && idx<24) && idx<25);

        const double coords[25][3] = {
          // The corners
          {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,0,1}, {0,1,1},
          // The edge midpoints
          {0.5,0,0}, {0.5,0.5,0}, {0,0.5,0},
          {0,0,0.5}, {1,0,0.5}, {0,1,0.5},
          {0.5,0,1}, {0.5,0.5,1}, {0,0.5,1},
          // A dummy for a midpoint of one of the two triangular faces,
          // which doesn't exist, because triangular faces have no midnodes
          {0,0,0},
          // The midnodes of the three quadrilateral faces
          {0.5,0,0.5}, {0.5,0.5,0.5}, {0,0.5,0.5},
          // Second triangular face
          {0,0,0},
          // Padding due to suboptimal implementation in UG
          {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
          // Side midpoint for the second triangular face
          {0.333333333333333333, 0.333333333333333333, 0.5}
        };

        for (int j=0; j<dim; j++)
          cornerCoordinates[duneIdx][j] = coords[idx][j];
        break;
      }
      case UG::D3::HEXAHEDRON : {

        assert(idx<27);
        const double coords[27][3] = {
          // The corners
          {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1},
          // The edge midpoints
          {0.5,0,0}, {1,0.5,0}, {0.5,1,0}, {0,0.5,0},
          {0,0,0.5}, {1,0,0.5}, {1,1,0.5}, {0,1,0.5},
          {0.5,0,1}, {1,0.5,1}, {0.5,1,1}, {0,0.5,1},
          // The face midpoints
          {0.5,0.5,0}, {0.5,0,0.5}, {1,0.5,0.5}, {0.5,1,0.5}, {0,0.5,0.5}, {0.5,0.5,1},
          // The element midpoint
          {0.5,0.5,0.5}
        };

        for (int j=0; j<dim; j++)
          cornerCoordinates[duneIdx][j] = coords[idx][j];
        break;
      }
      }
    }

  }

  return LocalGeometry(LocalGeometryImpl(type(), cornerCoordinates));
}

template class Dune::UGGridEntity<2,2, const Dune::UGGrid<2> >;
template class Dune::UGGridEntity<3,3, const Dune::UGGrid<3> >;

template class Dune::UGGridEntity<0,2, const Dune::UGGrid<2> >;
template class Dune::UGGridEntity<0,3, const Dune::UGGrid<3> >;

template class Dune::UGGridEntity<1,2, const Dune::UGGrid<2> >;
template class Dune::UGGridEntity<2,3, const Dune::UGGrid<3> >;

template int Dune::UGGridEntity<0, 2, const Dune::UGGrid<2> >::count<0>() const;
template int Dune::UGGridEntity<0, 2, const Dune::UGGrid<2> >::count<1>() const;
template int Dune::UGGridEntity<0, 2, const Dune::UGGrid<2> >::count<2>() const;

template int Dune::UGGridEntity<0, 3, const Dune::UGGrid<3> >::count<0>() const;
template int Dune::UGGridEntity<0, 3, const Dune::UGGrid<3> >::count<1>() const;
template int Dune::UGGridEntity<0, 3, const Dune::UGGrid<3> >::count<2>() const;
template int Dune::UGGridEntity<0, 3, const Dune::UGGrid<3> >::count<3>() const;


template Dune::Grid<2, 2, double, Dune::UGGridFamily<2> >::Codim<0>::Entity
Dune::UGGridEntity<0, 2, const Dune::UGGrid<2> >::subEntity<0>(int) const;

template Dune::Grid<2, 2, double, Dune::UGGridFamily<2> >::Codim<1>::Entity
Dune::UGGridEntity<0, 2, const Dune::UGGrid<2> >::subEntity<1>(int) const;

template Dune::Grid<2, 2, double, Dune::UGGridFamily<2> >::Codim<2>::Entity
Dune::UGGridEntity<0, 2, const Dune::UGGrid<2> >::subEntity<2>(int) const;

template Dune::Grid<3, 3, double, Dune::UGGridFamily<3> >::Codim<0>::Entity
Dune::UGGridEntity<0, 3, const Dune::UGGrid<3> >::subEntity<0>(int) const;

template Dune::Grid<3, 3, double, Dune::UGGridFamily<3> >::Codim<1>::Entity
Dune::UGGridEntity<0, 3, const Dune::UGGrid<3> >::subEntity<1>(int) const;

template Dune::Grid<3, 3, double, Dune::UGGridFamily<3> >::Codim<2>::Entity
Dune::UGGridEntity<0, 3, const Dune::UGGrid<3> >::subEntity<2>(int) const;

template Dune::Grid<3, 3, double, Dune::UGGridFamily<3> >::Codim<3>::Entity
Dune::UGGridEntity<0, 3, const Dune::UGGrid<3> >::subEntity<3>(int) const;
