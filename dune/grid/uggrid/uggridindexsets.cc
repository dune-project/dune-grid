// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridindexsets.hh>

namespace Dune {

template <class GridImp>
void UGGridLevelIndexSet<GridImp>::update(const GridImp& grid, int level, std::vector<unsigned int>* nodePermutation) {

  // Commit the index set to a specific level of a specific grid
  grid_ = &grid;
  level_ = level;

  // ///////////////////////////////////
  //   clear index for codim dim-1 and 1
  // ///////////////////////////////////

  for (const auto& element : elements(grid.levelGridView(level_))) {
    typename UG_NS<dim>::Element* target_ = element.impl().target_;
    // codim dim-1
    for (unsigned int i=0; i<element.subEntities(dim-1); i++)
    {
      GeometryType gt = element.type();
      auto ref_el = referenceElement<double,dim>(gt);
      int a = ref_el.subEntity(i,dim-1,0,dim);
      int b = ref_el.subEntity(i,dim-1,1,dim);
      int& index = UG_NS<dim>::levelIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(target_,
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(a,gt)),
                                                              UG_NS<dim>::Corner(target_,
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(b,gt))));
      index = -1;
    }
    /** \todo codim 1 (faces) */
    if (dim==3)
      for (unsigned int i=0; i<element.subEntities(1); i++)
        UG_NS<dim>::levelIndex(UG_NS<dim>::SideVector(target_,i)) = std::numeric_limits<UG::UINT>::max();

  }

  // ///////////////////////////////
  //   Init the codim<dim indices
  // ///////////////////////////////
  numSimplices_ = 0;
  numPyramids_  = 0;
  numPrisms_    = 0;
  numCubes_     = 0;
  numEdges_     = 0;
  numTriFaces_  = 0;
  numQuadFaces_ = 0;

  for (const auto& element : elements(grid.levelGridView(level_))) {

    typename UG_NS<dim>::Element* target = element.impl().target_;

    // codim 0 (elements)
    GeometryType eType = element.type();
    if (eType.isSimplex()) {
      UG_NS<dim>::levelIndex(target) = numSimplices_++;
    } else if (eType.isPyramid()) {
      UG_NS<dim>::levelIndex(target) = numPyramids_++;
    } else if (eType.isPrism()) {
      UG_NS<dim>::levelIndex(target) = numPrisms_++;
    } else if (eType.isCube()) {
      UG_NS<dim>::levelIndex(target) = numCubes_++;
    } else {
      DUNE_THROW(GridError, "Found the GeometryType " << element.type()
                                                      << ", which should never occur in a UGGrid!");
    }

    // codim dim-1 (edges)
    for (unsigned int i=0; i<element.subEntities(dim-1); i++)
    {
      auto ref_el = referenceElement<double,dim>(eType);
      int a = ref_el.subEntity(i,dim-1,0,dim);
      int b = ref_el.subEntity(i,dim-1,1,dim);
      int& index = UG_NS<dim>::levelIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(target,
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(a,eType)),
                                                              UG_NS<dim>::Corner(target,
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUG(b,eType))));
      if (index<0)
        index = numEdges_++;
    }

    // codim 1 (faces): todo
    if (dim==3)
      for (unsigned int i=0; i<element.subEntities(1); i++)
      {
        UG::UINT& index = UG_NS<dim>::levelIndex(UG_NS<dim>::SideVector(target,UGGridRenumberer<dim>::facesDUNEtoUG(i,eType)));
        if (index == std::numeric_limits<UG::UINT>::max()) {             // not visited yet
          GeometryType gtType = referenceElement<double,dim>(eType).type(i,1);
          if (gtType.isSimplex()) {
            index = numTriFaces_++;
          } else if (gtType.isCube()) {
            index = numQuadFaces_++;
          } else {
            std::cout << "face geometry type is " << gtType << std::endl;
            DUNE_THROW(GridError, "wrong geometry type in face");
          }
        }
      }
  }

  // Update the list of geometry types present
  myTypes_[0].resize(0);
  if (numSimplices_ > 0)
    myTypes_[0].push_back(GeometryTypes::simplex(dim));
  if (numPyramids_ > 0)
    myTypes_[0].push_back(GeometryTypes::pyramid);
  if (numPrisms_ > 0)
    myTypes_[0].push_back(GeometryTypes::prism);
  if (numCubes_ > 0)
    myTypes_[0].push_back(GeometryTypes::cube(dim));

  myTypes_[dim-1].resize(0);
  myTypes_[dim-1].push_back(GeometryTypes::line);

  if (dim==3) {
    myTypes_[1].resize(0);
    if (numTriFaces_ > 0)
      myTypes_[1].push_back(GeometryTypes::triangle);
    if (numQuadFaces_ > 0)
      myTypes_[1].push_back(GeometryTypes::quadrilateral);
  }

  // //////////////////////////////
  //   Init the vertex indices
  // //////////////////////////////

  numVertices_ = 0;

  if (nodePermutation!=0 and level==0)
  {
    for (const auto& vertex : vertices(grid.levelGridView(level_)))
      UG_NS<dim>::levelIndex(vertex.impl().target_) = (*nodePermutation)[numVertices_++];
  }
  else
  {
    for (const auto& vertex : vertices(grid.levelGridView(level_)))
      UG_NS<dim>::levelIndex(vertex.impl().target_) = numVertices_++;
  }

  /*    for (; vIt!=vEndIt; ++vIt)
          UG_NS<dim>::levelIndex(vIt->impl().target_) = numVertices_++;*/

  myTypes_[dim].resize(0);
  myTypes_[dim].push_back(GeometryTypes::vertex);
}

template <class GridImp>
void UGGridLeafIndexSet<GridImp>::update(std::vector<unsigned int>* nodePermutation) {

  // //////////////////////////////////////////////////////
  // Handle codim 1 and dim-1: levelwise from top to bottom
  // //////////////////////////////////////////////////////

  // first loop : clear indices
  for (int level_=grid_.maxLevel(); level_>=0; level_--)
  {
    for (const auto& element : elements(grid_.levelGridView(level_)))
    {
      // get pointer to UG object
      typename UG_NS<dim>::Element* target_ = element.impl().target_;

      // codim dim-1
      for (unsigned int i=0; i<element.subEntities(dim-1); i++)
      {
        GeometryType gt = element.type();
        auto ref_el = referenceElement<double,dim>(gt);
        int a = ref_el.subEntity(i,dim-1,0,dim);
        int b = ref_el.subEntity(i,dim-1,1,dim);
        int& index = UG_NS<dim>::leafIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(target_,
                                                                                  UGGridRenumberer<dim>::verticesDUNEtoUG(a,gt)),
                                                               UG_NS<dim>::Corner(target_,
                                                                                  UGGridRenumberer<dim>::verticesDUNEtoUG(b,gt))));
        index = -1;
      }

      // codim 1 (faces): todo
      if (dim==3)
        for (unsigned int i=0; i<element.subEntities(1); i++)
        {
          GeometryType gt = element.type();
          UG::UINT& index = UG_NS<dim>::leafIndex(UG_NS<dim>::SideVector(target_,UGGridRenumberer<dim>::facesDUNEtoUG(i,gt)));
          index = std::numeric_limits<UG::UINT>::max();
        }

      // reset the isLeaf information of the nodes
      for (unsigned int i=0; i<element.subEntities(dim); i++)
      {
        // no need to renumber, we are visiting all the corners
        UG_NS<dim>::Corner(target_,i)->isLeaf = false;
      }
    }
  }

  // init counters
  numEdges_     = 0;
  numTriFaces_  = 0;
  numQuadFaces_ = 0;

  // second loop : set indices
  for (int level_=grid_.maxLevel(); level_>=0; level_--)
  {

    // used to compute the coarsest level with leaf elements
    bool containsLeafElements = false;

    // this value is used in the parallel case, when a local grid may not contain any elements at all
    coarsestLevelWithLeafElements_ = 0;

    for (const auto& element : elements(grid_.levelGridView(level_)))
    {
      // we need only look at leaf elements
      if (!element.isLeaf())
        continue;
      else
        containsLeafElements = true;

      // get pointer to UG object
      typename UG_NS<dim>::Element* target_ = element.impl().target_;

      // codim dim-1 (edges)
      for (unsigned int i=0; i<element.subEntities(dim-1); i++)
      {
        GeometryType gt = element.type();
        auto ref_el = referenceElement<double,dim>(gt);
        int a = ref_el.subEntity(i,dim-1,0,dim);
        int b = ref_el.subEntity(i,dim-1,1,dim);
        int& index = UG_NS<dim>::leafIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(target_,
                                                                                  UGGridRenumberer<dim>::verticesDUNEtoUG(a,gt)),
                                                               UG_NS<dim>::Corner(target_,
                                                                                  UGGridRenumberer<dim>::verticesDUNEtoUG(b,gt))));
        if (index<0)
        {
          // get new index and assign
          index = numEdges_++;
          // write index through to coarser grids
          typename UG_NS<dim>::Element* father_ = UG_NS<dim>::EFather(target_);
          while (father_!=0)
          {
            if (!UG_NS<dim>::hasCopy(father_)) break;                                 // handle only copies
            UG_NS<dim>::leafIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(father_,
                                                                         UGGridRenumberer<dim>::verticesDUNEtoUG(a,gt)),
                                                      UG_NS<dim>::Corner(father_,
                                                                         UGGridRenumberer<dim>::verticesDUNEtoUG(b,gt)))) = index;
            father_ = UG_NS<dim>::EFather(father_);
          }
        }
      }

      // codim 1 (faces): todo
      if (dim==3)
        for (unsigned int i=0; i<element.subEntities(1); i++)
        {
          GeometryType gt = element.type();
          UG::UINT& index = UG_NS<dim>::leafIndex(UG_NS<dim>::SideVector(target_,UGGridRenumberer<dim>::facesDUNEtoUG(i,gt)));
          if (index==std::numeric_limits<UG::UINT>::max())                       // not visited yet
          {
            // get new index and assign
            GeometryType gtType = referenceElement<double,dim>(gt).type(i,1);
            if (gtType.isSimplex())
              index = numTriFaces_++;
            else if (gtType.isCube())
              index = numQuadFaces_++;
            else {
              std::cout << "face geometry type is " << gtType << std::endl;
              DUNE_THROW(GridError, "wrong geometry type in face");
            }
            // write index through to coarser grid
            typename UG_NS<dim>::Element* father_ = UG_NS<dim>::EFather(target_);
            while (father_!=0)
            {
              if (!UG_NS<dim>::hasCopy(father_)) break;                                   // handle only copies
              UG_NS<dim>::leafIndex(UG_NS<dim>::SideVector(father_,UGGridRenumberer<dim>::facesDUNEtoUG(i,gt))) = index;
              father_ = UG_NS<dim>::EFather(father_);
            }
          }
        }

      // set the isLeaf information of the nodes based on the leaf elements
      for (unsigned int i=0; i< element.subEntities(dim); i++)
      {
        // no need to renumber, we are visiting all the corners
        typename UG_NS<dim>::Node* theNode = UG_NS<dim>::Corner(target_,i);
        theNode->isLeaf = true;
      }
    }

    if (containsLeafElements)
      coarsestLevelWithLeafElements_ = level_;

  }

  // Update the list of geometry types present
  myTypes_[dim-1].resize(0);
  myTypes_[dim-1].push_back(GeometryTypes::line);

  if (dim==3) {

    myTypes_[1].resize(0);
    if (numTriFaces_ > 0)
      myTypes_[1].push_back(GeometryTypes::triangle);
    if (numQuadFaces_ > 0)
      myTypes_[1].push_back(GeometryTypes::quadrilateral);

  }

  // ///////////////////////////////
  //   Init the element indices
  // ///////////////////////////////
  numSimplices_ = 0;
  numPyramids_  = 0;
  numPrisms_    = 0;
  numCubes_     = 0;

  for (const auto& element : elements(grid_.leafGridView())) {

    GeometryType eType = element.type();

    if (eType.isSimplex())
      UG_NS<dim>::leafIndex(element.impl().target_) = numSimplices_++;
    else if (eType.isPyramid())
      UG_NS<dim>::leafIndex(element.impl().target_) = numPyramids_++;
    else if (eType.isPrism())
      UG_NS<dim>::leafIndex(element.impl().target_) = numPrisms_++;
    else if (eType.isCube())
      UG_NS<dim>::leafIndex(element.impl().target_) = numCubes_++;
    else {
      DUNE_THROW(GridError, "Found the GeometryType " << eType
                                                      << ", which should never occur in a UGGrid!");
    }
  }

  // Update the list of geometry types present
  myTypes_[0].resize(0);
  if (numSimplices_ > 0)
    myTypes_[0].push_back(GeometryTypes::simplex(dim));
  if (numPyramids_ > 0)
    myTypes_[0].push_back(GeometryTypes::pyramid);
  if (numPrisms_ > 0)
    myTypes_[0].push_back(GeometryTypes::prism);
  if (numCubes_ > 0)
    myTypes_[0].push_back(GeometryTypes::cube(dim));


  // //////////////////////////////
  //   Init the vertex indices
  // //////////////////////////////
  // leaf index in node writes through to vertex !
  numVertices_ = 0;

  if (nodePermutation!=0 and grid_.maxLevel()==0)
  {
    for (const auto& vertex : vertices(grid_.leafGridView()))
      UG_NS<dim>::leafIndex(vertex.impl().target_) = (*nodePermutation)[numVertices_++];
  }
  else
  {
    for (const auto& vertex : vertices(grid_.leafGridView()))
      UG_NS<dim>::leafIndex(vertex.impl().target_) = numVertices_++;
  }

  myTypes_[dim].resize(0);
  myTypes_[dim].push_back(GeometryTypes::vertex);

}

// Explicit template instantiations to compile the stuff in this file

template class UGGridLevelIndexSet<const UGGrid<2> >;
template class UGGridLevelIndexSet<const UGGrid<3> >;

template class UGGridLeafIndexSet<const UGGrid<2> >;
template class UGGridLeafIndexSet<const UGGrid<3> >;

} /* namespace Dune */
