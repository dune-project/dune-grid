// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridindexsets.hh>

template <class GridImp>
void Dune::UGGridLevelIndexSet<GridImp>::update(const GridImp& grid, int level, std::vector<unsigned int>* nodePermutation) {

  // Commit the index set to a specific level of a specific grid
  grid_ = &grid;
  level_ = level;

  // ///////////////////////////////////
  //   clear index for codim dim-1 and 1
  // ///////////////////////////////////

  typename GridImp::Traits::template Codim<0>::LevelIterator eIt    = grid_->template lbegin<0>(level_);
  typename GridImp::Traits::template Codim<0>::LevelIterator eEndIt = grid_->template lend<0>(level_);

  for (; eIt!=eEndIt; ++eIt) {
    typename UG_NS<dim>::Element* target_ = grid_->getRealImplementation(*eIt).target_;
    // codim dim-1
    for (int i=0; i<eIt->template count<dim-1>(); i++)
    {
      GeometryType gt = eIt->type();
      int a=GenericReferenceElements<double,dim>::general(gt).subEntity(i,dim-1,0,dim);
      int b=GenericReferenceElements<double,dim>::general(gt).subEntity(i,dim-1,1,dim);
      int& index = UG_NS<dim>::levelIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(target_,
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUGNew(a,gt)),
                                                              UG_NS<dim>::Corner(target_,
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUGNew(b,gt))));
      index = -1;
    }
    /** \todo codim 1 (faces) */
    if (dim==3)
      for (int i=0; i<eIt->template count<1>(); i++)
        UG_NS<dim>::levelIndex(UG_NS<dim>::SideVector(target_,i)) = -1;

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

  eIt    = grid_->template lbegin<0>(level_);
  eEndIt = grid_->template lend<0>(level_);

  for (; eIt!=eEndIt; ++eIt) {

    typename UG_NS<dim>::Element* target = grid_->getRealImplementation(*eIt).target_;

    // codim 0 (elements)
    GeometryType eType = eIt->type();
    if (eType.isSimplex()) {
      UG_NS<dim>::levelIndex(target) = numSimplices_++;
    } else if (eType.isPyramid()) {
      UG_NS<dim>::levelIndex(target) = numPyramids_++;
    } else if (eType.isPrism()) {
      UG_NS<dim>::levelIndex(target) = numPrisms_++;
    } else if (eType.isCube()) {
      UG_NS<dim>::levelIndex(target) = numCubes_++;
    } else {
      DUNE_THROW(GridError, "Found the GeometryType " << eIt->type()
                                                      << ", which should never occur in a UGGrid!");
    }

    // codim dim-1 (edges)
    for (int i=0; i<eIt->template count<dim-1>(); i++)
    {
      int a=GenericReferenceElements<double,dim>::general(eType).subEntity(i,dim-1,0,dim);
      int b=GenericReferenceElements<double,dim>::general(eType).subEntity(i,dim-1,1,dim);
      int& index = UG_NS<dim>::levelIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(target,
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUGNew(a,eType)),
                                                              UG_NS<dim>::Corner(target,
                                                                                 UGGridRenumberer<dim>::verticesDUNEtoUGNew(b,eType))));
      if (index<0)
        index = numEdges_++;
    }

    // codim 1 (faces): todo
    if (dim==3)
      for (int i=0; i<eIt->template count<1>(); i++)
      {
        int& index = UG_NS<dim>::levelIndex(UG_NS<dim>::SideVector(target,UGGridRenumberer<dim>::facesDUNEtoUGNew(i,eType)));
        if (index<0) {             // not visited yet
          GeometryType gtType = GenericReferenceElements<double,dim>::general(eType).type(i,1);
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
    myTypes_[0].push_back(GeometryType(GeometryType::simplex,dim));
  if (numPyramids_ > 0)
    myTypes_[0].push_back(GeometryType(GeometryType::pyramid,dim));
  if (numPrisms_ > 0)
    myTypes_[0].push_back(GeometryType(GeometryType::prism,dim));
  if (numCubes_ > 0)
    myTypes_[0].push_back(GeometryType(GeometryType::cube,dim));

  myTypes_[dim-1].resize(0);
  myTypes_[dim-1].push_back(GeometryType(1));

  if (dim==3) {
    myTypes_[1].resize(0);
    if (numTriFaces_ > 0)
      myTypes_[1].push_back(GeometryType(GeometryType::simplex,dim-1));
    if (numQuadFaces_ > 0)
      myTypes_[1].push_back(GeometryType(GeometryType::cube,dim-1));
  }

  // //////////////////////////////
  //   Init the vertex indices
  // //////////////////////////////

  typename GridImp::Traits::template Codim<dim>::LevelIterator vIt    = grid_->template lbegin<dim>(level_);
  typename GridImp::Traits::template Codim<dim>::LevelIterator vEndIt = grid_->template lend<dim>(level_);

  numVertices_ = 0;

  if (nodePermutation!=0 and level==0)
  {
    for (; vIt!=vEndIt; ++vIt)
      UG_NS<dim>::levelIndex(grid_->getRealImplementation(*vIt).target_) = (*nodePermutation)[numVertices_++];
  }
  else
  {
    for (; vIt!=vEndIt; ++vIt)
      UG_NS<dim>::levelIndex(grid_->getRealImplementation(*vIt).target_) = numVertices_++;
  }

  /*    for (; vIt!=vEndIt; ++vIt)
          UG_NS<dim>::levelIndex(grid_->getRealImplementation(*vIt).target_) = numVertices_++;*/

  myTypes_[dim].resize(0);
  myTypes_[dim].push_back(GeometryType(GeometryType::cube,0));
}

template <class GridImp>
void Dune::UGGridLeafIndexSet<GridImp>::update(std::vector<unsigned int>* nodePermutation) {

  // //////////////////////////////////////////////////////
  // Handle codim 1 and dim-1: levelwise from top to bottom
  // //////////////////////////////////////////////////////

  // first loop : clear indices
  for (int level_=grid_.maxLevel(); level_>=0; level_--)
  {
    typename GridImp::Traits::template Codim<0>::LevelIterator eIt    = grid_.template lbegin<0>(level_);
    typename GridImp::Traits::template Codim<0>::LevelIterator eEndIt = grid_.template lend<0>(level_);

    for (; eIt!=eEndIt; ++eIt)
    {
      // get pointer to UG object
      typename UG_NS<dim>::Element* target_ = grid_.getRealImplementation(*eIt).target_;

      // codim dim-1
      for (int i=0; i<eIt->template count<dim-1>(); i++)
      {
        GeometryType gt = eIt->type();
        int a=GenericReferenceElements<double,dim>::general(gt).subEntity(i,dim-1,0,dim);
        int b=GenericReferenceElements<double,dim>::general(gt).subEntity(i,dim-1,1,dim);
        int& index = UG_NS<dim>::leafIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(target_,
                                                                                  UGGridRenumberer<dim>::verticesDUNEtoUGNew(a,gt)),
                                                               UG_NS<dim>::Corner(target_,
                                                                                  UGGridRenumberer<dim>::verticesDUNEtoUGNew(b,gt))));
        index = -1;
      }

      // codim 1 (faces): todo
      if (dim==3)
        for (int i=0; i<eIt->template count<1>(); i++)
        {
          GeometryType gt = eIt->type();
          int& index = UG_NS<dim>::leafIndex(UG_NS<dim>::SideVector(target_,UGGridRenumberer<dim>::facesDUNEtoUGNew(i,gt)));
          index = -1;
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

    typename GridImp::Traits::template Codim<0>::LevelIterator eIt    = grid_.template lbegin<0>(level_);
    typename GridImp::Traits::template Codim<0>::LevelIterator eEndIt = grid_.template lend<0>(level_);

    for (; eIt!=eEndIt; ++eIt)
    {
      // we need only look at leaf elements
      if (!eIt->isLeaf())
        continue;
      else
        containsLeafElements = true;

      // get pointer to UG object
      typename UG_NS<dim>::Element* target_ = grid_.getRealImplementation(*eIt).target_;

      // codim dim-1 (edges)
      for (int i=0; i<eIt->template count<dim-1>(); i++)
      {
        GeometryType gt = eIt->type();
        int a=GenericReferenceElements<double,dim>::general(gt).subEntity(i,dim-1,0,dim);
        int b=GenericReferenceElements<double,dim>::general(gt).subEntity(i,dim-1,1,dim);
        int& index = UG_NS<dim>::leafIndex(UG_NS<dim>::GetEdge(UG_NS<dim>::Corner(target_,
                                                                                  UGGridRenumberer<dim>::verticesDUNEtoUGNew(a,gt)),
                                                               UG_NS<dim>::Corner(target_,
                                                                                  UGGridRenumberer<dim>::verticesDUNEtoUGNew(b,gt))));
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
                                                                         UGGridRenumberer<dim>::verticesDUNEtoUGNew(a,gt)),
                                                      UG_NS<dim>::Corner(father_,
                                                                         UGGridRenumberer<dim>::verticesDUNEtoUGNew(b,gt)))) = index;
            father_ = UG_NS<dim>::EFather(father_);
          }
        }
      }

      // codim 1 (faces): todo
      if (dim==3)
        for (int i=0; i<eIt->template count<1>(); i++)
        {
          GeometryType gt = eIt->type();
          int& index = UG_NS<dim>::leafIndex(UG_NS<dim>::SideVector(target_,UGGridRenumberer<dim>::facesDUNEtoUGNew(i,gt)));
          if (index<0)                       // not visited yet
          {
            // get new index and assign
            GeometryType gtType = GenericReferenceElements<double,dim>::general(gt).type(i,1);
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
              UG_NS<dim>::leafIndex(UG_NS<dim>::SideVector(father_,UGGridRenumberer<dim>::facesDUNEtoUGNew(i,gt))) = index;
              father_ = UG_NS<dim>::EFather(father_);
            }
          }
        }
    }

    if (containsLeafElements)
      coarsestLevelWithLeafElements_ = level_;

  }

  // Update the list of geometry types present
  myTypes_[dim-1].resize(0);
  myTypes_[dim-1].push_back(GeometryType(GeometryType::cube,1));

  if (dim==3) {

    myTypes_[1].resize(0);
    if (numTriFaces_ > 0)
      myTypes_[1].push_back(GeometryType(GeometryType::simplex,dim-1));
    if (numQuadFaces_ > 0)
      myTypes_[1].push_back(GeometryType(GeometryType::cube,dim-1));

  }

  // ///////////////////////////////
  //   Init the element indices
  // ///////////////////////////////
  numSimplices_ = 0;
  numPyramids_  = 0;
  numPrisms_    = 0;
  numCubes_     = 0;

  typename GridImp::Traits::template Codim<0>::LeafIterator eIt    = grid_.template leafbegin<0>();
  typename GridImp::Traits::template Codim<0>::LeafIterator eEndIt = grid_.template leafend<0>();

  for (; eIt!=eEndIt; ++eIt) {

    GeometryType eType = eIt->type();

    if (eType.isSimplex())
      UG_NS<dim>::leafIndex(grid_.getRealImplementation(*eIt).target_) = numSimplices_++;
    else if (eType.isPyramid())
      UG_NS<dim>::leafIndex(grid_.getRealImplementation(*eIt).target_) = numPyramids_++;
    else if (eType.isPrism())
      UG_NS<dim>::leafIndex(grid_.getRealImplementation(*eIt).target_) = numPrisms_++;
    else if (eType.isCube())
      UG_NS<dim>::leafIndex(grid_.getRealImplementation(*eIt).target_) = numCubes_++;
    else {
      DUNE_THROW(GridError, "Found the GeometryType " << eType
                                                      << ", which should never occur in a UGGrid!");
    }
  }

  // Update the list of geometry types present
  myTypes_[0].resize(0);
  if (numSimplices_ > 0)
    myTypes_[0].push_back(GeometryType(GeometryType::simplex,dim));
  if (numPyramids_ > 0)
    myTypes_[0].push_back(GeometryType(GeometryType::pyramid,dim));
  if (numPrisms_ > 0)
    myTypes_[0].push_back(GeometryType(GeometryType::prism,dim));
  if (numCubes_ > 0)
    myTypes_[0].push_back(GeometryType(GeometryType::cube,dim));


  // //////////////////////////////
  //   Init the vertex indices
  // //////////////////////////////
  typename GridImp::Traits::template Codim<dim>::LeafIterator vIt    = grid_.template leafbegin<dim>();
  typename GridImp::Traits::template Codim<dim>::LeafIterator vEndIt = grid_.template leafend<dim>();

  // leaf index in node writes through to vertex !
  numVertices_ = 0;

  if (nodePermutation!=0 and grid_.maxLevel()==0)
  {
    for (; vIt!=vEndIt; ++vIt)
      UG_NS<dim>::leafIndex(grid_.getRealImplementation(*vIt).target_) = (*nodePermutation)[numVertices_++];
  }
  else
  {
    for (; vIt!=vEndIt; ++vIt)
      UG_NS<dim>::leafIndex(grid_.getRealImplementation(*vIt).target_) = numVertices_++;
  }

  myTypes_[dim].resize(0);
  myTypes_[dim].push_back(GeometryType(0));

}

// Explicit template instantiations to compile the stuff in this file

template class Dune::UGGridLevelIndexSet<const Dune::UGGrid<2> >;
template class Dune::UGGridLevelIndexSet<const Dune::UGGrid<3> >;

template class Dune::UGGridLeafIndexSet<const Dune::UGGrid<2> >;
template class Dune::UGGridLeafIndexSet<const Dune::UGGrid<3> >;
