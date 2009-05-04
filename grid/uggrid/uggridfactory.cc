// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/uggrid/uggridfactory.hh>

template <int dimworld>
Dune::GridFactory<Dune::UGGrid<dimworld> >::
GridFactory()
{
  grid_ = new Dune::UGGrid<dimworld>;

  factoryOwnsGrid_ = true;

  createBegin();
}

template <int dimworld>
Dune::GridFactory<Dune::UGGrid<dimworld> >::
GridFactory(UGGrid<dimworld>* grid)
{
  grid_ = grid;

  factoryOwnsGrid_ = false;

  createBegin();
}

template <int dimworld>
Dune::GridFactory<Dune::UGGrid<dimworld> >::
~GridFactory()
{
  if (grid_ && factoryOwnsGrid_)
    delete grid_;
}

template <int dimworld>
void Dune::GridFactory<Dune::UGGrid<dimworld> >::
insertVertex(const Dune::FieldVector<typename Dune::GridFactory<Dune::UGGrid<dimworld> >::ctype,dimworld>& pos)
{
  grid_->vertexPositions_.push_back(pos);
}

template <int dimworld>
void Dune::GridFactory<Dune::UGGrid<dimworld> >::
insertElement(const GeometryType& type,
              const std::vector<unsigned int>& vertices)
{
  grid_->insertElement(type, vertices);
}

template <int dimworld>
void Dune::GridFactory<Dune::UGGrid<dimworld> >::
insertBoundarySegment(const std::vector<unsigned int> vertices,
                      const BoundarySegment<dimworld>* boundarySegment)
{
  grid_->insertBoundarySegment(vertices, boundarySegment);
}

template <int dimworld>
Dune::UGGrid<dimworld>* Dune::GridFactory<Dune::UGGrid<dimworld> >::
createGrid()
{
  // Prevent a crash when this method is called twice in a row
  // You never know who may do this...
  if (grid_==NULL)
    return NULL;

  // finalize grid creation
  grid_->createEnd();

  // hand it over and delete the member pointer
  Dune::UGGrid<dimworld>* tmp = grid_;
  grid_ = NULL;
  return tmp;
}

template <int dimworld>
void Dune::GridFactory<Dune::UGGrid<dimworld> >::
createBegin()
{
  // ///////////////////////////////////////////////////////
  //   Clean up existing grid structure if there is one
  // ///////////////////////////////////////////////////////
  for (unsigned int i=0; i<grid_->boundarySegments_.size(); i++)
    delete grid_->boundarySegments_[i];

  // Delete the UG multigrid if there is one (== createEnd() has already
  // been called once for this object)
  if (grid_->multigrid_) {
    // Set UG's currBVP variable to the BVP corresponding to this
    // grid.  This is necessary if we have more than one UGGrid in use.
    // DisposeMultiGrid will crash if we don't do this
    //UG_NS<dim>::Set_Current_BVP(grid_->multigrid_->theBVP);
    // set the multigrid's bvp pointer to NULL to make sure the BVP
    // is not deleted
    grid_->multigrid_->theBVP = NULL;
    UG_NS<dimworld>::DisposeMultiGrid(grid_->multigrid_);
    grid_->multigrid_ = NULL;
  }

  // Delete levelIndexSets if there are any
  for (unsigned int i=0; i<grid_->levelIndexSets_.size(); i++)
    if (grid_->levelIndexSets_[i])
      delete grid_->levelIndexSets_[i];

  grid_->levelIndexSets_.resize(0);

  // //////////////////////////////////////////////////////////
  //   Clear all buffers used during coarse grid creation
  // //////////////////////////////////////////////////////////
  grid_->boundarySegments_.resize(0);
  grid_->boundarySegmentVertices_.resize(0);
  grid_->elementTypes_.resize(0);
  grid_->elementVertices_.resize(0);
  grid_->vertexPositions_.resize(0);

  // //////////////////////////////////////////////////////////
  //   Delete the UG domain, if it exists
  // //////////////////////////////////////////////////////////
  std::string domainName = grid_->name_ + "_Domain";
  UG_NS<dimworld>::RemoveDomain(domainName.c_str());
}




// Explicit template instatiation
template class Dune::GridFactory<Dune::UGGrid<2> >;
template class Dune::GridFactory<Dune::UGGrid<3> >;
