// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/onedgrid/onedgridfactory.hh>

using namespace Dune;


Dune::GridFactory<Dune::OneDGrid >::
GridFactory() :
  factoryOwnsGrid_(true),
  vertexIndex_(0)
{
  grid_ = new OneDGrid;

  createBegin();
}

Dune::GridFactory<Dune::OneDGrid >::
GridFactory(OneDGrid* grid) :
  factoryOwnsGrid_(true),
  vertexIndex_(0)
{
  grid_ = grid;

  createBegin();
}

Dune::GridFactory<Dune::OneDGrid>::
~GridFactory()
{
  if (grid_ && factoryOwnsGrid_)
    delete grid_;
}

void Dune::GridFactory<Dune::OneDGrid>::
insertVertex(const Dune::FieldVector<GridFactory<OneDGrid >::ctype,1>& pos)
{
  vertexPositions_.insert(std::make_pair(pos, vertexIndex_++));
}

void Dune::GridFactory<Dune::OneDGrid>::
insertElement(const GeometryType& type,
              const std::vector<unsigned int>& vertices)
{
  if (type.dim() != 1)
    DUNE_THROW(GridError, "You cannot insert a " << type << " into a OneDGrid!");

  if (vertices.size() != 2)
    DUNE_THROW(GridError, "You cannot insert an element with " << vertices.size() << " vertices into a OneDGrid!");

  elements_.push_back(Dune::array<unsigned int,2>());
  elements_.back()[0] = vertices[0];
  elements_.back()[1] = vertices[1];

}

Dune::OneDGrid* Dune::GridFactory<Dune::OneDGrid>::
createGrid()
{
  // Prevent a crash when this method is called twice in a row
  // You never know who may do this...
  if (grid_==NULL)
    return NULL;

  // ////////////////////////////////////////////////////////
  //   Insert the vertices into the grid
  // ////////////////////////////////////////////////////////

  grid_->vertices.resize(1);
  grid_->elements.resize(1);


  VertexIterator vIt    = vertexPositions_.begin();
  VertexIterator vEndIt = vertexPositions_.end();

  for (; vIt!=vEndIt; ++vIt) {

    OneDEntityImp<0> newVertex(0,
                               vIt->first,
                               grid_->getNextFreeId(1));

    newVertex.leafIndex_ = vIt->second;
    newVertex.levelIndex_ = vIt->second;

    grid_->vertices[0].push_back(newVertex);

  }

  // ///////////////////////////////////////////////////////////////////
  //   Insert the elements into the grid
  //
  // This is a 1d grid and currently it has to be connected. Hence we actually
  // know where the elements are, even without being told explicitly.
  // The only thing of interest are the indices.
  // ///////////////////////////////////////////////////////////////////

  OneDGridList<OneDEntityImp<0> >::iterator it = grid_->vertices[0].begin();
  for (size_t i=0; i<vertexPositions_.size()-1; i++) {

    OneDEntityImp<1> newElement(0, grid_->getNextFreeId(0));
    newElement.vertex_[0] = it;
    it = it->succ_;
    newElement.vertex_[1] = it;

    grid_->elements[0].push_back(newElement);

  }

  // ///////////////////////////////////////////////////
  //   Create the index sets
  // ///////////////////////////////////////////////////

  grid_->setIndices();

  // ///////////////////////////////////////////////////
  // hand over the grid and delete the member pointer
  // ///////////////////////////////////////////////////

  Dune::OneDGrid* tmp = grid_;
  grid_ = NULL;
  return tmp;
}

void Dune::GridFactory<Dune::OneDGrid >::
createBegin()
{
  vertexPositions_.clear();

}
