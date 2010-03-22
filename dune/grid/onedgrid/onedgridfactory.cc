// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/onedgrid/onedgridfactory.hh>
#include <dune/grid/onedgrid/onedgridindexsets.hh>

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
  factoryOwnsGrid_(false),
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

void Dune::GridFactory<Dune::OneDGrid>::
insertBoundarySegment(const std::vector<unsigned int>& vertices)
{
  if (vertices.size() != 1)
    DUNE_THROW(GridError, "OneDGrid BoundarySegments must have exactly one vertex.");

  boundarySegments_.push_back(vertices[0]);
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

  grid_->entityImps_.resize(1);


  VertexIterator vIt    = vertexPositions_.begin();
  VertexIterator vEndIt = vertexPositions_.end();

  for (; vIt!=vEndIt; ++vIt) {

    OneDEntityImp<0> newVertex(0,
                               vIt->first,
                               grid_->getNextFreeId(1));

    newVertex.leafIndex_ = vIt->second;
    newVertex.levelIndex_ = vIt->second;

    grid_->vertices(0).push_back(newVertex);

  }

  // //////////////////////////////////////////////////////////////////
  //   Make an array with the vertex positions accessible by index
  //   We'll need that several times.
  // //////////////////////////////////////////////////////////////////
  std::vector<double> vertexPositionsByIndex(vertexPositions_.size());
  for (std::map<FieldVector<ctype,1>, unsigned int >::iterator it = vertexPositions_.begin();
       it != vertexPositions_.end();
       ++it)
    vertexPositionsByIndex[it->second] = it->first;

  // ///////////////////////////////////////////////////
  //   Set the numbering of the boundary segments
  // ///////////////////////////////////////////////////

  if (boundarySegments_.size() > 2)
    DUNE_THROW(GridError, "You cannot provide more than two boundary segments to a OneDGrid (it must be connected).");

  if (boundarySegments_.size() > 1
      && vertexPositionsByIndex[boundarySegments_[0]] > vertexPositions_.begin()->first[0])
    grid_->reversedBoundarySegmentNumbering_ = true;

  // ///////////////////////////////////////////////////////////////////
  //   Insert the elements into the grid
  //
  // This is a 1d grid and it has to be connected. Hence we actually
  // know where the elements are, even without being told explicitly.
  // The only thing of interest are the indices.
  // ///////////////////////////////////////////////////////////////////

  // first sort elements by increasing position.  That is how they are expected in the grid data structure
  std::map<double, std::pair<Dune::array<unsigned int, 2>, unsigned int> > elementsByPosition;

  for (size_t i=0; i<elements_.size(); i++)
    elementsByPosition.insert(std::make_pair(vertexPositionsByIndex[elements_[i][0]],     // order by position of left vertex
                                             std::make_pair(elements_[i], i)      // the element and its position in the insertion sequence
                                             ));


  OneDGridList<OneDEntityImp<0> >::iterator it = Dune::get<0>(grid_->entityImps_[0]).begin();
  std::map<double, std::pair<Dune::array<unsigned int, 2>, unsigned int> >::iterator eIt = elementsByPosition.begin();

  /** \todo Looping over the vertices to get all elements assumes that
      the grid is connected.
   */
  for (size_t i=0; i<vertexPositions_.size()-1; i++, ++eIt) {

    OneDEntityImp<1> newElement(0, grid_->getNextFreeId(0), grid_->reversedBoundarySegmentNumbering_);
    newElement.vertex_[0] = it;
    it = it->succ_;
    newElement.vertex_[1] = it;

    newElement.levelIndex_ = eIt->second.second;
    newElement.leafIndex_  = eIt->second.second;

    grid_->elements(0).push_back(newElement);

  }

  // ///////////////////////////////////////////////////
  //   Create the index sets
  // ///////////////////////////////////////////////////

  grid_->levelIndexSets_.resize(1);
  grid_->levelIndexSets_[0] = new OneDGridLevelIndexSet<const OneDGrid>(*grid_, 0);
  grid_->levelIndexSets_[0]->setSizesAndTypes(vertexPositions_.size(), elements_.size());

  grid_->leafIndexSet_.setSizesAndTypes(vertexPositions_.size(), elements_.size());

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
