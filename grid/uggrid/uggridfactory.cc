// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/grid/uggrid/uggridfactory.hh>

template <int dimworld>
Dune::GridFactory<Dune::UGGrid<dimworld> >::
GridFactory()
{
  grid_ = new Dune::UGGrid<dimworld>;

  factoryOwnsGrid_ = true;

  grid_->createBegin();
}

template <int dimworld>
Dune::GridFactory<Dune::UGGrid<dimworld> >::
GridFactory(UGGrid<dimworld>* grid)
{
  grid_ = grid;

  factoryOwnsGrid_ = false;

  grid_->createBegin();
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
  grid_->insertVertex(pos);
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



// Explicit template instatiation
template class Dune::GridFactory<Dune::UGGrid<2> >;
template class Dune::GridFactory<Dune::UGGrid<3> >;
