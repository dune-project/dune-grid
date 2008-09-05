// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include "../onedgrid.hh"

// ///////////////////////////////////////////////////////////////
//
//    OneDGridLevelIteratorFactory, a class used to simulate
//    specialization of member templates
//
// ///////////////////////////////////////////////////////////////

namespace Dune {

  template <int codim>
  struct OneDGridLevelIteratorFactory {};

  template <>
  struct OneDGridLevelIteratorFactory<1>
  {
    template <Dune::PartitionIteratorType PiType>
    static Dune::OneDGridLevelIterator<1,PiType, const Dune::OneDGrid>
    lbegin(const Dune::OneDGrid* g, int level) {

      return Dune::OneDGridLevelIterator<1,PiType, const Dune::OneDGrid>(g->vertices[level].begin);
    }

  };

  template <>
  struct OneDGridLevelIteratorFactory<0>
  {
    template <Dune::PartitionIteratorType PiType>
    static Dune::OneDGridLevelIterator<0,PiType, const Dune::OneDGrid>
    lbegin(const Dune::OneDGrid* g, int level) {

      return Dune::OneDGridLevelIterator<0,PiType, const Dune::OneDGrid>(g->elements[level].begin);
    }

  };

}


Dune::OneDGrid::OneDGrid(int numElements, const ctype& leftBoundary, const ctype& rightBoundary)
  : refinementType_(LOCAL),
    leafIndexSet_(*this),
    idSet_(*this),
    freeVertexIdCounter_(0),
    freeElementIdCounter_(0)
{
  if (numElements<1)
    DUNE_THROW(GridError, "Nonpositive number of elements requested!");

  if (leftBoundary >= rightBoundary)
    DUNE_THROW(GridError, "The left boundary coordinate has to be strictly less than the right boundary one!");
  // Init grid hierarchy
  vertices.resize(1);
  elements.resize(1);

  // Init vertex set
  for (int i=0; i<numElements+1; i++) {
    ctype newCoord = leftBoundary + i*(rightBoundary-leftBoundary) / numElements;

    OneDEntityImp<0>* newVertex = new OneDEntityImp<0>(0, newCoord);
    newVertex->id_ = getNextFreeId(1);
    vertices[0].insert_after(vertices[0].rbegin, newVertex);

  }

  // Init element set
  OneDEntityImp<0>* it = vertices[0].begin;
  for (int i=0; i<numElements; i++) {

    OneDEntityImp<1>* newElement = new OneDEntityImp<1>(0, getNextFreeId(0));
    newElement->vertex_[0] = it;
    it = it->succ_;
    newElement->vertex_[1] = it;

    elements[0].insert_after(elements[0].rbegin, newElement);

  }

  setIndices();
}

Dune::OneDGrid::OneDGrid(const std::vector<ctype>& coords)
  : refinementType_(LOCAL),
    leafIndexSet_(*this),
    idSet_(*this),
    freeVertexIdCounter_(0),
    freeElementIdCounter_(0)
{
  if (coords.size()<2)
    DUNE_THROW(GridError, "You have to provide at least two coordinates!");

  // Init grid hierarchy
  vertices.resize(1);
  elements.resize(1);

  // Init vertex set
  for (size_t i=0; i<coords.size(); i++) {
    OneDEntityImp<0>* newVertex = new OneDEntityImp<0>(0, coords[i], getNextFreeId(1));
    vertices[0].insert_after(vertices[0].rbegin, newVertex);
  }

  // Init element set
  OneDEntityImp<0>* it = vertices[0].begin;
  for (size_t i=0; i<coords.size()-1; i++) {

    OneDEntityImp<1>* newElement = new OneDEntityImp<1>(0, getNextFreeId(0));
    newElement->vertex_[0] = it;
    it = it->succ_;
    newElement->vertex_[1] = it;

    if (newElement->vertex_[0]->pos_ >= newElement->vertex_[1]->pos_)
      DUNE_THROW(GridError, "The coordinates have to be in ascending order!");

    elements[0].insert_after(elements[0].rbegin, newElement);

  }

  setIndices();
}


Dune::OneDGrid::~OneDGrid()
{
  // Delete all vertices
  for (unsigned int i=0; i<vertices.size(); i++) {

    OneDEntityImp<0>* v = vertices[i].begin;

    while (v) {

      OneDEntityImp<0>* vSucc = v->succ_;
      vertices[i].remove(v);
      delete(v);
      v = vSucc;

    }

  }

  // Delete all elements
  for (unsigned int i=0; i<elements.size(); i++) {

    OneDEntityImp<1>* e = elements[i].begin;

    while (e) {

      OneDEntityImp<1>* eSucc = e->succ_;
      elements[i].remove(e);
      delete(e);
      e = eSucc;

    }

  }

  // Delete levelIndexSets
  for (unsigned int i=0; i<levelIndexSets_.size(); i++)
    if (levelIndexSets_[i])
      delete levelIndexSets_[i];
}

template <int codim>
typename Dune::OneDGrid::Traits::template Codim<codim>::LevelIterator
Dune::OneDGrid::lbegin(int level) const
{
  if (level<0 || level>maxLevel())
    DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");

  return OneDGridLevelIteratorFactory<codim>::template lbegin<All_Partition>(this, level);
}

template <int codim>
typename Dune::OneDGrid::Traits::template Codim<codim>::LevelIterator
Dune::OneDGrid::lend(int level) const
{
  if (level<0 || level>maxLevel())
    DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");

  return OneDGridLevelIterator<codim,All_Partition, const Dune::OneDGrid>(0);
}

template <int codim, Dune::PartitionIteratorType PiType>
typename Dune::OneDGrid::Traits::template Codim<codim>::template Partition<PiType>::LevelIterator
Dune::OneDGrid::lbegin(int level) const
{
  if (level<0 || level>maxLevel())
    DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");

  return OneDGridLevelIteratorFactory<codim>::template lbegin<PiType>(this, level);
}

template <int codim, Dune::PartitionIteratorType PiType>
typename Dune::OneDGrid::Traits::template Codim<codim>::template Partition<PiType>::LevelIterator
Dune::OneDGrid::lend(int level) const
{
  if (level<0 || level>maxLevel())
    DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");

  return OneDGridLevelIterator<codim,PiType, const Dune::OneDGrid>(0);
}


template <int codim>
typename Dune::OneDGrid::Traits::template Codim<codim>::LeafIterator
Dune::OneDGrid::leafbegin() const
{
  return OneDGridLeafIterator<codim,All_Partition,const OneDGrid>(*this);
}

template <int codim>
typename Dune::OneDGrid::Traits::template Codim<codim>::LeafIterator
Dune::OneDGrid::leafend() const
{
  return OneDGridLeafIterator<codim,All_Partition, const Dune::OneDGrid>();
}

template <int codim, Dune::PartitionIteratorType PiType>
typename Dune::OneDGrid::Traits::template Codim<codim>::template Partition<PiType>::LeafIterator
Dune::OneDGrid::leafbegin() const
{
  return OneDGridLeafIterator<codim,PiType, const Dune::OneDGrid>(*this);
}

template <int codim, Dune::PartitionIteratorType PiType>
typename Dune::OneDGrid::Traits::template Codim<codim>::template Partition<PiType>::LeafIterator
Dune::OneDGrid::leafend() const
{
  return OneDGridLeafIterator<codim,PiType, const Dune::OneDGrid>();
}

Dune::OneDEntityImp<0>*
Dune::OneDGrid::getLeftUpperVertex(const OneDEntityImp<1>* eIt)
{
  OneDEntityImp<1>* l = eIt->pred_;

  if (!l)
    return 0;

  // return NULL if there is no geometrical left neighbor
  if (l->vertex_[1] != eIt->vertex_[0])
    return 0;

  // return NULL if that neighbor doesn't have sons
  if (l->isLeaf())
    return 0;

  // return the right vertex of the right son
  return l->sons_[1]->vertex_[1];

}

Dune::OneDEntityImp<0>*
Dune::OneDGrid::getRightUpperVertex(const OneDEntityImp<1>* eIt)
{
  OneDEntityImp<1>* r = eIt->succ_;

  if (!r)
    return 0;

  // return NULL if there is no geometrical right neighbor
  if (r->vertex_[0]!=eIt->vertex_[1])
    return 0;

  // return NULL if that neighbor doesn't have sons
  if (r->isLeaf())
    return 0;

  // return the left vertex of the left son
  return r->sons_[0]->vertex_[0];

}

Dune::OneDEntityImp<1>*
Dune::OneDGrid::getLeftNeighborWithSon(OneDEntityImp<1>* eIt)
{
  OneDEntityImp<1>* l = eIt;

  do {
    l = l->pred_;
  } while (l && l->isLeaf());

  return l;
}


bool Dune::OneDGrid::adapt()
{
  OneDEntityImp<1>* eIt;

  // for the return value:  true if the grid was changed
  bool changedGrid = false;

  // remove all elements that have been marked for coarsening
  for (int i=1; i<=maxLevel(); i++) {

    for (eIt = elements[i].begin; eIt!=NULL; ) {

      OneDEntityImp<1>* leftElementToBeDeleted  = eIt;
      OneDEntityImp<1>* rightElementToBeDeleted = eIt->succ_;

      assert(eIt->succ_);
      OneDEntityImp<1>* nextElement = eIt->succ_->succ_;

      if (leftElementToBeDeleted->markState_ == OneDEntityImp<1>::COARSEN && leftElementToBeDeleted->isLeaf()
          && rightElementToBeDeleted->markState_ == OneDEntityImp<1>::COARSEN && rightElementToBeDeleted->isLeaf()) {

        assert(rightElementToBeDeleted->isLeaf());

        // Is the left vertex obsolete?
        if (leftElementToBeDeleted->pred_==NULL
            || leftElementToBeDeleted->pred_->vertex_[1] != leftElementToBeDeleted->vertex_[0]) {

          // If the left vertex has a father remove the reference to this vertex at this father
          assert(leftElementToBeDeleted->father_->vertex_[0]->son_ == leftElementToBeDeleted->vertex_[0]);
          leftElementToBeDeleted->father_->vertex_[0]->son_ = NULL;

          vertices[i].remove(leftElementToBeDeleted->vertex_[0]);
          delete(leftElementToBeDeleted->vertex_[0]);
        }

        // Is the right vertex obsolete?
        if (rightElementToBeDeleted->succ_==NULL
            || rightElementToBeDeleted->succ_->vertex_[0] != rightElementToBeDeleted->vertex_[1]) {

          // If the left vertex has a father remove the reference to this vertex at this father
          assert(rightElementToBeDeleted->father_->vertex_[1]->son_ == rightElementToBeDeleted->vertex_[1]);
          rightElementToBeDeleted->father_->vertex_[1]->son_ = NULL;

          vertices[i].remove(rightElementToBeDeleted->vertex_[1]);
          delete(rightElementToBeDeleted->vertex_[1]);
        }

        // Delete vertex between left and right element to be deleted
        assert(leftElementToBeDeleted->vertex_[1] == rightElementToBeDeleted->vertex_[0]);
        vertices[i].remove(leftElementToBeDeleted->vertex_[1]);
        delete(leftElementToBeDeleted->vertex_[1]);

        // Remove references from the father element
        assert(rightElementToBeDeleted->father_->sons_[1] == rightElementToBeDeleted);
        leftElementToBeDeleted->father_->sons_[0]  = NULL;
        rightElementToBeDeleted->father_->sons_[1] = NULL;

        // Paranoia: make sure the father is not marked for refinement
        rightElementToBeDeleted->father_->markState_ = OneDEntityImp<1>::NONE;

        // Actually delete elements
        elements[i].remove(leftElementToBeDeleted);
        elements[i].remove(rightElementToBeDeleted);
        delete(leftElementToBeDeleted);
        delete(rightElementToBeDeleted);

        // The grid has been changed
        changedGrid = true;
      }

      // increment pointer
      eIt = nextElement;

    }

  }

  // /////////////////////////////////////////////////////////////////////////
  //  Check if one of the elements on the toplevel is marked for refinement
  //  In that case add another level
  // /////////////////////////////////////////////////////////////////////////
  bool toplevelRefinement = false;
  for (eIt = elements[maxLevel()].begin; eIt!=NULL; eIt=eIt->succ_)
    if (eIt->markState_ == OneDEntityImp<1>::REFINED) {
      toplevelRefinement = true;
      break;
    }

  if (toplevelRefinement) {
    List<OneDEntityImp<0> > newVertices;
    List<OneDEntityImp<1> > newElements;
    vertices.push_back(newVertices);
    elements.push_back(newElements);
  }

  // //////////////////////////////
  // refine all marked elements
  // //////////////////////////////
  int oldMaxlevel = (toplevelRefinement) ? maxLevel()-1 : maxLevel();
  for (int i=0; i<=oldMaxlevel; i++) {

    for (eIt = elements[i].begin; eIt!=NULL; eIt = eIt->succ_) {

      if (eIt->markState_ == OneDEntityImp<1>::REFINED
          && eIt->isLeaf()) {

        // Does the left vertex exist on the next-higher level?
        // If no create it
        OneDEntityImp<0>* leftUpperVertex = getLeftUpperVertex(eIt);

        if (leftUpperVertex==NULL)
          leftUpperVertex = new OneDEntityImp<0>(i+1,
                                                 eIt->vertex_[0]->pos_,
                                                 eIt->vertex_[0]->id_);

        eIt->vertex_[0]->son_ = leftUpperVertex;

        // Does the right vertex exist on the next-higher level?
        // If no create it
        OneDEntityImp<0>* rightUpperVertex = getRightUpperVertex(eIt);

        if (rightUpperVertex==NULL)
          rightUpperVertex = new OneDEntityImp<0>(i+1,
                                                  eIt->vertex_[1]->pos_,
                                                  eIt->vertex_[1]->id_);

        eIt->vertex_[1]->son_ = rightUpperVertex;

        // Create center vertex
        ctype p = 0.5*(eIt->vertex_[0]->pos_[0] + eIt->vertex_[1]->pos_[0]);

        OneDEntityImp<0>* centerVertex = new OneDEntityImp<0>(i+1, p, getNextFreeId(1));

        // //////////////////////////////////////
        // Insert new vertices into vertex list
        // //////////////////////////////////////

        OneDEntityImp<1>* leftNeighbor = getLeftNeighborWithSon(eIt);

        if (leftNeighbor!=NULL) {

          // leftNeighbor exists
          if ( leftNeighbor->sons_[1]->vertex_[1] != leftUpperVertex)
            vertices[i+1].insert_after(leftNeighbor->sons_[1]->vertex_[1], leftUpperVertex);

        } else {
          // leftNeighbor does not exist
          vertices[i+1].insert_before(vertices[i+1].begin, leftUpperVertex);

        }

        vertices[i+1].insert_after(leftUpperVertex, centerVertex);

        // Check if rightUpperVertex is already in the list
        OneDEntityImp<0>* succOfCenter = centerVertex->succ_;

        if (succOfCenter==NULL || succOfCenter != rightUpperVertex)
          vertices[i+1].insert_after(centerVertex, rightUpperVertex);

        // ///////////////////////
        // Create new elements
        // ///////////////////////
        OneDEntityImp<1>* newElement0 = new OneDEntityImp<1>(i+1, getNextFreeId(0));
        newElement0->vertex_[0] = leftUpperVertex;
        newElement0->vertex_[1] = centerVertex;
        newElement0->father_ = eIt;
        newElement0->adaptationState_ = OneDEntityImp<1>::REFINED;

        OneDEntityImp<1>* newElement1 = new OneDEntityImp<1>(i+1, getNextFreeId(0));
        newElement1->vertex_[0] = centerVertex;
        newElement1->vertex_[1] = rightUpperVertex;
        newElement1->father_ = eIt;
        newElement1->adaptationState_ = OneDEntityImp<1>::REFINED;

        // Insert new elements into element list
        if (leftNeighbor!=NULL)
          // leftNeighbor exists
          elements[i+1].insert_after(leftNeighbor->sons_[1], newElement0);
        else
          // leftNeighbor does not exist
          elements[i+1].insert_before(elements[i+1].begin, newElement0);

        elements[i+1].insert_after(newElement0, newElement1);

        // Mark the two new elements as the sons of the refined element
        eIt->sons_[0] = newElement0;
        eIt->sons_[1] = newElement1;

        // The grid has been modified
        changedGrid = true;

      }

    }

  }

  // delete uppermost level if it doesn't contain elements anymore
  if (elements[maxLevel()].size()==0) {
    assert(vertices[maxLevel()].size()==0);
    elements.pop_back();
    vertices.pop_back();
  }


  // If the refinement mode is 'COPY', fill the empty slots in the grid
  // by copying elements
  if (refinementType_ == COPY) {

    for (int i=0; i<maxLevel(); i++) {

      OneDEntityImp<1>* eIt;
      for (eIt = elements[i].begin; eIt!=NULL; eIt = eIt->succ_) {

        if (eIt->isLeaf()) {

          // Does the left vertex exist on the next-higher level?
          // If no create it
          OneDEntityImp<0>* leftUpperVertex = getLeftUpperVertex(eIt);

          if (leftUpperVertex==NULL)
            leftUpperVertex = new OneDEntityImp<0>(i+1, eIt->vertex_[0]->pos_, eIt->vertex_[0]->id_);

          eIt->vertex_[0]->son_ = leftUpperVertex;

          // Does the right vertex exist on the next-higher level?
          // If no create it
          OneDEntityImp<0>* rightUpperVertex = getRightUpperVertex(eIt);

          if (rightUpperVertex==NULL)
            rightUpperVertex = new OneDEntityImp<0>(i+1, eIt->vertex_[1]->pos_, eIt->vertex_[1]->id_);

          eIt->vertex_[1]->son_ = rightUpperVertex;

          // //////////////////////////////////////
          // Insert new vertices into vertex list
          // //////////////////////////////////////

          OneDEntityImp<1>* leftNeighbor = getLeftNeighborWithSon(eIt);

          if (leftNeighbor!=NULL) {

            // leftNeighbor exists
            if ( leftNeighbor->sons_[1]->vertex_[1] != leftUpperVertex)
              vertices[i+1].insert_after(leftNeighbor->sons_[1]->vertex_[1], leftUpperVertex);

          } else {
            // leftNeighbor does not exist
            vertices[i+1].insert_before(vertices[i+1].begin, leftUpperVertex);

          }

          // Check if rightUpperVertex is already in the list
          OneDEntityImp<0>* succOfLeft = leftUpperVertex->succ_;

          if (succOfLeft==NULL || succOfLeft != rightUpperVertex)
            vertices[i+1].insert_after(leftUpperVertex, rightUpperVertex);

          // /////////////////////////
          //   Create new element
          // /////////////////////////
          OneDEntityImp<1>* newElement = new OneDEntityImp<1>(i+1, eIt->id_);
          newElement->vertex_[0] = leftUpperVertex;
          newElement->vertex_[1] = rightUpperVertex;
          newElement->father_ = eIt;
          newElement->adaptationState_ = OneDEntityImp<1>::REFINED;

          // Insert new elements into element list
          if (leftNeighbor!=NULL)
            // leftNeighbor exists
            elements[i+1].insert_after(leftNeighbor->sons_[1], newElement);
          else
            // leftNeighbor does not exist
            elements[i+1].insert_before(elements[i+1].begin, newElement);

          // Mark the new element as the sons of the refined element
          eIt->sons_[0] = eIt->sons_[1] = newElement;

        }

      }

    }

  }

  // ////////////////////////////////////
  //   renumber vertices and elements
  // ////////////////////////////////////
  setIndices();

  return changedGrid;
}

bool Dune::OneDGrid::preAdapt()
{
  Codim<0>::LeafIterator eIt    = leafbegin<0>();
  Codim<0>::LeafIterator eEndIt = leafend<0>();

  for (; eIt!=eEndIt; ++eIt)
    if (getRealImplementation(*eIt).target_->markState_ != OneDEntityImp<1>::NONE)
      return true;

  return false;
}

void Dune::OneDGrid::postAdapt()
{
  for (int i=0; i<=maxLevel(); i++) {
    OneDEntityImp<1>* eIt;
    for (eIt = elements[i].begin; eIt!=NULL; eIt = eIt->succ_)
      eIt->markState_ = OneDEntityImp<1>::NONE;

  }

}

void Dune::OneDGrid::setIndices()
{
  // Add space for new LevelIndexSets if the grid hierarchy got higher
  // They are not created until they are actually requested
  for (int i=levelIndexSets_.size(); i<maxLevel()+1; i++)
    levelIndexSets_.push_back(0);

  // Delete old LevelIndexSets if the grid hierarchy got lower
  int excess = levelIndexSets_.size() - (maxLevel() + 1);
  for (int i=0; i<excess; i++) {
    if (levelIndexSets_.back())
      delete(levelIndexSets_.back());
    levelIndexSets_.pop_back();
  }

  for (int i=0; i<=maxLevel(); i++)
    if (levelIndexSets_[i])
      levelIndexSets_[i]->update();

  leafIndexSet_.update();

  idSet_.update();

}

void Dune::OneDGrid::globalRefine(int refCount)
{
  for (int i=0; i<refCount; i++) {

    // mark all entities for grid refinement
    Codim<0>::LeafIterator iIt    = leafbegin<0>();
    Codim<0>::LeafIterator iEndIt = leafend<0>();

    for (; iIt!=iEndIt; ++iIt)
      mark(1, *iIt);

    this->preAdapt();
    adapt();
    this->postAdapt();
  }
}

bool Dune::OneDGrid::mark(int refCount,
                          const Codim<0>::EntityPointer & ep )
{
  return this->mark( refCount, *ep);
}

bool Dune::OneDGrid::mark(int refCount,
                          const Codim<0>::Entity & e )
{
  if (refCount < 0) {

    if (getRealImplementation(e).target_->level_ == 0)
      return false;
    else {
      getRealImplementation(e).target_->markState_ = OneDEntityImp<1>::COARSEN;
      return true;
    }

  } else if (refCount > 0)
    getRealImplementation(e).target_->markState_ = OneDEntityImp<1>::REFINED;
  else
    getRealImplementation(e).target_->markState_ = OneDEntityImp<1>::NONE;

  return true;
}

int Dune::OneDGrid::getMark(const Codim<0>::EntityPointer & ep ) const
{
  return this->getMark( *ep );
}

int Dune::OneDGrid::getMark(const Codim<0>::Entity & e ) const
{
  if(getRealImplementation(e).target_->markState_ == OneDEntityImp<1>::COARSEN)
    return -1;
  else if(getRealImplementation(e).target_->markState_ == OneDEntityImp<1>::REFINED)
    return 1;
  return 0;
}



// /////////////////////////////////////////////////////////////////////////
//   Explicitly instantiate the OnedGrid member templates.
//   gcc-4.0 wants these instantiations after the method implementations
// /////////////////////////////////////////////////////////////////////////

template Dune::OneDGrid::Codim<0>::LevelIterator Dune::OneDGrid::lbegin<0>(int level) const;
template Dune::OneDGrid::Codim<1>::LevelIterator Dune::OneDGrid::lbegin<1>(int level) const;

template Dune::OneDGrid::Codim<0>::LevelIterator Dune::OneDGrid::lend<0>(int level) const;
template Dune::OneDGrid::Codim<1>::LevelIterator Dune::OneDGrid::lend<1>(int level) const;

template Dune::OneDGrid::Codim<0>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::OneDGrid::lbegin<0,Dune::Interior_Partition>(int level) const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::OneDGrid::lbegin<0,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::OneDGrid::lbegin<0,Dune::Overlap_Partition>(int level) const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDGrid::lbegin<0,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDGrid::lbegin<0,Dune::All_Partition>(int level) const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDGrid::lbegin<0,Dune::Ghost_Partition>(int level) const;

template Dune::OneDGrid::Codim<0>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::OneDGrid::lend<0,Dune::Interior_Partition>(int level) const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::OneDGrid::lend<0,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::OneDGrid::lend<0,Dune::Overlap_Partition>(int level) const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDGrid::lend<0,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDGrid::lend<0,Dune::All_Partition>(int level) const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDGrid::lend<0,Dune::Ghost_Partition>(int level) const;

template Dune::OneDGrid::Codim<1>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::OneDGrid::lbegin<1,Dune::Interior_Partition>(int level) const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::OneDGrid::lbegin<1,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::OneDGrid::lbegin<1,Dune::Overlap_Partition>(int level) const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDGrid::lbegin<1,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDGrid::lbegin<1,Dune::All_Partition>(int level) const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDGrid::lbegin<1,Dune::Ghost_Partition>(int level) const;

template Dune::OneDGrid::Codim<1>::Partition<Dune::Interior_Partition>::LevelIterator
Dune::OneDGrid::lend<1,Dune::Interior_Partition>(int level) const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LevelIterator
Dune::OneDGrid::lend<1,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::Overlap_Partition>::LevelIterator
Dune::OneDGrid::lend<1,Dune::Overlap_Partition>(int level) const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDGrid::lend<1,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDGrid::lend<1,Dune::All_Partition>(int level) const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDGrid::lend<1,Dune::Ghost_Partition>(int level) const;

template Dune::OneDGrid::Codim<0>::LeafIterator Dune::OneDGrid::leafbegin<0>() const;
template Dune::OneDGrid::Codim<1>::LeafIterator Dune::OneDGrid::leafbegin<1>() const;

template Dune::OneDGrid::Codim<0>::LeafIterator Dune::OneDGrid::leafend<0>() const;
template Dune::OneDGrid::Codim<1>::LeafIterator Dune::OneDGrid::leafend<1>() const;

template Dune::OneDGrid::Codim<0>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<0,Dune::Interior_Partition>() const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<0,Dune::InteriorBorder_Partition>() const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<0,Dune::Overlap_Partition>() const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<0,Dune::OverlapFront_Partition>() const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<0,Dune::All_Partition>() const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<0,Dune::Ghost_Partition>() const;

template Dune::OneDGrid::Codim<0>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::OneDGrid::leafend<0,Dune::Interior_Partition>() const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::OneDGrid::leafend<0,Dune::InteriorBorder_Partition>() const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::OneDGrid::leafend<0,Dune::Overlap_Partition>() const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDGrid::leafend<0,Dune::OverlapFront_Partition>() const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDGrid::leafend<0,Dune::All_Partition>() const;
template Dune::OneDGrid::Codim<0>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDGrid::leafend<0,Dune::Ghost_Partition>() const;

template Dune::OneDGrid::Codim<1>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<1,Dune::Interior_Partition>() const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<1,Dune::InteriorBorder_Partition>() const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<1,Dune::Overlap_Partition>() const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<1,Dune::OverlapFront_Partition>() const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<1,Dune::All_Partition>() const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDGrid::leafbegin<1,Dune::Ghost_Partition>() const;

template Dune::OneDGrid::Codim<1>::Partition<Dune::Interior_Partition>::LeafIterator
Dune::OneDGrid::leafend<1,Dune::Interior_Partition>() const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LeafIterator
Dune::OneDGrid::leafend<1,Dune::InteriorBorder_Partition>() const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::Overlap_Partition>::LeafIterator
Dune::OneDGrid::leafend<1,Dune::Overlap_Partition>() const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDGrid::leafend<1,Dune::OverlapFront_Partition>() const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDGrid::leafend<1,Dune::All_Partition>() const;
template Dune::OneDGrid::Codim<1>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDGrid::leafend<1,Dune::Ghost_Partition>() const;
