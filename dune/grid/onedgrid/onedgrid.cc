// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include "../onedgrid.hh"


Dune::OneDGrid::OneDGrid()
  : refinementType_(LOCAL),
    leafIndexSet_(*this),
    idSet_(*this),
    freeIdCounter_(0),
    reversedBoundarySegmentNumbering_(false)
{}

Dune::OneDGrid::OneDGrid(int numElements, const ctype& leftBoundary, const ctype& rightBoundary)
  : refinementType_(LOCAL),
    leafIndexSet_(*this),
    idSet_(*this),
    freeIdCounter_(0),
    reversedBoundarySegmentNumbering_(false)
{
  if (numElements<1)
    DUNE_THROW(GridError, "Nonpositive number of elements requested!");

  if (leftBoundary >= rightBoundary)
    DUNE_THROW(GridError, "The left boundary coordinate has to be strictly less than the right boundary one!");

  // Init grid hierarchy
  entityImps_.resize(1);

  // Init vertex set
  for (int i=0; i<numElements+1; i++) {
    ctype newCoord = leftBoundary + i*(rightBoundary-leftBoundary) / numElements;

    OneDEntityImp<0> newVertex(0, newCoord);
    newVertex.id_ = getNextFreeId();
    vertices(0).push_back(newVertex);

  }

  // Init element set
  OneDGridList<OneDEntityImp<0> >::iterator it = vertices(0).begin();
  for (int i=0; i<numElements; i++) {

    OneDEntityImp<1> newElement(0, getNextFreeId(), false);
    newElement.vertex_[0] = it;
    it = it->succ_;
    newElement.vertex_[1] = it;

    elements(0).push_back(newElement);

  }

  setIndices();
}

Dune::OneDGrid::OneDGrid(const std::vector<ctype>& coords)
  : refinementType_(LOCAL),
    leafIndexSet_(*this),
    idSet_(*this),
    freeIdCounter_(0),
    reversedBoundarySegmentNumbering_(false)
{
  if (coords.size()<2)
    DUNE_THROW(GridError, "You have to provide at least two coordinates!");

  // Init grid hierarchy
  entityImps_.resize(1);

  // Init vertex set
  for (size_t i=0; i<coords.size(); i++) {
    OneDEntityImp<0> newVertex(0, coords[i], getNextFreeId());
    vertices(0).push_back(newVertex);
  }

  // Init element set
  OneDGridList<OneDEntityImp<0> >::iterator it = vertices(0).begin();
  for (size_t i=0; i<coords.size()-1; i++) {

    OneDEntityImp<1> newElement(0, getNextFreeId(), false);
    newElement.vertex_[0] = it;
    it = it->succ_;
    newElement.vertex_[1] = it;

    if (newElement.vertex_[0]->pos_ >= newElement.vertex_[1]->pos_)
      DUNE_THROW(GridError, "The coordinates have to be in ascending order!");

    elements(0).push_back(newElement);

  }

  setIndices();
}


Dune::OneDGrid::~OneDGrid()
{
  // Delete all vertices
  for (unsigned int i=0; i<entityImps_.size(); i++) {

    OneDGridList<OneDEntityImp<0> >::iterator v = vertices(i).begin();

    while (v) {

      OneDGridList<OneDEntityImp<0> >::iterator vSucc = v->succ_;
      vertices(i).erase(v);
      v = vSucc;

    }

  }

  // Delete all elements
  for (unsigned int i=0; i<entityImps_.size(); i++) {

    OneDGridList<OneDEntityImp<1> >::iterator e = elements(i).begin();

    while (e) {

      OneDGridList<OneDEntityImp<1> >::iterator eSucc = e->succ_;
      elements(i).erase(e);
      e = eSucc;

    }

  }

  // Delete levelIndexSets
  for (unsigned int i=0; i<levelIndexSets_.size(); i++)
    if (levelIndexSets_[i])
      delete levelIndexSets_[i];
}

Dune::OneDGridList<Dune::OneDEntityImp<0> >::iterator
Dune::OneDGrid::getLeftUpperVertex(const OneDEntityImp<1>* eIt)
{
  OneDGridList<OneDEntityImp<1> >::iterator l = eIt->pred_;

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

Dune::OneDGridList<Dune::OneDEntityImp<0> >::iterator
Dune::OneDGrid::getRightUpperVertex(const OneDEntityImp<1>* eIt)
{
  OneDGridList<OneDEntityImp<1> >::iterator r = eIt->succ_;

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

Dune::OneDGridList<Dune::OneDEntityImp<1> >::iterator
Dune::OneDGrid::getLeftNeighborWithSon(OneDGridList<OneDEntityImp<1> >::iterator eIt)
{
  OneDGridList<OneDEntityImp<1> >::iterator l = eIt;

  do {
    l = l->pred_;
  } while (l && l->isLeaf());

  return l;
}


bool Dune::OneDGrid::adapt()
{
  // for the return value:  true if the grid was changed
  bool refinedGrid = false;

  // remove elements that have been marked for coarsening.
  // If one son of an element is marked for coarsening, and the other one is not,
  // then the element is not removed.
  for (int i=1; i<=maxLevel(); i++) {

    for (auto eIt = elements(i).begin(); eIt!=elements(i).end(); ) {

      OneDGridList<OneDEntityImp<1> >::iterator leftElementToBeDeleted  = eIt;
      OneDGridList<OneDEntityImp<1> >::iterator rightElementToBeDeleted = eIt->succ_;

      assert(eIt->succ_);
      OneDGridList<OneDEntityImp<1> >::iterator nextElement = eIt->succ_->succ_;

      if (leftElementToBeDeleted->markState_ ==
          OneDEntityImp<1>::COARSEN && leftElementToBeDeleted->isLeaf()
          && rightElementToBeDeleted->markState_ ==
          OneDEntityImp<1>::COARSEN && rightElementToBeDeleted->isLeaf()) {

        assert(rightElementToBeDeleted->isLeaf());

        // Is the left vertex obsolete?
        if (leftElementToBeDeleted->pred_==NULL
            || leftElementToBeDeleted->pred_->vertex_[1] != leftElementToBeDeleted->vertex_[0]) {

          // If the left vertex has a father remove the reference to this vertex at this father
          assert(leftElementToBeDeleted->father_->vertex_[0]->son_ == leftElementToBeDeleted->vertex_[0]);
          leftElementToBeDeleted->father_->vertex_[0]->son_ = NULL;

          vertices(i).erase(leftElementToBeDeleted->vertex_[0]);
        }

        // Is the right vertex obsolete?
        if (rightElementToBeDeleted->succ_==NULL
            || rightElementToBeDeleted->succ_->vertex_[0] != rightElementToBeDeleted->vertex_[1]) {

          // If the left vertex has a father remove the reference to this vertex at this father
          assert(rightElementToBeDeleted->father_->vertex_[1]->son_ == rightElementToBeDeleted->vertex_[1]);
          rightElementToBeDeleted->father_->vertex_[1]->son_ = NULL;

          vertices(i).erase(rightElementToBeDeleted->vertex_[1]);
        }

        // Delete vertex between left and right element to be deleted
        assert(leftElementToBeDeleted->vertex_[1] == rightElementToBeDeleted->vertex_[0]);
        vertices(i).erase(leftElementToBeDeleted->vertex_[1]);

        // Remove references from the father element
        assert(rightElementToBeDeleted->father_->sons_[1] == rightElementToBeDeleted);
        leftElementToBeDeleted->father_->sons_[0]  = NULL;
        rightElementToBeDeleted->father_->sons_[1] = NULL;

        // Paranoia: make sure the father is not marked for refinement
        rightElementToBeDeleted->father_->markState_ = OneDEntityImp<1>::DO_NOTHING;

        // Actually delete elements
        elements(i).erase(leftElementToBeDeleted);
        elements(i).erase(rightElementToBeDeleted);
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
  for (auto eIt =  elements(maxLevel()).begin();
       eIt != elements(maxLevel()).end();
       eIt=eIt->succ_)
    if (eIt->markState_ == OneDEntityImp<1>::REFINE) {
      toplevelRefinement = true;
      break;
    }

  if (toplevelRefinement) {
    OneDGridList<OneDEntityImp<0> > newVertices;
    OneDGridList<OneDEntityImp<1> > newElements;
    entityImps_.push_back(std::tuple<OneDGridList<OneDEntityImp<0> >, OneDGridList<OneDEntityImp<1> > >(newVertices, newElements));
  }

  // //////////////////////////////
  // refine all marked elements
  // //////////////////////////////
  int oldMaxlevel = (toplevelRefinement) ? maxLevel()-1 : maxLevel();
  for (int i=0; i<=oldMaxlevel; i++) {

    for (auto eIt = elements(i).begin(); eIt!=elements(i).end(); eIt = eIt->succ_) {

      if (eIt->markState_ == OneDEntityImp<1>::REFINE
          && eIt->isLeaf()) {

        OneDGridList<OneDEntityImp<1> >::iterator leftNeighbor = getLeftNeighborWithSon(eIt);

        // ////////////////////////////////////////////////////////////
        // Does the left vertex exist on the next-higher level?
        // ////////////////////////////////////////////////////////////

        OneDGridList<OneDEntityImp<0> >::iterator leftUpperVertex = getLeftUpperVertex(eIt);

        // If no create it
        if (leftUpperVertex==NULL) {

          OneDEntityImp<0> newLeftUpperVertex(i+1,
                                              eIt->vertex_[0]->pos_,
                                              eIt->vertex_[0]->id_);

          // Insert new vertex into vertex list
          leftUpperVertex = vertices(i+1).insert((leftNeighbor)
                                                 ? leftNeighbor->sons_[1]->vertex_[1]->succ_
                                                 : vertices(i+1).begin(),
                                                 newLeftUpperVertex);

        }

        eIt->vertex_[0]->son_ = leftUpperVertex;

        // //////////////////////////////////
        // Create new center vertex
        // //////////////////////////////////

        ctype p = 0.5*(eIt->vertex_[0]->pos_[0] + eIt->vertex_[1]->pos_[0]);

        OneDEntityImp<0> centerVertex(i+1, p, getNextFreeId());

        OneDGridList<OneDEntityImp<0> >::iterator centerVertexIterator = vertices(i+1).insert(leftUpperVertex->succ_, centerVertex);

        // ////////////////////////////////////////////////////////////
        // Does the right vertex exist on the next-higher level?
        // If no create it
        // ////////////////////////////////////////////////////////////
        OneDGridList<OneDEntityImp<0> >::iterator rightUpperVertex = getRightUpperVertex(eIt);

        if (rightUpperVertex==NULL) {
          OneDEntityImp<0> newRightUpperVertex(i+1,
                                               eIt->vertex_[1]->pos_,
                                               eIt->vertex_[1]->id_);

          rightUpperVertex = vertices(i+1).insert(centerVertexIterator->succ_, newRightUpperVertex);

        }

        eIt->vertex_[1]->son_ = rightUpperVertex;

        // ///////////////////////
        // Create new elements
        // ///////////////////////
        OneDEntityImp<1> newElement0(i+1, getNextFreeId(), reversedBoundarySegmentNumbering_);
        newElement0.vertex_[0] = leftUpperVertex;
        newElement0.vertex_[1] = centerVertexIterator;
        newElement0.father_ = eIt;
        newElement0.isNew_ = true;

        OneDEntityImp<1> newElement1(i+1, getNextFreeId(), reversedBoundarySegmentNumbering_);
        newElement1.vertex_[0] = centerVertexIterator;
        newElement1.vertex_[1] = rightUpperVertex;
        newElement1.father_ = eIt;
        newElement1.isNew_ = true;

        // Insert new elements into element list
        if (leftNeighbor!=NULL)
          // leftNeighbor exists
          eIt->sons_[0] = elements(i+1).insert(leftNeighbor->sons_[1]->succ_, newElement0);
        else
          // leftNeighbor does not exist
          eIt->sons_[0] = elements(i+1).insert(elements(i+1).begin(), newElement0);

        eIt->sons_[1] = elements(i+1).insert(eIt->sons_[0]->succ_, newElement1);

        // The grid has been modified
        refinedGrid = true;

      }

    }

  }

  // delete uppermost level if it doesn't contain elements anymore
  if (std::get<1>(entityImps_.back()).size()==0) {
    assert(std::get<0>(entityImps_.back()).size()==0);
    entityImps_.pop_back();
  }


  // If the refinement mode is 'COPY', fill the empty slots in the grid
  // by copying elements
  if (refinementType_ == COPY) {

    for (int i=0; i<maxLevel(); i++) {

      for (auto eIt = elements(i).begin(); eIt!=elements(i).end(); eIt = eIt->succ_) {

        if (eIt->isLeaf()) {

          OneDGridList<OneDEntityImp<1> >::iterator leftNeighbor = getLeftNeighborWithSon(eIt);

          // Does the left vertex exist on the next-higher level?
          // If no create it
          OneDGridList<OneDEntityImp<0> >::iterator leftUpperVertex = getLeftUpperVertex(eIt);

          if (leftUpperVertex==NULL) {

            OneDEntityImp<0> newLeftUpperVertex(i+1, eIt->vertex_[0]->pos_, eIt->vertex_[0]->id_);

            // Insert new vertex into vertex list
            leftUpperVertex = vertices(i+1).insert((leftNeighbor)
                                                   ? leftNeighbor->sons_[1]->vertex_[1]->succ_
                                                   : vertices(i+1).begin(),
                                                   newLeftUpperVertex);

          }

          eIt->vertex_[0]->son_ = leftUpperVertex;

          // Does the right vertex exist on the next-higher level?
          // If no create it
          OneDGridList<OneDEntityImp<0> >::iterator rightUpperVertex = getRightUpperVertex(eIt);

          if (rightUpperVertex==NULL) {

            OneDEntityImp<0> newRightUpperVertex(i+1, eIt->vertex_[1]->pos_, eIt->vertex_[1]->id_);

            // Insert new vertex into list
            rightUpperVertex = vertices(i+1).insert(leftUpperVertex->succ_, newRightUpperVertex);
          }

          eIt->vertex_[1]->son_ = rightUpperVertex;

          // /////////////////////////
          //   Create new element
          // /////////////////////////
          OneDEntityImp<1> newElement(i+1, eIt->id_, reversedBoundarySegmentNumbering_);
          newElement.vertex_[0] = leftUpperVertex;
          newElement.vertex_[1] = rightUpperVertex;
          newElement.father_ = eIt;
          newElement.isNew_ = true;

          // Insert new elements into element list
          OneDGridList<OneDEntityImp<1> >::iterator newElementIterator;
          if (leftNeighbor!=NULL)
            // leftNeighbor exists
            newElementIterator = elements(i+1).insert(leftNeighbor->sons_[1]->succ_, newElement);
          else
            // leftNeighbor does not exist
            newElementIterator = elements(i+1).insert(elements(i+1).begin(), newElement);

          // Mark the new element as the sons of the refined element
          eIt->sons_[0] = eIt->sons_[1] = newElementIterator;

        }

      }

    }

  }

  // ////////////////////////////////////
  //   renumber vertices and elements
  // ////////////////////////////////////
  setIndices();

  return refinedGrid;
}

bool Dune::OneDGrid::preAdapt()
{
  for (const auto& element : Dune::elements(leafGridView()))
    if (element.impl().target_->markState_ == OneDEntityImp<1>::COARSEN)
      return true;

  return false;
}

void Dune::OneDGrid::postAdapt()
{
  for (int i=0; i<=maxLevel(); i++) {
    OneDGridList<OneDEntityImp<1> >::iterator eIt;
    for (eIt = elements(i).begin(); eIt!=elements(i).end(); eIt = eIt->succ_)
    {
      eIt->isNew_ = false ;
      eIt->markState_ = OneDEntityImp<1>::DO_NOTHING;
    }

  }

}

void Dune::OneDGrid::setIndices()
{
  // Add space for new LevelIndexSets if the grid hierarchy got higher
  // They are not created until they are actually requested
  for (int i=levelIndexSets_.size(); i<maxLevel()+1; i++)
    levelIndexSets_.push_back( (OneDGridLevelIndexSet< const OneDGrid > *) 0 );

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

  // IdSets don't need updating
}

void Dune::OneDGrid::globalRefine(int refCount)
{
  for (int i=0; i<refCount; i++) {

    // mark all entities for grid refinement
    for (const auto& element : Dune::elements(leafGridView()))
      mark(1, element);

    this->preAdapt();
    adapt();
    this->postAdapt();
  }
}

bool Dune::OneDGrid::mark(int refCount,
                          const Codim<0>::Entity & e )
{
  // don't mark non-leaf entities
  if( ! e.isLeaf() ) return false ;

  if (refCount < 0) {

    if (e.impl().target_->level_ == 0)
      return false;
    else {
      e.impl().target_->markState_ = OneDEntityImp<1>::COARSEN;
      return true;
    }

  } else if (refCount > 0)
    e.impl().target_->markState_ = OneDEntityImp<1>::REFINE;
  else
    e.impl().target_->markState_ = OneDEntityImp<1>::DO_NOTHING;

  return true;
}

int Dune::OneDGrid::getMark(const Codim<0>::Entity & e ) const
{
  if(e.impl().target_->markState_ == OneDEntityImp<1>::COARSEN)
    return -1;
  else if(e.impl().target_->markState_ == OneDEntityImp<1>::REFINE)
    return 1;
  return 0;
}
