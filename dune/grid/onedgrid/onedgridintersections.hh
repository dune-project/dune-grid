// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GRID_INTERSECTIONS_HH
#define DUNE_ONE_D_GRID_INTERSECTIONS_HH

/** \file
 * \brief The OneDGridLevelIntersection and OneDGridLeafIntersection classes
 */

#include <dune/grid/onedgrid/onedgridentity.hh>

namespace Dune {

  /** \brief Intersection between two neighboring elements on a level grid */
  template<class GridImp>
  class OneDGridLevelIntersection
  {
    constexpr static int dim = GridImp::dimension;
    constexpr static int dimworld = GridImp::dimensionworld;

    // The corresponding iterator needs to access all members
    friend class OneDGridLevelIntersectionIterator<GridImp>;

    template<typename,typename>
    friend class Intersection;

    OneDGridLevelIntersection()
      : center_(nullptr)
      , neighbor_(-1) // marker for invalid intersection
    {}

    //! Constructor for a given grid entity and a given neighbor
    OneDGridLevelIntersection(OneDEntityImp<1>* center, int nb)
      : center_(center), neighbor_(nb)
    {}

    /** \brief Constructor creating the 'one-after-last'-iterator */
    OneDGridLevelIntersection(OneDEntityImp<1>* center)
      : center_(center), neighbor_(2)
    {}

    typedef typename GridImp::Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim< 1 >::LocalGeometryImpl LocalGeometryImpl;

  public:

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;

    //! equality
    bool equals(const OneDGridLevelIntersection<GridImp>& other) const {
      return (center_ == other.center_) && (neighbor_ == other.neighbor_);
    }

    OneDEntityImp<1>* target() const {
      const bool isValid = center_ && neighbor_>=0 && neighbor_<2;

      if (!isValid)
        return center_;
      else if (neighbor_==0)
        return center_->pred_;
      else
        return center_->succ_;

    }

    //! return true if intersection is with boundary.
    bool boundary () const {

      // Check whether we're on the left boundary
      if (neighbor_==0) {

        // If there's an element to the left we can't be on the boundary
        if (center_->pred_)
          return false;

        const OneDEntityImp<1>* ancestor = center_;

        while (ancestor->level_!=0) {

          // Check if we're the left son of our father
          if (ancestor != ancestor->father_->sons_[0])
            return false;

          ancestor = ancestor->father_;
        }

        // We have reached level 0.  If there is no element of the left
        // we're truly on the boundary
        return !ancestor->pred_;
      }

      // ////////////////////////////////
      //   Same for the right boundary
      // ////////////////////////////////
      // If there's an element to the right we can't be on the boundary
      if (center_->succ_)
        return false;

      const OneDEntityImp<1>* ancestor = center_;

      while (ancestor->level_!=0) {

        // Check if we're the left son of our father
        if (ancestor != ancestor->father_->sons_[1])
          return false;

        ancestor = ancestor->father_;
      }

      // We have reached level 0.  If there is no element of the left
      // we're truly on the boundary
      return !ancestor->succ_;

    }

    //! return true if across the edge a neighbor on this level exists
    bool neighbor () const {
      assert(neighbor_ >= 0 && neighbor_ < 2);

      return (neighbor_==0)
             ? center_->pred_ && center_->pred_->vertex_[1] == center_->vertex_[0]
             : center_->succ_ && center_->succ_->vertex_[0] == center_->vertex_[1];

    }

    //! return true if intersection is conform.
    bool conforming () const {
      return true;
    }

    //! return the Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    Entity inside() const
    {
      return Entity(OneDGridEntity<0,dim,GridImp>(center_));
    }

    //! return the Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    Entity outside() const
    {
      assert(neighbor());
      return Entity(OneDGridEntity<0,dim,GridImp>(target()));
    }

    //! return index of the boundary segment
    int boundarySegmentIndex () const {
      // It is hardwired here that the domain is connected, i.e., the boundary consists of two points
      return ((neighbor_==0 && center_->reversedBoundarySegmentNumbering_==false)
              || (neighbor_==1 && center_->reversedBoundarySegmentNumbering_==true)) ? 0 : 1;
    }

    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( LocalGeometryImpl( (indexInInside() == 0) ? 0 : 1 ) );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const
    {
      return LocalGeometry( LocalGeometryImpl( (indexInInside() == 0) ? 1 : 0 ) );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const
    {
      return Geometry( GeometryImpl( center_->vertex_[neighbor_]->pos_ ) );
    }

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return GeometryTypes::vertex;
    }

    //! local index of codim 1 entity in self where intersection is contained in
    int indexInInside () const
    {
      return neighbor_;
    }

    //! local index of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const
    {
      // If numberInSelf is 0 then numberInNeighbor is 1 and vice versa
      return 1-neighbor_;
    }

    //! return outer normal
    FieldVector<typename GridImp::ctype, dimworld> outerNormal ([[maybe_unused]] const FieldVector<typename GridImp::ctype, dim-1>&  local ) const {
      return centerUnitOuterNormal();
    }

    //! Return outer normal scaled with the integration element
    FieldVector<typename GridImp::ctype, dimworld> integrationOuterNormal ([[maybe_unused]] const FieldVector<typename GridImp::ctype, dim-1>& local ) const {
      return centerUnitOuterNormal();
    }

    //! return unit outer normal
    FieldVector<typename GridImp::ctype, dimworld> unitOuterNormal ([[maybe_unused]] const FieldVector<typename GridImp::ctype, dim-1>& local ) const {
      return centerUnitOuterNormal();
    }

    //! return unit outer normal at center of intersection
    FieldVector<typename GridImp::ctype, dimworld> centerUnitOuterNormal () const {
      return FieldVector<typename GridImp::ctype, dimworld>(2 * neighbor_ - 1);
    }

  private:
    //**********************************************************
    //  private methods
    //**********************************************************

    OneDEntityImp<1>* center_;

    /** \brief Count on which neighbor we are lookin' at.  Can be only 0 or 1. */
    int neighbor_;

  };


  /** \brief Intersection between two neighboring elements on a leaf grid */
  template<class GridImp>
  class OneDGridLeafIntersection
  {
    constexpr static int dim = GridImp::dimension;
    constexpr static int dimworld = GridImp::dimensionworld;

    // The corresponding iterator needs to access all members
    friend class OneDGridLeafIntersectionIterator<GridImp>;

    template<typename,typename>
    friend class Intersection;

    OneDGridLeafIntersection()
      : center_(nullptr)
      , neighbor_(-1) // marker for invalid intersection
    {}

    //! Constructor for a given grid entity and a given neighbor
    OneDGridLeafIntersection(OneDEntityImp<1>* center, int nb)
      : center_(center), neighbor_(nb)
    {}

    /** \brief Constructor creating the 'one-after-last'-iterator */
    OneDGridLeafIntersection(OneDEntityImp<1>* center)
      : center_(center), neighbor_(2)
    {}

    typedef typename GridImp::Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim< 1 >::LocalGeometryImpl LocalGeometryImpl;

  public:

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;

    //! equality
    bool equals(const OneDGridLeafIntersection<GridImp>& other) const {
      return (center_ == other.center_) && (neighbor_ == other.neighbor_);
    }

    OneDEntityImp<1>* target() const {
      const bool isValid = center_ && neighbor_>=0 && neighbor_<2;

      if (!isValid)
        return center_;

      if (neighbor_==0) {

        // Get left leaf neighbor
        if (center_->pred_ && center_->pred_->vertex_[1] == center_->vertex_[0]) {

          OneDEntityImp<1>* leftLeafNeighbor = center_->pred_;
          while (!leftLeafNeighbor->isLeaf()) {
            assert (leftLeafNeighbor->sons_[1] != NULL);
            leftLeafNeighbor = leftLeafNeighbor->sons_[1];
          }
          return leftLeafNeighbor;

        } else {

          OneDEntityImp<1>* ancestor = center_;
          while (ancestor->father_) {
            ancestor = ancestor->father_;
            if (ancestor->pred_ && ancestor->pred_->vertex_[1] == ancestor->vertex_[0]) {
              assert(ancestor->pred_->isLeaf());
              return ancestor->pred_;
            }
          }

          DUNE_THROW(GridError, "Programming error, apparently we're on the left boundary, neighbor_==2 should not occur!");
        }

      } else {

        // Get right leaf neighbor
        if (center_->succ_ && center_->succ_->vertex_[0] == center_->vertex_[1]) {

          OneDEntityImp<1>* rightLeafNeighbor = center_->succ_;
          while (!rightLeafNeighbor->isLeaf()) {
            assert (rightLeafNeighbor->sons_[0] != NULL);
            rightLeafNeighbor = rightLeafNeighbor->sons_[0];
          }
          return rightLeafNeighbor;

        } else {

          OneDEntityImp<1>* ancestor = center_;
          while (ancestor->father_) {
            ancestor = ancestor->father_;
            if (ancestor->succ_ && ancestor->succ_->vertex_[0] == ancestor->vertex_[1]) {
              assert(ancestor->succ_->isLeaf());
              return ancestor->succ_;
            }
          }

          DUNE_THROW(GridError, "Programming error, apparently we're on the right boundary, neighbor_==3 should not occur!");
        }

      }

    }

    //! return true if intersection is with boundary.
    bool boundary () const {

      // Check whether we're on the left boundary
      if (neighbor_==0) {

        // If there's an element to the left we can't be on the boundary
        if (center_->pred_)
          return false;

        const OneDEntityImp<1>* ancestor = center_;

        while (ancestor->level_!=0) {

          // Check if we're the left son of our father
          if (ancestor != ancestor->father_->sons_[0])
            return false;

          ancestor = ancestor->father_;
        }

        // We have reached level 0.  If there is no element of the left
        // we're truly on the boundary
        return !ancestor->pred_;
      }

      // ////////////////////////////////
      //   Same for the right boundary
      // ////////////////////////////////
      // If there's an element to the right we can't be on the boundary
      if (center_->succ_)
        return false;

      const OneDEntityImp<1>* ancestor = center_;

      while (ancestor->level_!=0) {

        // Check if we're the left son of our father
        if (ancestor != ancestor->father_->sons_[1])
          return false;

        ancestor = ancestor->father_;
      }

      // We have reached level 0.  If there is no element of the left
      // we're truly on the boundary
      return !ancestor->succ_;

    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return !boundary();
    }

    //! return true if intersection is conform.
    bool conforming () const {
      return true;
    }

    //! return Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    Entity inside() const
    {
      return Entity(OneDGridEntity<0,dim,GridImp>(center_));
    }

    //! return the Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    Entity outside() const
    {
      return Entity(OneDGridEntity<0,dim,GridImp>(target()));
    }

    //! return index of the boundary segment
    int boundarySegmentIndex () const {
      // It is hardwired here that the domain is connected, i.e., the boundary consists of two points
      return ((neighbor_==0 && center_->reversedBoundarySegmentNumbering_==false)
              || (neighbor_==1 && center_->reversedBoundarySegmentNumbering_==true)) ? 0 : 1;
    }

    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( LocalGeometryImpl( (indexInInside() == 0) ? 0 : 1 ) );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const
    {
      return LocalGeometry( LocalGeometryImpl( (indexInInside() == 0) ? 1 : 0 ) );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const
    {
      return Geometry( GeometryImpl( center_->vertex_[neighbor_%2]->pos_ ) );
    }

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return GeometryTypes::vertex;
    }

    //! local index of codim 1 entity in self where intersection is contained in
    int indexInInside () const
    {
      return neighbor_ % 2;
    }

    //! local index of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const
    {
      // If numberInSelf is 0 then numberInNeighbor is 1 and vice versa
      return 1-(neighbor_ % 2);
    }

    //! return outer normal
    FieldVector<typename GridImp::ctype, dimworld> outerNormal ([[maybe_unused]] const FieldVector<typename GridImp::ctype, dim-1>& local) const {
      return centerUnitOuterNormal();
    }

    //! Return outer normal scaled with the integration element
    FieldVector<typename GridImp::ctype, dimworld> integrationOuterNormal ([[maybe_unused]] const FieldVector<typename GridImp::ctype, dim-1>& local) const
    {
      return centerUnitOuterNormal();
    }

    //! return unit outer normal
    FieldVector<typename GridImp::ctype, dimworld> unitOuterNormal ([[maybe_unused]] const FieldVector<typename GridImp::ctype, dim-1>& local) const {
      return centerUnitOuterNormal();
    }

    //! return unit outer normal at center of intersection
    FieldVector<typename GridImp::ctype, dimworld> centerUnitOuterNormal () const {
      return FieldVector<typename GridImp::ctype, dimworld>(2 * neighbor_ - 1);
    }

  private:
    //**********************************************************
    //  private methods
    //**********************************************************

    OneDEntityImp<1>* center_;

    /** \brief Count on which neighbor we are lookin' at

       0,1 are the level neighbors, 2 and 3 are the leaf neighbors,
       if they differ from the level neighbors. */
    int neighbor_;

  };

}  // namespace Dune

#endif
