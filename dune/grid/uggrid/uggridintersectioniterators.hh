// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGINTERSECTIONIT_HH
#define DUNE_UGINTERSECTIONIT_HH

#include <dune/common/sllist.hh>
#include <dune/common/shared_ptr.hh>
/** \file
 * \brief The UGGridIntersectionIterator class
 */

namespace Dune {

  //**********************************************************************
  //
  // --UGGridIntersectionIterator
  // --IntersectionIterator
  /** \brief Iterator over all element neighbors
   * \ingroup UGGrid
     Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
     a neighbor is an entity of codimension 0 which has a common entity of codimension 1
     These neighbors are accessed via a IntersectionIterator. This allows the implement
     non-matching meshes. The number of neigbors may be different from the number
     of an element!
   */
  template<class GridImp>
  class UGGridLevelIntersectionIterator
  {

    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    friend class UGGridEntity<0,dim,GridImp>;

    // The type used to store coordinates
    typedef typename GridImp::ctype UGCtype;

  public:

    typedef Dune::Intersection<const GridImp, Dune::UGGridLevelIntersection> Intersection;

    /** The default Constructor makes empty Iterator
     */
    UGGridLevelIntersectionIterator(typename UG_NS<dim>::Element* center, int nb)
      : intersection_(UGGridLevelIntersection<GridImp>(center,nb))
    {}

    //! equality
    bool equals(const UGGridLevelIntersectionIterator<GridImp>& other) const {
      return GridImp::getRealImplementation(intersection_).equals(GridImp::getRealImplementation(other.intersection_));
    }

    //! prefix increment
    void increment() {
      GridImp::getRealImplementation(intersection_).neighborCount_++;
    }

    //! \brief dereferencing
    const Intersection & dereference() const {
      return intersection_;
    }

  private:

    mutable MakeableInterfaceObject<Intersection> intersection_;

  };


  template<class GridImp>
  class UGGridLeafIntersectionIterator
  {

    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    friend class UGGridEntity<0,dim,GridImp>;

    // The type used to store coordinates
    typedef typename GridImp::ctype UGCtype;

  public:

    typedef Dune::Intersection<const GridImp, Dune::UGGridLeafIntersection> Intersection;

    UGGridLeafIntersectionIterator(typename UG_NS<dim>::Element* center, int nb)
      : intersection_(UGGridLeafIntersection<GridImp>(center,nb))
    {}

    //! equality
    bool equals(const UGGridLeafIntersectionIterator<GridImp>& other) const {
      return GridImp::getRealImplementation(intersection_).equals(GridImp::getRealImplementation(other.intersection_));
    }

    //! prefix increment
    void increment() {

      GridImp::getRealImplementation(intersection_).subNeighborCount_++;

      if (GridImp::getRealImplementation(intersection_).subNeighborCount_ >= GridImp::getRealImplementation(intersection_).leafSubFaces_.size() ) {

        GridImp::getRealImplementation(intersection_).neighborCount_++;
        GridImp::getRealImplementation(intersection_).subNeighborCount_ = 0;

        if (GridImp::getRealImplementation(intersection_).neighborCount_ < UG_NS<dim>::Sides_Of_Elem(GridImp::getRealImplementation(intersection_).center_))
          GridImp::getRealImplementation(intersection_).constructLeafSubfaces();

      }

    }

    //! \brief dereferencing
    const Intersection & dereference() const {
      return intersection_;
    }

  private:

    mutable MakeableInterfaceObject<Intersection> intersection_;

  };

}  // namespace Dune

#endif
