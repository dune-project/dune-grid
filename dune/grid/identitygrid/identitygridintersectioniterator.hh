// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRID_INTERSECTIONITERATOR_HH
#define DUNE_IDENTITYGRID_INTERSECTIONITERATOR_HH

#include "identitygridintersections.hh"
#include "identitygridentity.hh"

#include <dune/grid/common/intersection.hh>

/** \file
 * \brief The IdentityGridLeafIntersectionIterator and IdentityGridLevelIntersectionIterator classes
 */

namespace Dune {

  /** \brief Iterator over all element neighbors
   * \ingroup IdentityGrid
   * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
   * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
   * These neighbors are accessed via a IntersectionIterator. This allows the implement
   * non-matching meshes. The number of neighbors may be different from the number
   * of an element!
   */
  template<class GridImp>
  class IdentityGridLeafIntersectionIterator
  {

    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::HostGridType::template Codim<0>::Entity::LeafIntersectionIterator HostLeafIntersectionIterator;

  public:

    typedef Dune::Intersection<const GridImp, Dune::IdentityGridLeafIntersection<GridImp> > Intersection;

    IdentityGridLeafIntersectionIterator(const GridImp* identityGrid,
                                         const HostLeafIntersectionIterator& hostIterator)
      : intersection_(IdentityGridLeafIntersection<GridImp>(identityGrid, hostIterator))
    {}

    //! The Destructor
    ~IdentityGridLeafIntersectionIterator() {};

    //! equality
    bool equals(const IdentityGridLeafIntersectionIterator<GridImp>& other) const {
      return GridImp::getRealImplementation(intersection_).hostIterator_
             == GridImp::getRealImplementation(other.intersection_).hostIterator_;
    }


    //! prefix increment
    void increment() {
      ++GridImp::getRealImplementation(intersection_).hostIterator_;
    }

    //! \brief dereferencing
    const Intersection & dereference() const {
      return intersection_;
    }

  private:
    //**********************************************************
    //  private data
    //**********************************************************

    /** \brief The actual intersection
     */
    mutable MakeableInterfaceObject<Intersection> intersection_;
  };




  //! \todo Please doc me !
  template<class GridImp>
  class IdentityGridLevelIntersectionIterator
  {
    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::HostGridType::template Codim<0>::Entity::LevelIntersectionIterator HostLevelIntersectionIterator;

  public:

    typedef Dune::Intersection<const GridImp, Dune::IdentityGridLevelIntersection<GridImp> > Intersection;

    IdentityGridLevelIntersectionIterator(const GridImp* identityGrid,
                                          const HostLevelIntersectionIterator& hostIterator)
      : intersection_(IdentityGridLevelIntersection<GridImp>(identityGrid,hostIterator))
    {}

    //! equality
    bool equals(const IdentityGridLevelIntersectionIterator<GridImp>& other) const {
      return GridImp::getRealImplementation(intersection_).hostIterator_ == GridImp::getRealImplementation(other.intersection_).hostIterator_;
    }


    //! prefix increment
    void increment() {
      ++GridImp::getRealImplementation(intersection_).hostIterator_;
    }

    //! \brief dereferencing
    const Intersection & dereference() const {
      return intersection_;
    }

  private:

    /** \brief The actual intersection
     */
    mutable MakeableInterfaceObject<Intersection> intersection_;

  };


}  // namespace Dune

#endif
