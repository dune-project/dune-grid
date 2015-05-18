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

    typedef typename GridImp::HostGridType::LeafGridView::IntersectionIterator HostLeafIntersectionIterator;

  public:

    typedef Dune::Intersection<const GridImp, Dune::IdentityGridLeafIntersection<GridImp> > Intersection;

    IdentityGridLeafIntersectionIterator()
    {}

    IdentityGridLeafIntersectionIterator(const GridImp* identityGrid,
                                         const HostLeafIntersectionIterator& hostIterator)
      : identityGrid_(identityGrid)
      , hostIterator_(hostIterator)
    {}

    //! equality
    bool equals(const IdentityGridLeafIntersectionIterator& other) const {
      return hostIterator_ == other.hostIterator_;
    }


    //! prefix increment
    void increment() {
      ++hostIterator_;
    }

    //! \brief dereferencing
    Intersection dereference() const {
      return IdentityGridLeafIntersection<GridImp>(identityGrid_,*hostIterator_);
    }

  private:
    //**********************************************************
    //  private data
    //**********************************************************

    const GridImp* identityGrid_;
    HostLeafIntersectionIterator hostIterator_;
  };




  //! \todo Please doc me !
  template<class GridImp>
  class IdentityGridLevelIntersectionIterator
  {
    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::HostGridType::LevelGridView::IntersectionIterator HostLevelIntersectionIterator;

  public:

    typedef Dune::Intersection<const GridImp, Dune::IdentityGridLevelIntersection<GridImp> > Intersection;

    IdentityGridLevelIntersectionIterator()
    {}

    IdentityGridLevelIntersectionIterator(const GridImp* identityGrid,
                                          const HostLevelIntersectionIterator& hostIterator)
      : identityGrid_(identityGrid)
      , hostIterator_(hostIterator)
    {}

    //! equality
    bool equals(const IdentityGridLevelIntersectionIterator<GridImp>& other) const {
      return hostIterator_ == other.hostIterator_;
    }


    //! prefix increment
    void increment() {
      ++hostIterator_;
    }

    //! \brief dereferencing
    Intersection dereference() const {
      return IdentityGridLevelIntersection<GridImp>(identityGrid_,*hostIterator_);
    }

  private:


    const GridImp* identityGrid_;
    HostLevelIntersectionIterator hostIterator_;

  };


}  // namespace Dune

#endif
