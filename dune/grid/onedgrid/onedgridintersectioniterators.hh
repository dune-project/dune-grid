// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GRID_INTERSECTION_ITERATORS_HH
#define DUNE_ONE_D_GRID_INTERSECTION_ITERATORS_HH

/** \file
 * \brief The OneDGridLevelIntersectionIterator and OneDGridLeafIntersectionIterator classes
 */

#include <dune/grid/onedgrid/onedgridintersections.hh>

namespace Dune {

  /** \brief The iterator over level intersections */
  template<class GridImp>
  class OneDGridLevelIntersectionIterator
  {
    constexpr static int dim = GridImp::dimension;
    constexpr static int dimworld = GridImp::dimensionworld;

    friend class OneDGridEntity<0,dim,GridImp>;

    template<typename, typename, typename>
    friend class Dune::IntersectionIterator;

    OneDGridLevelIntersectionIterator()
    {}

    //! Constructor for a given grid entity and a given neighbor
    OneDGridLevelIntersectionIterator(OneDEntityImp<1>* center, int nb)
      : intersection_(OneDGridLevelIntersection<GridImp>(center,nb))
    {}

    /** \brief Constructor creating the 'one-after-last'-iterator */
    OneDGridLevelIntersectionIterator(OneDEntityImp<1>* center)
      : intersection_(OneDGridLevelIntersection<GridImp>(center))
    {}

  public:

    typedef Dune::Intersection< GridImp, Dune::OneDGridLevelIntersection< GridImp > > Intersection;

    //! equality
    bool equals(const OneDGridLevelIntersectionIterator<GridImp>& other) const {
      return intersection_.impl().equals(other.intersection_.impl());
    }

    //! prefix increment
    void increment() {
      intersection_.impl().neighbor_++;
    }

    //! \brief dereferencing
    const Intersection & dereference() const
    {
      return intersection_;
    }

  private:

    Intersection intersection_;
  };


  /** \brief Iterator over OneDGrid leaf intersections */
  template<class GridImp>
  class OneDGridLeafIntersectionIterator
  {
    constexpr static int dim = GridImp::dimension;
    constexpr static int dimworld = GridImp::dimensionworld;

    friend class OneDGridEntity<0,dim,GridImp>;

    template<typename, typename, typename>
    friend class Dune::IntersectionIterator;

    OneDGridLeafIntersectionIterator()
    {}

    //! Constructor for a given grid entity and a given neighbor
    OneDGridLeafIntersectionIterator(OneDEntityImp<1>* center, int nb)
      : intersection_(OneDGridLeafIntersection<GridImp>(center,nb))
    {}

    /** \brief Constructor creating the 'one-after-last'-iterator */
    OneDGridLeafIntersectionIterator(OneDEntityImp<1>* center)
      : intersection_(OneDGridLeafIntersection<GridImp>(center))
    {}

  public:

    typedef Dune::Intersection< GridImp, Dune::OneDGridLeafIntersection< GridImp > > Intersection;

    //! equality
    bool equals(const OneDGridLeafIntersectionIterator<GridImp>& other) const {
      return intersection_.impl().equals(other.intersection_.impl());
    }

    //! prefix increment
    void increment() {
      intersection_.impl().neighbor_++;
    }

    //! \brief dereferencing
    const Intersection & dereference() const
    {
      return intersection_;
    }

  private:

    Intersection intersection_;

  };

}  // namespace Dune

#endif
