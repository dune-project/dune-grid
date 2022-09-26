// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGHIERITERATOR_HH
#define DUNE_UGHIERITERATOR_HH

/** \file
 * \brief The UGGridHierarchicIterator class
 */

#include <stack>

#include <dune/grid/uggrid/uggridentity.hh>

namespace Dune {

  //**********************************************************************
  //
  // --UGGridHierarchicIterator
  // --HierarchicIterator
  /** \brief Iterator over the descendants of an entity.
   * \ingroup UGGrid
     Mesh entities of codimension 0 ("elements") allow to visit all entities of
     codimension 0 obtained through nested, hierarchic refinement of the entity.
     Iteration over this set of entities is provided by the HierarchicIterator,
     starting from a given entity.
   */

  template<class GridImp>
  class UGGridHierarchicIterator
  {

    friend class UGGridEntity<0,GridImp::dimension,GridImp>;

  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;

    //! iterate only over elements
    constexpr static int codimension = 0;

    //! the default Constructor
    UGGridHierarchicIterator(int maxLevel, const GridImp* gridImp)
      : maxlevel_(maxLevel),
        gridImp_(gridImp)
    {}

    void increment();

    //! dereferencing
    const Entity& dereference() const {return entity_;}

    //! equality
    bool equals(const UGGridHierarchicIterator<GridImp>& other) const {
      return entity_ == other.entity_;
    }

  private:
    //! The makeable entity that the iterator is pointing to
    Entity entity_;

    //! max level to go down
    int maxlevel_;

    std::stack<typename UG_NS<GridImp::dimension>::Element*> elementStack_;

    const GridImp* gridImp_;
  };

}  // end namespace Dune

#endif
