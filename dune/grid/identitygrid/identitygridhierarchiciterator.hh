// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDHIERITERATOR_HH
#define DUNE_IDENTITYGRIDHIERITERATOR_HH

/** \file
 * \brief The IdentityGridHierarchicIterator class
 */

namespace Dune {


  //**********************************************************************
  //
  /** \brief Iterator over the descendants of an entity.
   * \ingroup IdentityGrid
     Mesh entities of codimension 0 ("elements") allow to visit all entities of
     codimension 0 obtained through nested, hierarchic refinement of the entity.
     Iteration over this set of entities is provided by the HierarchicIterator,
     starting from a given entity.
   */
  template<class GridImp>
  class IdentityGridHierarchicIterator
  {

    // Type of the corresponding HierarchicIterator in the host grid
    typedef typename GridImp::HostGridType::template Codim<0>::Entity::HierarchicIterator HostGridHierarchicIterator;

  public:

    constexpr static int codimension = 0;

    typedef typename GridImp::template Codim<0>::Entity Entity;

    //! the default Constructor
    explicit IdentityGridHierarchicIterator(const GridImp* identityGrid, const Entity& startEntity, int maxLevel) :
      identityGrid_(identityGrid),
      hostHierarchicIterator_(startEntity.impl().hostEntity_.hbegin(maxLevel))
    {}


    //! \todo Please doc me !
    explicit IdentityGridHierarchicIterator(const GridImp* identityGrid, const Entity& startEntity, int maxLevel, [[maybe_unused]] bool endDummy) :
      identityGrid_(identityGrid),
      hostHierarchicIterator_(startEntity.impl().hostEntity_.hend(maxLevel))
    {}


    //! \todo Please doc me !
    void increment()
    {
      ++hostHierarchicIterator_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity{{identityGrid_,*hostHierarchicIterator_}};
    }

    //! equality
    bool equals(const IdentityGridHierarchicIterator& i) const {
      return hostHierarchicIterator_ == i.hostHierarchicIterator_;
    }

  private:
    const GridImp* identityGrid_;

    HostGridHierarchicIterator hostHierarchicIterator_;

  };


}  // end namespace Dune

#endif
