// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDHIERITERATOR_HH
#define DUNE_IDENTITYGRIDHIERITERATOR_HH

#include "identitygridentitypointer.hh"
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
  class IdentityGridHierarchicIterator :
    public Dune::IdentityGridEntityPointer<0,GridImp,typename GridImp::HostGridType::template Codim<0>::Entity::HierarchicIterator>
  {

    // Type of the corresponding HierarchicIterator in the host grid
    typedef typename GridImp::HostGridType::template Codim<0>::Entity::HierarchicIterator HostGridHierarchicIterator;

    typedef Dune::IdentityGridEntityPointer<0,GridImp,typename GridImp::HostGridType::template Codim<0>::Entity::HierarchicIterator> Base;

  public:

    typedef typename Base::Entity Entity;

    //! the default Constructor
    explicit IdentityGridHierarchicIterator(const GridImp* identityGrid, const Entity& startEntity, int maxLevel) :
      Base(identityGrid, GridImp::getRealImplementation(startEntity).hostEntity_.hbegin(maxLevel))
    {}


    //! \todo Please doc me !
    explicit IdentityGridHierarchicIterator(const GridImp* identityGrid, const Entity& startEntity, int maxLevel, bool endDummy) :
      Base(identityGrid, GridImp::getRealImplementation(startEntity).hostEntity_.hend(maxLevel))
    {}


    //! \todo Please doc me !
    void increment()
    {
      ++this->hostEntityPointer_;
    }

  };


}  // end namespace Dune

#endif
