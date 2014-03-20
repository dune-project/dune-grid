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
    public Dune::IdentityGridEntityPointer <0,GridImp>
  {
  public:

    typedef typename GridImp::template Codim<0>::Entity Entity;

    typedef IdentityGridEntity <0, GridImp::dimension, GridImp> IdentityGridElement;


    //! the default Constructor
    explicit IdentityGridHierarchicIterator(const GridImp* identityGrid, const IdentityGridElement& startEntity, int maxLevel) :
      IdentityGridEntityPointer<0,GridImp>(identityGrid, startEntity.hostEntity_->hbegin(maxLevel)),
      identityGrid_(identityGrid),
      hostGridHierarchicIterator_(startEntity.hostEntity_->hbegin(maxLevel)),
      hostGridHierarchicEndIterator_(startEntity.hostEntity_->hend(maxLevel))
    {
      this->virtualEntity_.setToTarget(hostGridHierarchicIterator_);
    }


    //! \todo Please doc me !
    explicit IdentityGridHierarchicIterator(const GridImp* identityGrid, const IdentityGridElement& startEntity, int maxLevel, bool endDummy) :
      IdentityGridEntityPointer<0,GridImp>(identityGrid, startEntity.hostEntity_->hend(maxLevel)),
      identityGrid_(identityGrid),
      hostGridHierarchicIterator_(startEntity.hostEntity_->hbegin(maxLevel)),
      hostGridHierarchicEndIterator_(startEntity.hostEntity_->hend(maxLevel))
    {}


    //! \todo Please doc me !
    void increment()
    {
      ++hostGridHierarchicIterator_;
      this->virtualEntity_.setToTarget(hostGridHierarchicIterator_);
    }


  private:

    // Type of the corresponding HierarchicIterator in the host grid
    typedef typename GridImp::HostGridType::template Codim<0>::Entity::HierarchicIterator HostGridHierarchicIterator;

    enum {dim = GridImp::HostGridType::dimension};


    // The level index of the host entity that we are pointing to
    //! \todo Please doc me !
    unsigned int hostLevelIndex() const {
      return identityGrid_->hostgrid_->levelIndexSet(hostGridHierarchicIterator_.level()).index(*hostGridHierarchicIterator_);
    }


    const GridImp* identityGrid_;

    HostGridHierarchicIterator hostGridHierarchicIterator_;

    HostGridHierarchicIterator hostGridHierarchicEndIterator_;
  };


}  // end namespace Dune

#endif
