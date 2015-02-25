// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDLEVELITERATOR_HH
#define DUNE_IDENTITYGRIDLEVELITERATOR_HH

#include "identitygridentitypointer.hh"

/** \file
 * \brief The IdentityGridLevelIterator class and its specializations
 */

namespace Dune {




  //**********************************************************************
  //
  // --IdentityGridLevelIterator
  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup IdentityGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class IdentityGridLevelIterator :
    public Dune::IdentityGridEntityPointer<codim,GridImp,typename GridImp::HostGridType::Traits::template Codim<codim>::template Partition<pitype>::LevelIterator>
  {

    typedef typename GridImp::HostGridType::Traits::template Codim<codim>::template Partition<pitype>::LevelIterator HostGridLevelIterator;
    typedef Dune::IdentityGridEntityPointer<codim,GridImp,HostGridLevelIterator> Base;

  public:

    //! Constructor
    explicit IdentityGridLevelIterator(const GridImp* identityGrid, int level)
      : Base(identityGrid, identityGrid->hostgrid_->levelGridView(level).template begin<codim,pitype>())
    {}


    /** \brief Constructor which create the end iterator
        \param endDummy      Here only to distinguish it from the other constructor
        \param identityGrid  pointer to IdentityGrid instance
        \param level         grid level on which the iterator shall be created
     */
    explicit IdentityGridLevelIterator(const GridImp* identityGrid, int level, bool endDummy)
      : Base(identityGrid, identityGrid->hostgrid_->levelGridView(level).template end<codim,pitype>())
    {}


    //! prefix increment
    void increment() {
      ++this->hostEntityPointer_;
    }

  };


}  // namespace Dune

#endif
