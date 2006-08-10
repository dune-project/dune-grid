// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRIDLEVELITERATOR_HH
#define DUNE_UGGRIDLEVELITERATOR_HH

#include <dune/grid/uggrid/uggridentitypointer.hh>

/** \file
 * \brief The UGGridLevelIterator class
 */

namespace Dune {

  //**********************************************************************
  //
  // --UGGridLevelIterator
  // --LevelIterator
  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup UGGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class UGGridLevelIterator :
    public Dune::UGGridEntityPointer <codim,GridImp>,
    public LevelIteratorDefaultImplementation <codim,pitype,GridImp,UGGridLevelIterator>
  {
    enum {dim = GridImp::dimension};

    friend class UGGridEntity<codim,GridImp::dimension,GridImp>;
    friend class UGGridEntity<0,    GridImp::dimension,GridImp>;

  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor
    explicit UGGridLevelIterator()
    {
      this->virtualEntity_.setToTarget(NULL);
    }

    //! Constructor
    explicit UGGridLevelIterator(typename UG_NS<dim>::template Entity<codim>::T* target)
    {
      this->virtualEntity_.setToTarget(target);
    }

    //! prefix increment
    void increment() {
      assert(this->level() == UG_NS<dim>::myLevel(this->virtualEntity_.getTarget()));
      this->virtualEntity_.setToTarget(UG_NS<dim>::succ(this->virtualEntity_.getTarget()));
    }

  };

}  // namespace Dune

#endif
