// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGHIERITERATOR_HH
#define DUNE_UGHIERITERATOR_HH

/** \file
 * \brief The UGGridHierarchicIterator class
 */

#include <stack>

#include <dune/grid/uggrid/uggridentitypointer.hh>

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
  class UGGridHierarchicIterator :
    public Dune::UGGridEntityPointer <0,GridImp>
  {

    friend class UGGridEntity<0,GridImp::dimension,GridImp>;

  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;

    //! the default Constructor
    UGGridHierarchicIterator(int maxLevel)
      : maxlevel_(maxLevel)
    {}

    void increment();

    //! max level to go down
    int maxlevel_;

    std::stack<typename UG_NS<GridImp::dimension>::Element*> elementStack_;

  };

  // Include class method definitions
#include "uggridhieriterator.cc"

}  // end namespace Dune

#endif
