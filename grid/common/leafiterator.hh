// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_LEAFITERATOR_HH
#define DUNE_GRID_LEAFITERATOR_HH

#include <dune/grid/common/entitypointer.hh>

namespace Dune
{

  /**********************************************************************/
  /**
     @brief Enables iteration over all leaf entities
     of a codimension zero of a grid.
     See also the documentation of Dune::EntityPointer.

     @ingroup GIEntityPointer
   */
  template<int codim, PartitionIteratorType pitype, class GridImp,
      template<int,PartitionIteratorType,class> class LeafIteratorImp>
  class LeafIterator :
    public EntityPointer<GridImp, LeafIteratorImp<codim,pitype,GridImp> >
  {
  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;

    /** @brief Preincrement operator. */
    LeafIterator& operator++()
    {
      this->realIterator.increment();
      return *this;
    }

    //===========================================================
    /** @name Implementor interface
     */
    //@{
    //===========================================================

    /** @brief copy constructor from LevelIteratorImp */
    LeafIterator (const LeafIteratorImp<codim, pitype, const GridImp> & i) :
      EntityPointer<GridImp, LeafIteratorImp<codim, pitype, GridImp> >(i) {};
    //@}
  };

}

#endif // DUNE_GRID_LEAFITERATOR_HH
