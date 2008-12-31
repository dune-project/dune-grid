// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_LEAFITERATOR_HH
#define DUNE_ALBERTA_LEAFITERATOR_HH

#include <dune/grid/common/leafiterator.hh>

#include <dune/grid/albertagrid/treeiterator.hh>

namespace Dune
{

  // AlbertaGridLeafIterator
  // -----------------------

  //! LeafIterator which is just a hull for the LevelIterator
  template< int codim, PartitionIteratorType pitype, class GridImp >
  class AlbertaGridLeafIterator
    : public AlbertaGridTreeIterator< codim, GridImp, true >
  {
    typedef AlbertaGridLeafIterator< codim, pitype, GridImp > This;
    typedef AlbertaGridTreeIterator< codim, GridImp, true > Base;

  public:
    typedef typename GridImp::template Codim< codim >::Entity Entity;

    //! Constructor making end iterator
    AlbertaGridLeafIterator ( const GridImp &grid, int level )
      : Base( grid, level )
    {}

    //! Constructor making begin iterator
    AlbertaGridLeafIterator ( const GridImp &grid,
                              const AlbertaMarkerVector *vec,
                              int level )
      : Base( grid, vec, level )
    {}

    //! increment the iterator
    void increment ()
    {
      Base::increment();
    }
  };

}

#endif
