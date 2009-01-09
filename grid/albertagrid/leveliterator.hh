// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_LEVELITERATOR_HH
#define DUNE_ALBERTA_LEVELITERATOR_HH

#include <dune/grid/common/leveliterator.hh>

#include <dune/grid/albertagrid/treeiterator.hh>

namespace Dune
{

  // AlbertaGridLevelIterator
  // ------------------------

  //! the same as TreeIterator
  template< int codim, PartitionIteratorType pitype, class GridImp >
  class AlbertaGridLevelIterator
    : public AlbertaGridTreeIterator< codim, GridImp, false >
  {
    typedef AlbertaGridLevelIterator< codim, pitype, GridImp > This;
    typedef AlbertaGridTreeIterator< codim, GridImp, false > Base;

  public:
    typedef typename Base::Entity Entity;
    typedef typename Base::MarkerVector MarkerVector;

    //! Constructor making end iterator
    AlbertaGridLevelIterator ( const GridImp &grid, int level )
      : Base( grid, level )
    {}

    //! Constructor making begin iterator
    AlbertaGridLevelIterator ( const GridImp &grid,
                               const MarkerVector *vec,
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
