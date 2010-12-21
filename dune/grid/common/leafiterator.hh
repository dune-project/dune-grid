// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_LEAFITERATOR_HH
#define DUNE_GRID_LEAFITERATOR_HH

#include <cstddef>
#include <iterator>

#include <dune/grid/common/entityiterator.hh>
#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  /**********************************************************************/
  /**
     @brief Enables iteration over all leaf entities
     of a codimension zero of a grid.
     See also the documentation of Dune::EntityPointer.

     \note The LeafIterator interface is deprecated. Use the EntityIterator
           interface instead.

     @ingroup GIEntityPointer
   */
  template<int codim, PartitionIteratorType pitype, class GridImp,
      template<int,PartitionIteratorType,class> class LeafIteratorImp>
  class LeafIterator
    : public EntityIterator< codim, GridImp, LeafIteratorImp< codim, pitype, GridImp > >
  {
    typedef EntityIterator< codim, GridImp, LeafIteratorImp< codim, pitype, GridImp > > Base;

  public:
    /** @brief Preincrement operator. */
    LeafIterator& operator++()
    {
      ++static_cast< Base & >( *this );
      this->realIterator.increment();
      return *this;
    }

    //===========================================================
    /** @name Implementor interface
     */
    //@{
    //===========================================================

    /** @brief copy constructor from LevelIteratorImp */
    LeafIterator (const LeafIteratorImp<codim, pitype, const GridImp> & i) DUNE_DEPRECATED
      : Base( i )
    {}
    //@}
  };

}

namespace std {

  template
  < int codim, Dune::PartitionIteratorType pitype, class GridImp,
      template<int,Dune::PartitionIteratorType,class> class LeafIteratorImp>
  struct iterator_traits<Dune::LeafIterator<codim, pitype, GridImp,
          LeafIteratorImp> > {
    typedef ptrdiff_t difference_type;
    typedef const typename Dune::LeafIterator<codim, pitype, GridImp,
        LeafIteratorImp>::Entity value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef forward_iterator_tag iterator_category;
  };

} // namespace std

#endif // DUNE_GRID_LEAFITERATOR_HH
