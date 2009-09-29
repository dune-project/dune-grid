// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ITERATOR_HH
#define DUNE_GEOGRID_ITERATOR_HH

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  // Internal Forward Declarations
  // -----------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLeafIterator;

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLevelIterator;



  // GeometryGridLeafIterator
  // ------------------------

  template< int codim, PartitionIteratorType pitype,
      class HostGrid, class CoordFunction >
  class GeometryGridLeafIterator
  < codim, pitype, const GeometryGrid< HostGrid, CoordFunction > >
    : public GeometryGridEntityPointer
      < codim, const GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    typedef typename HostGrid :: template Codim< codim > :: LeafIterator
    HostIterator;

    HostIterator hostIterator_;

  public:
    typedef GeometryGridEntityPointer< codim, const Grid > Base;

  protected:
    using Base :: setToTarget;

  public:
    GeometryGridLeafIterator ( const Grid &grid,
                               const HostIterator &hostIterator )
      : Base( &grid, hostIterator ),
        hostIterator_( hostIterator )
    {}

    void increment ()
    {
      ++hostIterator_;
      setToTarget( hostIterator_ );
    }
  };



  // GeometryGridLeafIterator
  // ------------------------

  template< int codim, PartitionIteratorType pitype,
      class HostGrid, class CoordFunction >
  class GeometryGridLevelIterator
  < codim, pitype, const GeometryGrid< HostGrid, CoordFunction > >
    : public GeometryGridEntityPointer
      < codim, const GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    typedef typename HostGrid :: template Codim< codim > :: LevelIterator
    HostIterator;

    HostIterator hostIterator_;

  public:
    typedef GeometryGridEntityPointer< codim, const Grid > Base;

  protected:
    using Base :: setToTarget;

  public:
    GeometryGridLevelIterator ( const Grid &grid,
                                const HostIterator &hostIterator )
      : Base( &grid, hostIterator ),
        hostIterator_( hostIterator )
    {}

    void increment ()
    {
      ++hostIterator_;
      setToTarget( hostIterator_ );
    }
  };

}

#endif
