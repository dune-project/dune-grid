// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_LEAFITERATOR_HH
#define DUNE_GEOGRID_LEAFITERATOR_HH

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



  // GeometryGridLeafIterator
  // ------------------------

  template< int codim, PartitionIteratorType pitype, class HostGrid, class CoordFunction >
  class GeometryGridLeafIterator< codim, pitype, const GeometryGrid< HostGrid, CoordFunction > >
    : public GeometryGridEntityPointer< codim, const GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    enum { dimension = Grid :: dimension };

    typedef typename HostGrid :: template Codim< codim > :: LeafIterator HostLeafIterator;

    HostLeafIterator hostIterator_;

  public:
    typedef GeometryGridEntityPointer< codim, const Grid > Base;

  protected:
    using Base :: setToTarget;

  public:
    explicit GeometryGridLeafIterator ( const Grid *grid )
      : Base( grid, grid->hostGrid().template leafbegin< codim >() ),
        hostIterator_( grid->hostGrid().template leafbegin< codim >() )
    {
      setToTarget( hostIterator_ );
    }

    GeometryGridLeafIterator( const Grid *grid, bool endDummy )
      : Base( grid, grid->hostGrid().template leafend< codim >() ),
        hostIterator_( grid->hostGrid().template leafend< codim >() )
    {}

    void increment ()
    {
      ++hostIterator_;
      setToTarget( hostIterator_ );
    }
  };

}  // namespace Dune

#endif
