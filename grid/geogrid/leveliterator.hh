// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_LEVELITERATOR_HH
#define DUNE_GEOGRID_LEVELITERATOR_HH

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  // Internal Forward Declarations
  // -----------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLevelIterator;



  // GeometryGridLeafIterator
  // ------------------------

  template< int codim, PartitionIteratorType pitype, class HostGrid, class CoordFunction >
  class GeometryGridLevelIterator< codim, pitype, const GeometryGrid< HostGrid, CoordFunction > >
    : public GeometryGridEntityPointer< codim, const GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    enum { dimension = Grid :: dimension };

    typedef typename HostGrid :: template Codim< codim > :: LevelIterator HostLevelIterator;

    HostLevelIterator hostIterator_;

  public:
    typedef GeometryGridEntityPointer< codim, const Grid > Base;

    explicit GeometryGridLevelIterator ( const Grid *grid, int level )
      : Base( grid, grid->hostGrid().template lbegin< codim >( level ) ),
        hostIterator_( grid->hostGrid().template lbegin< codim >( level ) )
    {
      this->virtualEntity_.setToTarget( hostIterator_ );
    }

    GeometryGridLevelIterator( const Grid *grid, int level, bool endDummy )
      : Base( grid, grid->hostGrid().template lend< codim >( level ) ),
        hostIterator_( grid->hostGrid().template lend< codim >( level ) )
    {}

    void increment ()
    {
      ++hostIterator_;
      this->virtualEntity_.setToTarget( hostIterator_ );
    }
  };

}  // namespace Dune

#endif
