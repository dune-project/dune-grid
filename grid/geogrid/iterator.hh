// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ITERATOR_HH
#define DUNE_GEOGRID_ITERATOR_HH

#include <dune/grid/geogrid/entitypointer.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  // Internal Forward Declarations
  // -----------------------------

  template< int codim, PartitionIteratorType pitype, class Grid,
      bool fake = !(Capabilities :: hasHostEntity< Grid, codim > :: v) >
  class GeometryGridLeafIterator;

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLeafIteratorAdapter;

  template< int codim, PartitionIteratorType pitype, class Grid,
      bool fake = !(Capabilities :: hasHostEntity< Grid, codim > :: v ) >
  class GeometryGridLevelIterator;

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLevelIteratorAdapter;

  template< class Grid >
  class GeometryGridHierarchicIterator;



  // GeometryGridLeafIterator
  // ------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLeafIterator< codim, pitype, Grid, false >
    : public GeometryGridEntityPointer< codim, Grid, false >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;
    typedef typename Traits :: HostGrid HostGrid;

  public:
    typedef GeometryGridEntityPointer< codim, Grid, false > Base;

    typedef typename HostGrid :: template Codim< codim >
    :: template Partition< pitype > :: LeafIterator
    HostIterator;

  private:
    HostIterator hostIterator_;

  protected:
    using Base :: setToTarget;

  public:
    GeometryGridLeafIterator ( const Grid &grid,
                               const HostIterator &hostIterator )
      : Base( grid, hostIterator ),
        hostIterator_( hostIterator )
    {}

    void increment ()
    {
      ++hostIterator_;
      setToTarget( hostIterator_ );
    }
  };



  // GeometryGridLeafIteratorAdapter
  // -------------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLeafIteratorAdapter
    : public GeometryGridLeafIterator< codim, pitype, Grid >
  {
    typedef GeometryGridLeafIterator< codim, pitype, Grid > BaseType;

  public:
    typedef typename BaseType :: HostIterator HostIterator;

    GeometryGridLeafIteratorAdapter ( const Grid &grid,
                                      const HostIterator &hostIterator )
      : BaseType( grid, hostIterator )
    {}
  };



  // GeometryGridLeafIterator
  // ------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLevelIterator< codim, pitype, Grid, false >
    : public GeometryGridEntityPointer< codim, Grid, false >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;
    typedef typename Traits :: HostGrid HostGrid;

  public:
    typedef GeometryGridEntityPointer< codim, Grid, false > Base;

    typedef typename HostGrid :: template Codim< codim >
    :: template Partition< pitype > :: LevelIterator
    HostIterator;

  private:
    HostIterator hostIterator_;

  protected:
    using Base :: setToTarget;

  public:
    GeometryGridLevelIterator ( const Grid &grid,
                                const HostIterator &hostIterator )
      : Base( grid, hostIterator ),
        hostIterator_( hostIterator )
    {}

    void increment ()
    {
      ++hostIterator_;
      setToTarget( hostIterator_ );
    }
  };



  // GeometryGridLeafIteratorAdapter
  // -------------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class GeometryGridLevelIteratorAdapter
    : public GeometryGridLevelIterator< codim, pitype, Grid >
  {
    typedef GeometryGridLevelIterator< codim, pitype, Grid > BaseType;

  public:
    typedef typename BaseType :: HostIterator HostIterator;

    GeometryGridLevelIteratorAdapter ( const Grid &grid,
                                       const HostIterator &hostIterator )
      : BaseType( grid, hostIterator )
    {}
  };



  // GeometryGridHierarchicIterator
  // ------------------------------

  template< class Grid >
  class GeometryGridHierarchicIterator
    : public GeometryGridEntityPointer< 0, Grid >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;
    typedef typename Traits :: HostGrid HostGrid;

  public:
    typedef GeometryGridEntityPointer< 0, Grid > Base;

    typedef typename HostGrid :: Traits :: HierarchicIterator HostIterator;

  private:
    HostIterator hostIterator_;

  public:
    GeometryGridHierarchicIterator ( const Grid &grid,
                                     const HostIterator &hostIterator )
      : Base( grid, hostIterator ),
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
