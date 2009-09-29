// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INDEXSETS_HH
#define DUNE_GEOGRID_INDEXSETS_HH

#include <vector>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/indexidset.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Grid >
    class LevelIndexSet;

    template< class Grid >
    class LeafIndexSet;



    // LevelIteratorProvider
    // ---------------------

    template< class Grid >
    class LevelIteratorProvider
    {
      typedef typename remove_const< Grid > :: type :: Traits Traits;

    public:
      template< int codim >
      struct Codim
      {
        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename Traits :: template Codim< codim >
          :: template Partition< pitype > :: LevelIterator
          Iterator;
        };
      };

    private:
      const Grid &grid_;
      const int level_;

    public:
      LevelIteratorProvider ( const Grid &grid, int level )
        : grid_( grid ),
          level_( level )
      {}

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim > :: template Partition< pitype > :: Iterator
      begin () const
      {
        return grid_->template lbegin< codim, pitype >( level_ );
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim > :: template Partition< pitype > :: Iterator
      end () const
      {
        return grid_->template lend< codim, pitype >( level_ );
      }
    };



    // LevelIndexSet
    // -------------

    template< class HostGrid, class CoordFunction >
    class LevelIndexSet< const GeometryGrid< HostGrid, CoordFunction > >
      : public IndexSet
        < const GeometryGrid< HostGrid, CoordFunction >,
            LevelIndexSet< const GeometryGrid< HostGrid, CoordFunction > >,
            LevelIteratorProvider< const GeometryGrid< HostGrid, CoordFunction > > >
    {
      typedef GeometryGrid< HostGrid, CoordFunction > Grid;

      typedef typename HostGrid :: Traits :: LevelIndexSet HostIndexSet;

      typedef LevelIteratorProvider< const Grid > IteratorProvider;

      const HostIndexSet *hostIndexSet_;
      const IteratorProvider iteratorProvider_;

    public:
      enum { dimension = Grid :: dimension };

      typedef unsigned int IndexType;

      typedef IndexSet< Grid, LevelIndexSet< Grid >, IteratorProvider >
      Base;

      LevelIndexSet ( const Grid &grid, int level )
        : hostIndexSet_( &(grid.hostGrid().levelIndexSet( level )) ),
          iteratorProvider_( grid, level )
      {}

      template< int codim >
      IndexType index ( const typename Grid :: template Codim< codim > :: Entity &entity ) const
      {
        return Grid :: getRealImplementation( entity ).index( hostIndexSet() );
      }

      template< class Entity >
      IndexType index ( const Entity &entity ) const
      {
        return index< Entity :: codimension >( entity );
      }

      template< int codim >
      IndexType subIndex ( const typename Grid :: template Codim< 0 > :: Entity &entity, int i ) const
      {
        typedef typename HostGrid :: template Codim< 0 > :: Entity HostEntity;
        const HostEntity &hostEntity = Grid :: template getHostEntity< 0 >( entity );
        return hostIndexSet().template subIndex< codim >( hostEntity, i );
      }

      IndexType size ( GeometryType type ) const
      {
        return hostIndexSet().size( type );
      }

      int size ( int codim ) const
      {
        return hostIndexSet().size( codim );
      }

      template< int codim >
      bool contains ( const typename Grid :: template Codim< codim > :: Entity &entity ) const
      {
        typedef typename HostGrid :: template Codim< codim > :: Entity HostEntity;
        const HostEntity &hostEntity
          = Grid :: template getHostEntity< codim >( entity );
        return hostIndexSet().contains( hostEntity );
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return contains< Entity :: codimension >( entity );
      }

      const std :: vector< GeometryType > &geomTypes ( int codim ) const
      {
        return hostIndexSet().geomTypes( codim );
      }

#if INDEXSET_HAS_ITERATORS
      template< int codim, PartitionIteratorType pitype >
      typename Base :: template Codim< codim > :: template Partition< pitype > :: Iterator
      begin () const
      {
        return iteratorProvider_.begin();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Base :: template Codim< codim > :: template Partition< pitype > :: Iterator
      end () const
      {
        return iteratorProvider_.end();
      }
#endif

      void update ( const HostIndexSet &hostIndexSet )
      {
        hostIndexSet_ = &hostIndexSet;
      }

    private:
      const HostIndexSet &hostIndexSet () const
      {
        assert( hostIndexSet_ != 0 );
        return *hostIndexSet_;
      }
    };




    // LeafIteratorProvider
    // --------------------

    template< class Grid >
    class LeafIteratorProvider
    {
      typedef typename remove_const< Grid > :: type :: Traits Traits;

    public:
      template< int codim >
      struct Codim
      {
        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename Traits :: template Codim< codim >
          :: template Partition< pitype > :: LeafIterator
          Iterator;
        };
      };

    private:
      const Grid &grid_;

    public:
      LeafIteratorProvider ( const Grid &grid )
        : grid_( grid )
      {}

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim > :: template Partition< pitype > :: Iterator
      begin () const
      {
        return grid_->template leafbegin< codim, pitype >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim > :: template Partition< pitype > :: Iterator
      end () const
      {
        return grid_->template leafend< codim, pitype >();
      }
    };



    // LeafIndexSet
    // ------------

    template< class HostGrid, class CoordFunction >
    class LeafIndexSet< const GeometryGrid< HostGrid, CoordFunction > >
      : public IndexSet
        < const GeometryGrid< HostGrid, CoordFunction >,
            LeafIndexSet< const GeometryGrid< HostGrid, CoordFunction > >,
            LeafIteratorProvider< const GeometryGrid< HostGrid, CoordFunction > > >
    {
      typedef GeometryGrid< HostGrid, CoordFunction > Grid;

      typedef typename HostGrid :: Traits :: LeafIndexSet HostIndexSet;

      typedef LeafIteratorProvider< const Grid > IteratorProvider;

      const HostIndexSet *hostIndexSet_;
      const IteratorProvider iteratorProvider_;

    public:
      enum { dimension = Grid :: dimension };

      typedef unsigned int IndexType;

      typedef IndexSet< Grid, LeafIndexSet< Grid >, IteratorProvider >
      Base;

      LeafIndexSet ( const Grid &grid )
        : hostIndexSet_( &(grid.hostGrid().leafIndexSet()) ),
          iteratorProvider_( grid )
      {}

      template< int codim >
      IndexType index ( const typename Grid :: template Codim< codim > :: Entity &entity ) const
      {
        return Grid :: getRealImplementation( entity ).index( hostIndexSet() );
      }

      template< class Entity >
      IndexType index ( const Entity &entity ) const
      {
        return index< Entity :: codimension >( entity );
      }

      template< int codim >
      IndexType subIndex ( const typename Grid :: template Codim< 0 > :: Entity &entity, int i ) const
      {
        typedef typename HostGrid :: template Codim< 0 > :: Entity HostEntity;
        const HostEntity &hostEntity = Grid :: template getHostEntity< 0 >( entity );
        return hostIndexSet().template subIndex< codim >( hostEntity, i );
      }

      IndexType size ( GeometryType type ) const
      {
        return hostIndexSet().size( type );
      }

      int size ( int codim ) const
      {
        return hostIndexSet().size( codim );
      }

      template< int codim >
      bool contains ( const typename Grid :: template Codim< codim > :: Entity &entity ) const
      {
        typedef typename HostGrid :: template Codim< codim > :: Entity HostEntity;
        const HostEntity &hostEntity
          = Grid :: template getHostEntity< codim >( entity );
        return hostIndexSet().contains( hostEntity );
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return contains< Entity :: codimension >( entity );
      }

      const std :: vector< GeometryType > &geomTypes ( int codim ) const
      {
        return hostIndexSet().geomTypes( codim );
      }

  #if INDEXSET_HAS_ITERATORS
      template< int codim, PartitionIteratorType pitype >
      typename Base :: template Codim< codim > :: template Partition< pitype > :: Iterator
      begin () const
      {
        return iteratorProvider_.begin();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Base :: template Codim< codim > :: template Partition< pitype > :: Iterator
      end () const
      {
        return iteratorProvider.end();
      }
  #endif

      void update ( const HostIndexSet &hostIndexSet )
      {
        hostIndexSet_ = &hostIndexSet;
      }

    private:
      const HostIndexSet &hostIndexSet () const
      {
        assert( hostIndexSet_ != 0 );
        return *hostIndexSet_;
      }
    };

  }

}

#endif
