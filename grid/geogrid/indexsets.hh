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

    // Forward Declarations
    // --------------------

    template< class Grid >
    class LevelIndexSet;

    template< class Grid >
    class LeafIndexSet;

    template< class Grid, class HostIdSet >
    class IdSet;



    // LevelIndexSetTypes
    // ------------------

    template< class Grid >
    struct LevelIndexSetTypes
    {
      template< int codim >
      struct Codim
      {
        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename Grid :: template Codim< codim >
          :: template Partition< pitype > :: LevelIterator
          Iterator;
        };
      };
    };



    template< class HostGrid, class CoordFunction >
    class LevelIndexSet< const GeometryGrid< HostGrid, CoordFunction > >
      : public IndexSet
        < const GeometryGrid< HostGrid, CoordFunction >,
            LevelIndexSet< const GeometryGrid< HostGrid, CoordFunction > >,
            LevelIndexSetTypes< const GeometryGrid< HostGrid, CoordFunction > > >
    {
      typedef GeometryGrid< HostGrid, CoordFunction > Grid;

      typedef typename HostGrid :: Traits :: LevelIndexSet HostLevelIndexSet;

      const Grid *grid_;
      int level_;
      const HostLevelIndexSet *hostIndexSet_;

    public:
      enum { dimension = Grid :: dimension };

      typedef unsigned int IndexType;

      typedef IndexSet< Grid, LevelIndexSet< Grid >, LevelIndexSetTypes< Grid > >
      Base;

      LevelIndexSet ( const Grid &grid, int level )
        : grid_( &grid ),
          level_( level )
      {
        update();
      }

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
        const HostEntity &hostEntity = Grid :: template getHostEntity< codim >( entity );
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
        return grid_->template lbegin< codim, pitype >( level_ );
      }

      template< int codim, PartitionIteratorType pitype >
      typename Base :: template Codim< codim > :: template Partition< pitype > :: Iterator
      end () const
      {
        return grid_->template lend< codim, pitype >( level_ );
      }
  #endif

      void update ()
      {
        hostIndexSet_ = &(grid_->hostGrid().levelIndexSet( level_ ));
      }

    private:
      const HostLevelIndexSet &hostIndexSet () const
      {
        assert( hostIndexSet_ != 0 );
        return *hostIndexSet_;
      }
    };



    // LevelIndexSetTypes
    // ------------------

    template< class Grid >
    struct LeafIndexSetTypes
    {
      template< int codim >
      struct Codim
      {
        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename Grid :: template Codim< codim >
          :: template Partition< pitype > :: LeafIterator
          Iterator;
        };
      };
    };



    // LevelIndexSet
    // -------------

    template< class HostGrid, class CoordFunction >
    class LeafIndexSet< const GeometryGrid< HostGrid, CoordFunction > >
      : public IndexSet
        < const GeometryGrid< HostGrid, CoordFunction >,
            LeafIndexSet< const GeometryGrid< HostGrid, CoordFunction > >,
            LeafIndexSetTypes< const GeometryGrid< HostGrid, CoordFunction > > >
    {
      typedef GeometryGrid< HostGrid, CoordFunction > Grid;

      typedef typename HostGrid :: Traits :: LeafIndexSet HostLeafIndexSet;

      const Grid *grid_;
      const HostLeafIndexSet *hostIndexSet_;

    public:
      enum { dimension = Grid :: dimension };

      typedef unsigned int IndexType;

      typedef IndexSet< Grid, LeafIndexSet< Grid >, LeafIndexSetTypes< Grid > >
      Base;

      LeafIndexSet ( const Grid &grid )
        : grid_( &grid )
      {
        update();
      }

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
        const HostEntity &hostEntity = Grid :: template getHostEntity< codim >( entity );
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
        return grid_->template leafbegin< codim, pitype >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Base :: template Codim< codim > :: template Partition< pitype > :: Iterator
      end () const
      {
        return grid_->template leafend< codim, pitype >();
      }
  #endif

      void update ()
      {
        hostIndexSet_ = &(grid_->hostGrid().leafIndexSet());
      }

    private:
      const HostLeafIndexSet &hostIndexSet () const
      {
        assert( hostIndexSet_ != 0 );
        return *hostIndexSet_;
      }
    };



    // IdSet
    // -----

    template< class Grid, class HostIdSet >
    class IdSet
      : public Dune :: IdSet< Grid, IdSet< Grid, HostIdSet >, typename HostIdSet :: IdType >
    {
      typedef typename remove_const< Grid > :: type :: Traits Traits;

      const HostIdSet &hostIdSet_;

    public:
      typedef typename HostIdSet :: IdType IdType;

      IdSet ( const HostIdSet &hostIdSet )
        : hostIdSet_( hostIdSet )
      {}

      template< int codim >
      IdType id ( const typename Traits :: template Codim< codim > :: Entity &entity ) const
      {
        return Grid :: getRealImplementation( entity ).id( hostIdSet_ );
      }

      template< class Entity >
      IdType id ( const Entity &entity ) const
      {
        return id< Entity :: codimension >( entity );
      }

      template< int codim >
      IdType subId ( const typename Traits :: template Codim< 0 > :: Entity &entity, int i) const
      {
        return hostIdSet_.template subId< codim >( Grid :: template getHostEntity< 0 >( entity ), i );
      }

    private:
      IdSet ( const IdSet & );
      IdSet &operator= ( const IdSet & );
    };

  }

}

#endif
