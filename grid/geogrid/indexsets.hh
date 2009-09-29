// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_INDEXSETS_HH
#define DUNE_GEOGRID_INDEXSETS_HH

#include <vector>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  // Forward Declarations
  // --------------------

  template< class Grid >
  class GeometryGridLevelIndexSet;

  template< class Grid >
  class GeometryGridLeafIndexSet;

  template< class Grid, class HostIdSet >
  class GeometryGridIdSet;



  // GeometryGridLevelIndexSetTypes
  // ------------------------------

  template< class Grid >
  struct GeometryGridLevelIndexSetTypes
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
  class GeometryGridLevelIndexSet< const GeometryGrid< HostGrid, CoordFunction > >
    : public IndexSetDefaultImplementation
      < const GeometryGrid< HostGrid, CoordFunction >,
          GeometryGridLevelIndexSet< const GeometryGrid< HostGrid, CoordFunction > >,
          GeometryGridLevelIndexSetTypes< const GeometryGrid< HostGrid, CoordFunction > > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    typedef typename HostGrid :: Traits :: LevelIndexSet HostLevelIndexSet;

    const Grid *grid_;
    const HostLevelIndexSet *hostIndexSet_;

  public:
    enum { dimension = Grid :: dimension };

    typedef unsigned int IndexType;

    typedef IndexSet
    < Grid, GeometryGridLevelIndexSet< Grid >, GeometryGridLevelIndexSetTypes< Grid > >
    Base;

    GeometryGridLevelIndexSet ( const Grid &grid, int level )
    {
      update( grid, level );
    }

    template< int codim >
    IndexType index ( const typename Grid :: template Codim< codim > :: Entity &entity ) const
    {
      typedef typename HostGrid :: template Codim< codim > :: EntityPointer HostEntityPointer;
      HostEntityPointer hostEntity = Grid :: template getHostEntity< codim >( entity );
      return hostIndexSet().template index< codim >( *hostEntity );
    }

    template< class Entity >
    IndexType index ( const Entity &entity ) const
    {
      return index< Entity :: codimension >( entity );
    }

    template< int codim >
    IndexType subIndex ( const typename Grid :: template Codim< 0 > :: Entity &entity, int i ) const
    {
      typedef typename HostGrid :: template Codim< 0 > :: EntityPointer HostEntityPointer;
      HostEntityPointer hostEntity = Grid :: template getHostEntity< 0 >( entity );
      return hostIndexSet().template subIndex< codim >( *hostEntity, i );
    }

    IndexType size ( GeometryType type ) const
    {
      return hostIndexSet().size( type );
    }

    int size ( int codim ) const
    {
      return hostIndexSet().size( codim );
    }

    const std :: vector< GeometryType > &geomTypes ( int codim ) const
    {
      return hostIndexSet().geomTypes( codim );
    }

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

    void update ( const Grid &grid, int level )
    {
      grid_ = &grid;
      hostIndexSet_ = &(grid.hostGrid().levelIndexSet( level ));
    }

  private:
    const HostLevelIndexSet &hostIndexSet () const
    {
      return *hostIndexSet_;
    }
  };



  // GeometryGridLevelIndexSetTypes
  // ------------------------------

  template< class Grid >
  struct GeometryGridLeafIndexSetTypes
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



  // GeometryGridLevelIndexSet
  // -------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGridLeafIndexSet< const GeometryGrid< HostGrid, CoordFunction > >
    : public IndexSetDefaultImplementation
      < const GeometryGrid< HostGrid, CoordFunction >,
          GeometryGridLeafIndexSet< const GeometryGrid< HostGrid, CoordFunction > >,
          GeometryGridLeafIndexSetTypes< const GeometryGrid< HostGrid, CoordFunction > > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    typedef typename HostGrid :: Traits :: LeafIndexSet HostLeafIndexSet;

    const Grid *grid_;
    const HostLeafIndexSet *hostIndexSet_;

  public:
    enum { dimension = Grid :: dimension };

    typedef unsigned int IndexType;

    typedef IndexSet
    < Grid, GeometryGridLeafIndexSet< Grid >, GeometryGridLeafIndexSetTypes< Grid > >
    Base;

    GeometryGridLeafIndexSet ( const Grid &grid )
    {
      update( grid );
    }

    template< int codim >
    IndexType index ( const typename Grid :: template Codim< codim > :: Entity &entity ) const
    {
      typedef typename HostGrid :: template Codim< codim > :: EntityPointer HostEntityPointer;
      HostEntityPointer hostEntity = Grid :: template getHostEntity< codim >( entity );
      return hostIndexSet().template index< codim >( *hostEntity );
    }

    template< class Entity >
    IndexType index ( const Entity &entity ) const
    {
      return index< Entity :: codimension >( entity );
    }

    template< int codim >
    IndexType subIndex ( const typename Grid :: template Codim< 0 > :: Entity &entity, int i ) const
    {
      typedef typename HostGrid :: template Codim< 0 > :: EntityPointer HostEntityPointer;
      HostEntityPointer hostEntity = Grid :: template getHostEntity< 0 >( entity );
      return hostIndexSet().template subIndex< codim >( *hostEntity, i );
    }

    IndexType size ( GeometryType type ) const
    {
      return hostIndexSet().size( type );
    }

    int size ( int codim ) const
    {
      return hostIndexSet().size( codim );
    }

    const std :: vector< GeometryType > &geomTypes ( int codim ) const
    {
      return hostIndexSet().geomTypes( codim );
    }

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

    void update ( const Grid &grid )
    {
      grid_ = &grid;
      hostIndexSet_ = &(grid.hostGrid().leafIndexSet());
    }

  private:
    const HostLeafIndexSet &hostIndexSet () const
    {
      return *hostIndexSet_;
    }
  };



  // GeometryGridIdSet
  // -----------------

  template< class Grid, class HostIdSet >
  class GeometryGridIdSet
    : public IdSet
      < Grid, GeometryGridIdSet< Grid, HostIdSet >, typename HostIdSet :: IdType >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

    const HostIdSet &hostIdSet_;

  public:
    typedef typename HostIdSet :: IdType IdType;

    GeometryGridIdSet ( const HostIdSet &hostIdSet )
      : hostIdSet_( hostIdSet )
    {}

    template< int codim >
    IdType id ( const typename Traits :: template Codim< codim > :: Entity &entity ) const
    {
      return hostIdSet_.id( *(Grid :: template getHostEntity< codim >( entity )) );
    }

    template< class Entity >
    IdType id ( const Entity &entity ) const
    {
      return id< Entity :: codimension >( entity );
    }

    template< int codim >
    IdType subId ( const typename Traits :: template Codim< 0 > :: Entity &entity, int i) const
    {
      return hostIdSet_.subId( *(Grid :: template getHostEntity< 0 >( entity )), i );
    }

  private:
    GeometryGridIdSet ( const GeometryGridIdSet & );
    GeometryGridIdSet &operator= ( const GeometryGridIdSet & );
  };

}

#endif
