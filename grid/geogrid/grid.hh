// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_GRID_HH
#define DUNE_GEOGRID_GRID_HH

#include <string>

#include <dune/grid/common/grid.hh>

#include <dune/grid/geogrid/capabilities.hh>
#include <dune/grid/geogrid/entity.hh>
#include <dune/grid/geogrid/entitypointer.hh>
#include <dune/grid/geogrid/intersectioniterator.hh>
#include <dune/grid/geogrid/iterator.hh>
#include <dune/grid/geogrid/indexsets.hh>
#include <dune/grid/geogrid/datahandle.hh>

#include <dune/grid/genericgeometry/geometry.hh>

namespace Dune
{

  // Forward Declarations
  // --------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  // GenericGeometry :: GeometryTraits
  // ---------------------------------

  namespace GenericGeometry
  {

    template< class HostGrid, class CoordFunction >
    struct GeometryTraits< GeometryGrid< HostGrid, CoordFunction > >
      : public DefaultGeometryTraits
        < typename HostGrid :: ctype,
            HostGrid :: dimension, CoordFunction :: dimRange >
    {};

  }



  // GeometryGridFamily
  // ------------------

  template< class HostGrid, class CoordFunction >
  struct GeometryGridFamily
  {
    struct Traits
    {
      typedef GeometryGrid< HostGrid, CoordFunction > Grid;

      typedef typename HostGrid :: ctype ctype;

      dune_static_assert( (int)HostGrid :: dimensionworld == (int)CoordFunction :: dimDomain,
                          "HostGrid and CoordFunction are incompatible." );
      enum { dimension = HostGrid :: dimension };
      enum { dimensionworld = CoordFunction :: dimRange };

      typedef Intersection< const Grid, GeometryGridLeafIntersection >
      LeafIntersection;
      typedef Intersection< const Grid, GeometryGridLevelIntersection >
      LevelIntersection;

      typedef IntersectionIterator
      < const Grid, GeometryGridLeafIntersectionIterator,
          GeometryGridLeafIntersection >
      LeafIntersectionIterator;
      typedef IntersectionIterator
      < const Grid, GeometryGridLevelIntersectionIterator,
          GeometryGridLevelIntersection >
      LevelIntersectionIterator;

      typedef Dune :: HierarchicIterator
      < const Grid, GeometryGridHierarchicIterator >
      HierarchicIterator;

      template< int codim >
      struct Codim
      {
        typedef Dune :: Geometry
        < dimension-codim, dimensionworld, const Grid,
            Dune :: GenericGeometry :: Geometry >
        Geometry;
        typedef Dune :: Geometry
        < dimension-codim, dimension, const Grid,
            Dune :: GenericGeometry :: LocalGeometry >
        LocalGeometry;

        typedef Dune :: Entity
        < codim, dimension, const Grid, GeometryGridEntity >
        Entity;
        typedef Dune :: EntityPointer
        < const Grid, GeometryGridEntityPointer< codim, const Grid > >
        EntityPointer;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef Dune :: LeafIterator
          < codim, pitype, const Grid, GeometryGridLeafIterator >
          LeafIterator;
          typedef Dune :: LevelIterator
          < codim, pitype, const Grid, GeometryGridLevelIterator >
          LevelIterator;
        };

        typedef typename Partition< All_Partition > :: LeafIterator
        LeafIterator;
        typedef typename Partition< All_Partition > :: LevelIterator
        LevelIterator;
      };

      typedef GeometryGridLeafIndexSet< const Grid > LeafIndexSet;
      typedef GeometryGridLevelIndexSet< const Grid > LevelIndexSet;

      typedef GeometryGridIdSet< const Grid, typename HostGrid :: Traits :: GlobalIdSet >
      GlobalIdSet;
      typedef GeometryGridIdSet< const Grid, typename HostGrid :: Traits :: LocalIdSet >
      LocalIdSet;

      typedef typename HostGrid :: Traits :: CollectiveCommunication
      CollectiveCommunication;

      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef Dune :: GridView
        < DefaultLeafGridViewTraits< const Grid, pitype > >
        LeafGridView;
        typedef Dune :: GridView
        < DefaultLevelGridViewTraits< const Grid, pitype > >
        LevelGridView;
      };
    };
  };



  // GeometryGrid
  // ------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid
    : public GridDefaultImplementation
      < HostGrid :: dimension, CoordFunction :: dimRange, typename HostGrid :: ctype,
          GeometryGridFamily< HostGrid, CoordFunction > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    typedef GridDefaultImplementation
    < HostGrid :: dimension, CoordFunction :: dimRange, typename HostGrid :: ctype,
        GeometryGridFamily< HostGrid, CoordFunction > >
    Base;

    friend class GeometryGridLevelIndexSet< const Grid >;
    friend class GeometryGridLeafIndexSet< const Grid >;
    friend class GeometryGridHierarchicIterator< const Grid >;

    template< int, int, class > friend class GeometryGridEntity;
    template< int, class > friend class GeometryGridEntityPointer;
    template< int, PartitionIteratorType, class > friend class GeometryGridLevelIterator;
    template< int, PartitionIteratorType, class > friend class GeometryGridLeafIterator;
    template< class, class > friend class GeometryGridIntersection;
    template< class, class > friend class GeometryGridIdSet;

  public:
    typedef HostGrid HostGridType;

    typedef GeometryGridFamily< HostGrid, CoordFunction > GridFamily;
    typedef typename GridFamily :: Traits Traits;

    typedef typename Traits :: ctype ctype;

    typedef typename Traits :: HierarchicIterator HierarchicIterator;
    typedef typename Traits :: LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename Traits :: LevelIntersectionIterator LevelIntersectionIterator;

    typedef typename Traits :: LeafIndexSet LeafIndexSet;
    typedef typename Traits :: LevelIndexSet LevelIndexSet;

    typedef typename Traits :: GlobalIdSet GlobalIdSet;
    typedef typename Traits :: LocalIdSet LocalIdSet;

    typedef typename Traits :: CollectiveCommunication CollectiveCommunication;

    template< int codim >
    struct Codim
    {
      typedef typename Traits :: template Codim< codim > :: Entity Entity;
      typedef typename Traits :: template Codim< codim > :: EntityPointer EntityPointer;
      typedef typename Traits :: template Codim< codim > :: Geometry Geometry;
      typedef typename Traits :: template Codim< codim > :: LocalGeometry LocalGeometry;

      typedef typename Traits :: HierarchicIterator HierarchicIterator;
      typedef typename Traits :: LeafIntersectionIterator LeafIntersectionIterator;
      typedef typename Traits :: LevelIntersectionIterator LevelIntersectionIterator;

      typedef typename Traits :: LeafIndexSet LeafIndexSet;
      typedef typename Traits :: LevelIndexSet LevelIndexSet;

      typedef typename Traits :: GlobalIdSet GlobalIdSet;
      typedef typename Traits :: LocalIdSet LocalIdSet;

      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef typename Traits :: template Codim< codim >
        :: template Partition< pitype > :: LeafIterator
        LeafIterator;
        typedef typename Traits :: template Codim< codim >
        :: template Partition< pitype > :: LevelIterator
        LevelIterator;
      };

      typedef typename Partition< All_Partition > :: LeafIterator LeafIterator;
      typedef typename Partition< All_Partition > :: LevelIterator LevelIterator;
    };

  private:
    HostGrid *const hostGrid_;
    const CoordFunction &coordFunction_;
    mutable std :: vector< LevelIndexSet * > levelIndexSets_;
    mutable LeafIndexSet *leafIndexSet_;

    GlobalIdSet globalIdSet_;
    LocalIdSet localIdSet_;

  public:
    GeometryGrid ( HostGrid &hostGrid, const CoordFunction &coordFunction )
      : hostGrid_( &hostGrid ),
        coordFunction_( coordFunction ),
        levelIndexSets_( hostGrid.maxLevel()+1 ),
        leafIndexSet_( 0 ),
        globalIdSet_( hostGrid.globalIdSet() ),
        localIdSet_( hostGrid.localIdSet() )
    {
      for( int i = 0; i < hostGrid.maxLevel(); ++i )
        levelIndexSets_[ i ] = 0;
    }

    ~GeometryGrid ()
    {
      if( leafIndexSet_ != 0 )
        delete leafIndexSet_;

      for( unsigned int i = 0; i < levelIndexSets_.size(); ++i )
      {
        if( levelIndexSets_[ i ] != 0 )
          delete( levelIndexSets_[ i ] );
      }
    }

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    std :: string name () const
    {
      return std :: string( "GeometryGrid< " )
             + hostGrid().name() + std :: string( " >" );
    }


    //! Return maximum level defined in this grid. Levels are numbered
    //! 0 ... maxlevel with 0 the coarsest level.
    int maxLevel () const
    {
      return hostGrid().maxLevel();
    }

    template< int codim >
    typename Codim< codim > :: LevelIterator lbegin ( int level ) const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim > :: LevelIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      Impl impl( *this, hostGrid().template lbegin< codim >( level ) );
      return MakeableIterator( impl );
    }

    template< int codim >
    typename Codim< codim > :: LevelIterator lend ( int level ) const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim > :: LevelIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      Impl impl( *this, hostGrid().template lend< codim >( level ) );
      return MakeableIterator( impl );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition<pitype> :: LevelIterator
    lbegin ( int level ) const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim >
          :: template Partition< pitype > :: LevelIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      Impl impl( *this, hostGrid().template lbegin< codim, pitype >( level ) );
      return MakeableIterator( impl );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition<pitype> :: LevelIterator
    lend ( int level ) const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim >
          :: template Partition< pitype > :: LevelIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      Impl impl( *this, hostGrid().template lend< codim, pitype >( level ) );
      return MakeableIterator( impl );
    }

    template< int codim >
    typename Codim< codim > :: LeafIterator leafbegin () const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim > :: LeafIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      Impl impl( *this, hostGrid().template leafbegin< codim >() );
      return MakeableIterator( impl );
    }

    template< int codim >
    typename Codim< codim > :: LeafIterator leafend () const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim > :: LeafIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      Impl impl( *this, hostGrid().template leafend< codim >() );
      return MakeableIterator( impl );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: LeafIterator
    leafbegin () const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim >
          :: template Partition< pitype > :: LeafIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      Impl impl( *this, hostGrid().template leafbegin< codim, pitype >() );
      return MakeableIterator( impl );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: LeafIterator
    leafend () const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim >
          :: template Partition< pitype > :: LeafIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      Impl impl( *this, hostGrid().template leafend< codim, pitype >() );
      return MakeableIterator( impl );
    }

    int size ( int level, int codim ) const
    {
      return hostGrid().size( level, codim );
    }

    int size ( int codim ) const
    {
      return hostGrid().size( codim );
    }

    int size ( int level, GeometryType type ) const
    {
      return hostGrid().size( level, type );
    }

    int size ( GeometryType type ) const
    {
      return hostGrid().size( type );
    }

    const GlobalIdSet &globalIdSet () const
    {
      return globalIdSet_;
    }

    const LocalIdSet &localIdSet () const
    {
      return localIdSet_;
    }

    const LevelIndexSet &levelIndexSet ( int level ) const
    {
      if( (level < 0) || (level > maxLevel()) )
        DUNE_THROW( GridError, "levelIndexSet of nonexisting level " << level << " requested." );
      if( levelIndexSets_[ level ] == 0 )
        levelIndexSets_[ level ] = new LevelIndexSet( *this, level );
      return *levelIndexSets_[ level ];
    }

    const LeafIndexSet &leafIndexSet () const
    {
      if( leafIndexSet_ == 0 )
        leafIndexSet_ = new LeafIndexSet( *this );
      return *leafIndexSet_;
    }

    void globalRefine ( int refCount )
    {
      hostGrid().globalRefine( refCount );
      updateIndexSets();
    }

    bool mark( int refCount, const typename Codim< 0 > :: EntityPointer &entity )
    {
      return hostGrid().mark( refCount, getHostEntityPointer< 0 >( entity ) );
    }

    int getMark ( const typename Codim< 0 > :: EntityPointer &entity ) const
    {
      return hostGrid().getMark( getHostEntityPointer< 0 >( entity ) );
    }

    bool preAdapt ()
    {
      return hostGrid().preAdapt();
    }

    bool adapt ()
    {
      bool ret = hostGrid().adapt();
      updateIndexSets();
      return ret;
    }

    void postAdapt ()
    {
      return hostGrid().postAdapt();
    }

    unsigned int overlapSize ( int codim ) const
    {
      return hostGrid().overlapSize( codim );
    }

    unsigned int ghostSize( int codim ) const
    {
      return hostGrid().ghostSize( codim );
    }


    /** \brief Size of the overlap on a given level */
    unsigned int overlapSize ( int level, int codim ) const
    {
      return hostGrid().overlapSize( level, codim );
    }


    /** \brief Size of the ghost cell layer on a given level */
    unsigned int ghostSize ( int level, int codim ) const
    {
      return hostGrid().ghostSize( level, codim );
    }


    bool loadBalance ()
    {
      bool ret = hostGrid().loalBalance();
      updateIndexSets();
      return ret;
    }

    template< class DataHandle, class DataType >
    bool loadBalance ( CommDataHandleIF< DataHandle, DataType > &data )
    {
      typedef CommDataHandleIF< DataHandle, DataType > DataHandleIF;
      typedef GeometryGridCommDataHandle< Grid, DataHandleIF > WrappedDataHandle;

      WrappedDataHandle wrappedData( *this, data );
      bool ret = hostGrid().loadBalance( wrappedData );
      updateIndexSets();
      return ret;
    }

    template< class DataHandle, class DataType >
    void communicate ( CommDataHandleIF< DataHandle, DataType > &data,
                       InterfaceType iftype, CommunicationDirection dir,
                       int level ) const
    {
      typedef CommDataHandleIF< DataHandle, DataType > DataHandleIF;
      typedef GeometryGridCommDataHandle< Grid, DataHandleIF > WrappedDataHandle;

      WrappedDataHandle wrappedData( *this, data );
      hostGrid().communicate( wrappedData, iftype, dir, level );
    }

    template< class DataHandle, class DataType >
    void communicate ( CommDataHandleIF< DataHandle, DataType > &data,
                       InterfaceType iftype, CommunicationDirection dir ) const
    {
      typedef CommDataHandleIF< DataHandle, DataType > DataHandleIF;
      typedef GeometryGridCommDataHandle< Grid, DataHandleIF > WrappedDataHandle;

      WrappedDataHandle wrappedData( *this, data );
      hostGrid().communicate( wrappedData, iftype, dir );
    }

    const CollectiveCommunication &comm () const
    {
      return hostGrid().comm();
    }


    // **********************************************************
    // End of Interface Methods
    // **********************************************************

#if 0
    HostGrid &getHostGrid() const
    {
      return *hostGrid_;
    }
#endif

  private:
    //! Obtain the host grid wrapped by this GeometryGrid
    const HostGrid &hostGrid () const
    {
      return *hostGrid_;
    }

    //! Obtain the host grid wrapped by this GeometryGrid
    HostGrid &hostGrid ()
    {
      return *hostGrid_;
    }

    //! Obtain the coordinate function
    const CoordFunction &coordFunction () const
    {
      return coordFunction_;
    }

  protected:
    using Base :: getRealImplementation;

    template< int codim >
    static const typename HostGrid :: template Codim< codim > :: Entity &
    getHostEntity( const typename Codim< codim > :: Entity &entity )
    {
      return getRealImplementation( entity ).hostEntity();
    }

    template< int codim >
    static const typename HostGrid :: template Codim< codim > :: EntityPointer &
    getHostEntityPointer( const typename Codim< codim > :: EntityPointer &entity )
    {
      return getRealImplementation( entity ).hostEntityPointer();
    }

  private:
    void updateIndexSets ()
    {
      if( leafIndexSet_ != 0 )
        leafIndexSet_->update();

      const int newNumLevels = maxLevel()+1;
      const int oldNumLevels = levelIndexSets_.size();
      int updateLevels = std :: min( oldNumLevels, newNumLevels );

      for( int i = 0; i < updateLevels; ++i )
      {
        if( levelIndexSets_[ i ] != 0 )
          levelIndexSets_[ i ]->update();
      }

      for( int i = updateLevels; i < oldNumLevels; ++i )
      {
        if( levelIndexSets_[ i ] != 0 )
          delete levelIndexSets_[ i ];
      }

      levelIndexSets_.resize( newNumLevels );
      for( int i = updateLevels; i < newNumLevels; ++i )
        levelIndexSets_[ i ] = 0;
    }
  };

}

#endif
