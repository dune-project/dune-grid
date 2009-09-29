// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_GRID_HH
#define DUNE_GEOGRID_GRID_HH

#include <string>

#include <dune/common/static_assert.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid/geogrid/capabilities.hh>
#include <dune/grid/geogrid/entity.hh>
#include <dune/grid/geogrid/entitypointer.hh>
#include <dune/grid/geogrid/intersection.hh>
#include <dune/grid/geogrid/intersectioniterator.hh>
#include <dune/grid/geogrid/iterator.hh>
#include <dune/grid/geogrid/idset.hh>
#include <dune/grid/geogrid/indexsets.hh>
#include <dune/grid/geogrid/datahandle.hh>
#include <dune/grid/geogrid/backuprestore.hh>

#include <dune/grid/genericgeometry/geometry.hh>

#include <dune/grid/geogrid/identity.hh>

/** \brief namespace shared by all DUNE modules
 *
 *  This is DUNE's main name space. All DUNE modules provide their
 *  functionality though this name space (or a subspace of it).
 */
namespace Dune
{

  // Forward Declarations
  // --------------------

  template< class HostGrid >
  class DefaultCoordFunction;

  template< class HostGrid,
      class CoordFunction = DefaultCoordFunction< HostGrid > >
  class GeometryGrid;



  // DefaultCoordFunction
  // --------------------

  template< class HostGrid >
  class DefaultCoordFunction
    : public IdenticalCoordFunction
      < typename HostGrid :: ctype, HostGrid :: dimensionworld >
  {};



  // GenericGeometry :: GeometryTraits
  // ---------------------------------

  namespace GenericGeometry
  {

    template< class HostGrid, class CoordFunction >
    struct GlobalGeometryTraits< GeometryGrid< HostGrid, CoordFunction > >
      : public DefaultGeometryTraits
        < typename HostGrid :: ctype, HostGrid :: dimension, CoordFunction :: dimRange >
    {
      typedef GeometryGrid< HostGrid, CoordFunction > Grid;

      typedef DuneCoordTraits< typename HostGrid :: ctype > CoordTraits;

      static const int dimGrid = HostGrid :: dimension;
      static const int dimWorld = CoordFunction :: dimRange;

      static const bool hybrid = true;
      // static const GeometryType :: BasicType dunetype = GeometryType :: simplex;

      static const GeometryType :: BasicType linetype = GeometryType :: simplex;

      template< class Topology >
      struct Mapping
      {
        typedef GeoGrid :: CornerStorage< Topology, const Grid > CornerStorage;
        typedef CornerMapping< CoordTraits, Topology, dimWorld, CornerStorage > type;
      };

      struct Caching
      {
        static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
        static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
        static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
        static const EvaluationType evaluateNormal = ComputeOnDemand;
      };
    };

  }



  /** \brief namespace containing the implementations of GeometryGrid */
  namespace GeoGrid
  {

    // GeometryGridExportParams
    // ------------------------

    template< class HG, class CF >
    class ExportParams
    {
      static const bool isCoordFunction
        = isCoordFunctionInterface< typename CF :: Interface > :: value;
      dune_static_assert( isCoordFunction, "Invalid CoordFunction." );

    public:
      typedef HG HostGrid;
      typedef CF CoordFunction;
    };



    // GeometryGridFamily
    // ------------------

    template< class HostGrid, class CoordFunction >
    struct GridFamily
    {
      struct Traits
        : public ExportParams< HostGrid, CoordFunction >
      {
        typedef GeometryGrid< HostGrid, CoordFunction > Grid;

        typedef typename HostGrid :: ctype ctype;

        static const int dimension = HostGrid :: dimension;
        static const int dimensionworld = CoordFunction :: dimRange;

        typedef Dune :: Intersection< const Grid, GeoGrid :: LeafIntersection >
        LeafIntersection;
        typedef Dune :: Intersection< const Grid, GeoGrid :: LevelIntersection >
        LevelIntersection;

        typedef Dune :: IntersectionIterator
        < const Grid, GeoGrid :: LeafIntersectionIterator,
            GeoGrid :: LeafIntersection >
        LeafIntersectionIterator;
        typedef Dune :: IntersectionIterator
        < const Grid, GeoGrid :: LevelIntersectionIterator,
            GeoGrid :: LevelIntersection >
        LevelIntersectionIterator;

        typedef Dune :: HierarchicIterator
        < const Grid, GeoGrid :: HierarchicIterator >
        HierarchicIterator;

        template< int codim >
        struct Codim
        {
          typedef Dune :: Geometry
          < dimension-codim, dimensionworld, const Grid,
              Dune :: GenericGeometry :: Geometry >
          Geometry;
          typedef typename HostGrid :: template Codim< codim > :: LocalGeometry
          LocalGeometry;

          typedef GeoGrid :: EntityPointerTraits< codim, const Grid >
          EntityPointerTraits;
          typedef GeoGrid :: EntityPointer< EntityPointerTraits > EntityPointerImpl;
          typedef Dune :: EntityPointer< const Grid, EntityPointerImpl > EntityPointer;
          typedef typename EntityPointerTraits :: Entity Entity;

          template< PartitionIteratorType pitype >
          struct Partition
          {
            typedef Dune::LeafIterator
            < codim, pitype, const Grid, GeoGrid::LeafIterator >
            LeafIterator;
            typedef Dune::LevelIterator
            < codim, pitype, const Grid, GeoGrid::LevelIterator >
            LevelIterator;
          };

          typedef typename Partition< All_Partition >::LeafIterator
          LeafIterator;
          typedef typename Partition< All_Partition >::LevelIterator
          LevelIterator;
        };

        typedef GeoGrid::IndexSet
        < const Grid, typename HostGrid::Traits::LeafIndexSet >
        LeafIndexSet;
        typedef GeoGrid::IndexSet
        < const Grid, typename HostGrid::Traits::LevelIndexSet >
        LevelIndexSet;

        typedef GeoGrid::IdSet
        < const Grid, typename HostGrid::Traits::GlobalIdSet >
        GlobalIdSet;
        typedef GeoGrid::IdSet
        < const Grid, typename HostGrid::Traits::LocalIdSet >
        LocalIdSet;

        typedef typename HostGrid::Traits::CollectiveCommunication
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

  }



  // GeometryGrid
  // ------------

  /** \class GeometryGrid
   *  \brief grid wrapper replacing the geometries
   *
   *  GeometryGrid wraps another DUNE grid and replaces its geometry by the
   *  generic geometries from dune-grid. These are linear (respectively
   *  n-linear) DUNE geometries interpolating some given corners. These corners
   *  are obtained by mapping the corners of the host grid's geometry (also
   *  called host geometry) by a coordinate function.
   *
   *  An example of a coordinate function is given by the following code:
   *  \code
   *  class ExampleFunction
   *  : public Dune :: AnalyticalCoordFunction< double, 2, 3, ExampleFunction >
   *  {
   *    typedef ExampleFunction This;
   *    typedef Dune :: AnalyticalCoordFunction< double, 2, 3, This > Base;
   *
   *  public:
   *    typedef Base :: DomainVector DomainVector;
   *    typedef Base :: RangeVector RangeVector;
   *
   *    void evaluate ( const DomainVector &x, RangeVector &y ) const
   *    {
   *      y[ 0 ] = x[ 0 ];
   *      y[ 1 ] = x[ 1 ];
   *      y[ 2 ] = x[ 0 ] + x[ 1 ];
   *    }
   *  };
   *  \endcode
   *
   *  \note A dune-fem Function can be used as a coordinate function.
   *        The evaluation of discrete functions would be very expensive,
   *        though.
   *
   *  \tparam HostGrid       DUNE grid to be wrapped (called host grid)
   *  \tparam CoordFunction  coordinate function
   *
   *  \nosubgrouping
   */
  template< class HostGrid, class CoordFunction >
  class GeometryGrid
  /** \cond */
    : public GridDefaultImplementation
      < HostGrid :: dimension, CoordFunction :: dimRange,
          typename HostGrid :: ctype,
          GeoGrid :: GridFamily< HostGrid, CoordFunction > >,
      public GeoGrid :: ExportParams< HostGrid, CoordFunction >,
      public GeoGrid :: HierarchicIndexSetProvider< HostGrid, CoordFunction >,
      public GeoGrid :: BackupRestoreFacilities
      < GeometryGrid< HostGrid, CoordFunction >,
          Capabilities :: hasBackupRestoreFacilities< GeometryGrid< HostGrid, CoordFunction > > :: v >
      /** \endcond */
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    typedef GridDefaultImplementation
    < HostGrid :: dimension, CoordFunction :: dimRange,
        typename HostGrid :: ctype,
        GeoGrid :: GridFamily< HostGrid, CoordFunction > >
    Base;

    typedef GeoGrid :: HierarchicIndexSetProvider< HostGrid, CoordFunction >
    HierarchicIndexSetProvider;

    friend class GeoGrid :: HierarchicIterator< const Grid >;

    template< int, class, bool > friend class GeoGrid :: EntityBase;
    template< class, bool > friend class GeoGrid :: EntityPointer;
    template< class, class > friend class GeoGrid :: Intersection;
    template< class, class > friend class GeoGrid :: IdSet;
    template< class, class, class > friend class GeoGrid :: IndexSet;
    template< class > friend class HostGridAccess;

    template< int, PartitionIteratorType, class >
    friend class GeoGrid :: LevelIteratorTraits;
    template< int, PartitionIteratorType, class >
    friend class GeoGrid :: LeafIteratorTraits;

  public:
    /** \cond */
    typedef GeoGrid :: GridFamily< HostGrid, CoordFunction > GridFamily;
    /** \endcond */

    /** \name Traits
     *  \{ */

    //! type of the grid traits
    typedef typename GridFamily :: Traits Traits;

    /** \brief traits structure containing types for a codimension
     *
     *  \tparam codim  codimension
     *
     *  \nosubgrouping
     */
    template< int codim >
    struct Codim;

    /** \} */

    /** \name Iterator Types
     *  \{ */

    //! iterator over the grid hierarchy
    typedef typename Traits :: HierarchicIterator HierarchicIterator;
    //! iterator over intersections with other entities on the leaf level
    typedef typename Traits :: LeafIntersectionIterator LeafIntersectionIterator;
    //! iterator over intersections with other entities on the same level
    typedef typename Traits :: LevelIntersectionIterator LevelIntersectionIterator;

    /** \} */

    /** \name Index and Id Set Types
     *  \{ */

    /** \brief type of leaf index set
     *
     *  The index set assigns consecutive indices to the entities of the
     *  leaf grid. The indices are of integral type and can be used to access
     *  arrays.
     *
     *  The leaf index set is a model of Dune::IndexSet.
     */
    typedef typename Traits :: LeafIndexSet LeafIndexSet;

    /** \brief type of level index set
     *
     *  The index set assigns consecutive indices to the entities of a grid
     *  level. The indices are of integral type and can be used to access
     *  arrays.
     *
     *  The level index set is a model of Dune::IndexSet.
     */
    typedef typename Traits :: LevelIndexSet LevelIndexSet;

    /** \brief type of global id set
     *
     *  The id set assigns a unique identifier to each entity within the
     *  grid. This identifier is unique over all processes sharing this grid.
     *
     *  \note Id's are neither consecutive nor necessarily of an integral
     *        type.
     *
     *  The global id set is a model of Dune::IdSet.
     */
    typedef typename Traits :: GlobalIdSet GlobalIdSet;

    /** \brief type of local id set
     *
     *  The id set assigns a unique identifier to each entity within the
     *  grid. This identifier needs only to be unique over this process.
     *
     *  Though the local id set may be identical to the global id set, it is
     *  often implemented more efficiently.
     *
     *  \note Ids are neither consecutive nor necessarily of an integral
     *        type.
     *  \note Local ids need not be compatible with global ids. Also, no
     *        mapping from local ids to global ones needs to exist.
     *
     *  The global id set is a model of Dune::IdSet.
     */
    typedef typename Traits :: LocalIdSet LocalIdSet;

    /** \} */

    /** \name Miscellaneous Types
     * \{ */

    //! type of vector coordinates (e.g., double)
    typedef typename Traits :: ctype ctype;

    //! communicator with all other processes having some part of the grid
    typedef typename Traits :: CollectiveCommunication CollectiveCommunication;

    /** \} */

  private:
    HostGrid *const hostGrid_;
    CoordFunction &coordFunction_;
    mutable std :: vector< LevelIndexSet * > levelIndexSets_;
    mutable LeafIndexSet *leafIndexSet_;

    GlobalIdSet globalIdSet_;
    LocalIdSet localIdSet_;

  public:
    /** \name Construction and Destruction
     *  \{ */

    /** \brief constructor
     *
     *  The references to host grid and coordinate function are stored in the
     *  grid. Therefore, they must remain valid until the grid is destroyed.
     *
     *  \param hostGrid       reference to the grid to wrap
     *  \param coordFunction  reference to the coordinate function
     */
    GeometryGrid ( HostGrid &hostGrid, CoordFunction &coordFunction )
      : HierarchicIndexSetProvider( *this ),
        hostGrid_( &hostGrid ),
        coordFunction_( coordFunction ),
        levelIndexSets_( hostGrid.maxLevel()+1, (LevelIndexSet *) 0 ),
        leafIndexSet_( 0 ),
        globalIdSet_( hostGrid.globalIdSet() ),
        localIdSet_( hostGrid.localIdSet() )
    {}

    /** \brief destructor
     */
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

    /** \} */

    /** \name Grid Identification Methods
     *  \{ */

    /** \brief obtain a string naming the grid
     *
     *  \returns ''GeometryGrid\< \em host \em grid \em name \>''
     */
    std :: string name () const
    {
      return std :: string( "GeometryGrid< " )
             + hostGrid().name() + std :: string( " >" );
    }

    /** \} */


    /** \name Size Methods
     *  \{ */

    /** \brief obtain maximal grid level
     *
     *  Grid levels are numbered 0, ..., L, where L is the value returned by
     *  this method.
     *
     *  \returns maximal grid level
     */
    int maxLevel () const
    {
      return hostGrid().maxLevel();
    }

    /** \brief obtain number of entites on a level
     *
     *  \param[in]  level  level to consider
     *  \param[in]  codim  codimension to consider
     *
     *  \returns number of entities of codimension \em codim on grid level
     *           \em level.
     */
    int size ( int level, int codim ) const
    {
      return hostGrid().size( level, codim );
    }

    /** \brief obtain number of leaf entities
     *
     *  \param[in]  codim  codimension to consider
     *
     *  \returns number of leaf entities of codimension \em codim
     */
    int size ( int codim ) const
    {
      return hostGrid().size( codim );
    }

    /** \brief obtain number of entites on a level
     *
     *  \param[in]  level  level to consider
     *  \param[in]  type   geometry type to consider
     *
     *  \returns number of entities with a geometry of type \em type on grid
     *           level \em level.
     */
    int size ( int level, GeometryType type ) const
    {
      return hostGrid().size( level, type );
    }

    /** \brief obtain number of leaf entities
     *
     *  \param[in]  type   geometry type to consider
     *
     *  \returns number of leaf entities with a geometry of type \em type
     */
    int size ( GeometryType type ) const
    {
      return hostGrid().size( type );
    }

    /** \} */

    template< int codim >
    typename Codim< codim > :: LevelIterator lbegin ( int level ) const
    {
      typedef MakeableInterfaceObject< typename Codim< codim > :: LevelIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      return MakeableIterator( Impl( *this, level, Impl :: Traits :: begin ) );
    }

    template< int codim >
    typename Codim< codim > :: LevelIterator lend ( int level ) const
    {
      typedef MakeableInterfaceObject< typename Codim< codim > :: LevelIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      return MakeableIterator( Impl( *this, level, Impl :: Traits :: end ) );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition<pitype> :: LevelIterator
    lbegin ( int level ) const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim > :: template Partition< pitype > :: LevelIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      return MakeableIterator( Impl( *this, level, Impl :: Traits :: begin ) );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition<pitype> :: LevelIterator
    lend ( int level ) const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim > :: template Partition< pitype > :: LevelIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      return MakeableIterator( Impl( *this, level, Impl :: Traits :: end ) );
    }

    template< int codim >
    typename Codim< codim > :: LeafIterator leafbegin () const
    {
      typedef MakeableInterfaceObject< typename Codim< codim > :: LeafIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      return MakeableIterator( Impl( *this, Impl :: Traits :: begin ) );
    }

    template< int codim >
    typename Codim< codim > :: LeafIterator leafend () const
    {
      typedef MakeableInterfaceObject< typename Codim< codim > :: LeafIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      return MakeableIterator( Impl( *this, Impl :: Traits :: end ) );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: LeafIterator
    leafbegin () const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim > :: template Partition< pitype > :: LeafIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      return MakeableIterator( Impl( *this, Impl :: Traits :: begin ) );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: LeafIterator
    leafend () const
    {
      typedef MakeableInterfaceObject
      < typename Codim< codim > :: template Partition< pitype > :: LeafIterator >
      MakeableIterator;
      typedef typename MakeableIterator :: ImplementationType Impl;
      return MakeableIterator( Impl( *this, Impl :: Traits :: end ) );
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
      assert( levelIndexSets_.size() == (size_t)(maxLevel()+1) );
      if( (level < 0) || (level > maxLevel()) )
      {
        DUNE_THROW( GridError, "LevelIndexSet for nonexisting level " << level
                                                                      << " requested." );
      }

      if( levelIndexSets_[ level ] == 0 )
      {
        levelIndexSets_[ level ]
          = new LevelIndexSet( hostGrid().levelIndexSet( level ) );
      }
      assert( levelIndexSets_[ level ] != 0 );
      return *levelIndexSets_[ level ];
    }

    const LeafIndexSet &leafIndexSet () const
    {
      if( leafIndexSet_ == 0 )
        leafIndexSet_ = new LeafIndexSet( hostGrid().leafIndexSet() );
      assert( leafIndexSet_ != 0 );
      return *leafIndexSet_;
    }

    void globalRefine ( int refCount )
    {
      hostGrid().globalRefine( refCount );
      update();
    }

    bool mark ( int refCount, const typename Codim< 0 > :: Entity &entity )
    {
      return hostGrid().mark( refCount, getHostEntity< 0 >( entity ) );
    }

    int getMark ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return hostGrid().getMark( getHostEntity< 0 >( entity ) );
    }

    bool preAdapt ()
    {
      return hostGrid().preAdapt();
    }

    bool adapt ()
    {
      bool ret = hostGrid().adapt();
      update();
      return ret;
    }

    void postAdapt ()
    {
      return hostGrid().postAdapt();
    }

    /** \name Parallel Data Distribution and Communication Methods
     *  \{ */

    /** \brief obtain size of overlap region for the leaf grid
     *
     *  \param[in]  codim  codimension for with the information is desired
     */
    int overlapSize ( int codim ) const
    {
      return hostGrid().overlapSize( codim );
    }

    /** \brief obtain size of ghost region for the leaf grid
     *
     *  \param[in]  codim  codimension for with the information is desired
     */
    int ghostSize( int codim ) const
    {
      return hostGrid().ghostSize( codim );
    }

    /** \brief obtain size of overlap region for a grid level
     *
     *  \param[in]  level  grid level (0, ..., maxLevel())
     *  \param[in]  codim  codimension (0, ..., dimension)
     */
    int overlapSize ( int level, int codim ) const
    {
      return hostGrid().overlapSize( level, codim );
    }

    /** \brief obtain size of ghost region for a grid level
     *
     *  \param[in]  level  grid level (0, ..., maxLevel())
     *  \param[in]  codim  codimension (0, ..., dimension)
     */
    int ghostSize ( int level, int codim ) const
    {
      return hostGrid().ghostSize( level, codim );
    }

    /** \brief communicate information on a grid level
     *
     *  \param      datahandle  communication data handle (user defined)
     *  \param[in]  interface   communication interface (one of
     *                          InteriorBorder_InteriorBorder_Interface,
     *                          InteriorBorder_All_Interface,
     *                          Overlap_OverlapFront_Interface,
     *                          Overlap_All_Interface,
     *                          All_All_Interface)
     *  \param[in]  direction   communication direction (one of
     *                          ForwardCommunication or BackwardCommunication)
     *  \param[in]  level       grid level to communicate
     */
    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &datahandle,
                       InterfaceType interface,
                       CommunicationDirection direction,
                       int level ) const
    {
      typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
      typedef GeoGrid :: CommDataHandle< Grid, DataHandleIF > WrappedDataHandle;

      WrappedDataHandle wrappedDataHandle( *this, datahandle );
      hostGrid().communicate( wrappedDataHandle, interface, direction, level );
    }

    /** \brief communicate information on leaf entities
     *
     *  \param      datahandle  communication data handle (user defined)
     *  \param[in]  interface   communication interface (one of
     *                          InteriorBorder_InteriorBorder_Interface,
     *                          InteriorBorder_All_Interface,
     *                          Overlap_OverlapFront_Interface,
     *                          Overlap_All_Interface,
     *                          All_All_Interface)
     *  \param[in]  direction   communication direction (one of
     *                          ForwardCommunication, BackwardCommunication)
     */
    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &datahandle,
                       InterfaceType interface,
                       CommunicationDirection direction ) const
    {
      typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
      typedef GeoGrid :: CommDataHandle< Grid, DataHandleIF > WrappedDataHandle;

      WrappedDataHandle wrappedDataHandle( *this, datahandle );
      hostGrid().communicate( wrappedDataHandle, interface, direction );
    }

    /** \brief obtain CollectiveCommunication object
     *
     *  The CollectiveCommunication object should be used to globally
     *  communicate information between all processes sharing this grid.
     *
     *  \note The CollectiveCommunication object returned is identical to the
     *        one returned by the host grid.
     */
    const CollectiveCommunication &comm () const
    {
      return hostGrid().comm();
    }

    /** \brief rebalance the load each process has to handle
     *
     *  A parallel grid is redistributed such that each process has about
     *  the same load (e.g., the same number of leaf entites).
     *
     *  \note DUNE does not specify, how the load is measured.
     *
     *  \returns \b true, if the grid has changed.
     */
    bool loadBalance ()
    {
      const bool gridChanged= hostGrid().loadBalance();
      if( gridChanged )
        update();
      return gridChanged;
    }

    /** \brief rebalance the load each process has to handle
     *
     *  A parallel grid is redistributed such that each process has about
     *  the same load (e.g., the same number of leaf entites).
     *
     *  The data handle is used to communicate the data associated with
     *  entities that move from one process to another.
     *
     *  \note DUNE does not specify, how the load is measured.
     *
     *  \param  datahandle  communication data handle (user defined)
     *
     *  \returns \b true, if the grid has changed.
     */
    template< class DataHandle, class Data >
    bool loadBalance ( CommDataHandleIF< DataHandle, Data > &datahandle )
    {
      typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
      typedef GeoGrid :: CommDataHandle< Grid, DataHandleIF > WrappedDataHandle;

      WrappedDataHandle wrappedDataHandle( *this, datahandle );
      const bool gridChanged = hostGrid().loadBalance( wrappedDataHandle );
      if( gridChanged )
        update();
      return gridChanged;
    }

    /** \} */

    /** \name Miscellaneous Methods
     *  \{ */

    const HostGrid &hostGrid () const
    {
      return *hostGrid_;
    }

    HostGrid &hostGrid ()
    {
      return *hostGrid_;
    }

    /** \brief update grid caches
     *
     *  This method has to be called whenever the underlying host grid changes.
     *
     *  \note If you adapt the host grid through this geometry grid's
     *        adaptation or load balancing methods, update is automatically
     *        called.
     */
    void update ()
    {
      // adapt the coordinate function
      GeoGrid :: AdaptCoordFunction< typename CoordFunction :: Interface >
      :: adapt( coordFunction_ );

      const int newNumLevels = maxLevel()+1;
      const int oldNumLevels = levelIndexSets_.size();

      for( int i = newNumLevels; i < oldNumLevels; ++i )
      {
        if( levelIndexSets_[ i ] != 0 )
          delete levelIndexSets_[ i ];
      }
      levelIndexSets_.resize( newNumLevels, (LevelIndexSet *)0 );
    }

    /** \} */

  protected:
    using Base :: getRealImplementation;

    const CoordFunction &coordFunction () const
    {
      return coordFunction_;
    }

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
  };



  template< class HostGrid, class CoordFunction >
  template< int codim >
  struct GeometryGrid< HostGrid, CoordFunction > :: Codim
    : public Base :: template Codim< codim >
  {
    /** \name Entity and Entity Pointer Types
     *  \{ */

    /** \brief type of entity
     *
     *  The entity is a model of Dune::Entity.
     */
    typedef typename Traits :: template Codim< codim > :: Entity Entity;

    /** \brief type of entity pointer
     *
     *  The entity pointer is a model of Dune::EntityPointer.
     */
    typedef typename Traits :: template Codim< codim > :: EntityPointer EntityPointer;

    /** \} */

    /** \name Geometry Types
     *  \{ */

    /** \brief type of world geometry
     *
     *  Models the geomtry mapping of the entity, i.e., the mapping from the
     *  reference element into world coordinates.
     *
     *  The geometry is a model of Dune::Geometry, implemented through the
     *  generic geometries provided by dune-grid.
     */
    typedef typename Traits :: template Codim< codim > :: Geometry Geometry;

    /** \brief type of local geometry
     *
     *  Models the geomtry mapping into the reference element of dimension
     *  \em dimension.
     *
     *  The local geometry is a model of Dune::Geometry, implemented through
     *  the generic geometries provided by dune-grid.
     */
    typedef typename Traits :: template Codim< codim > :: LocalGeometry LocalGeometry;

    /** \} */

    /** \name Iterator Types
     *  \{ */

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

    /** \brief type of level iterator
     *
     *  This iterator enumerates the entites of codimension \em codim of a
     *  grid level.
     *
     *  The level iterator is a model of Dune::LevelIterator.
     */
    typedef typename Partition< All_Partition > :: LeafIterator LeafIterator;

    /** \brief type of leaf iterator
     *
     *  This iterator enumerates the entites of codimension \em codim of the
     *  leaf grid.
     *
     *  The leaf iterator is a model of Dune::LeafIterator.
     */
    typedef typename Partition< All_Partition > :: LevelIterator LevelIterator;

    /** \} */
  };


}

#endif
