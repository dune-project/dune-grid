// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_GRID_HH
#define DUNE_GEOGRID_GRID_HH

#include <dune/common/nullptr.hh>
#include <dune/common/deprecated.hh>

#include <dune/grid/common/grid.hh>

#include <dune/grid/geometrygrid/backuprestore.hh>
#include <dune/grid/geometrygrid/capabilities.hh>
#include <dune/grid/geometrygrid/datahandle.hh>
#include <dune/grid/geometrygrid/gridfamily.hh>
#include <dune/grid/geometrygrid/identity.hh>
#include <dune/grid/geometrygrid/persistentcontainer.hh>

namespace Dune
{

  // DefaultCoordFunction
  // --------------------

  template< class HostGrid >
  class DefaultCoordFunction
    : public IdenticalCoordFunction< typename HostGrid::ctype, HostGrid::dimensionworld >
  {};



  // GeometryGrid
  // ------------

  /** \class GeometryGrid
   *  \brief grid wrapper replacing the geometries
   *  \ingroup GeoGrid
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
  template< class HostGrid, class CoordFunction = DefaultCoordFunction< HostGrid >, class Allocator = std::allocator< void > >
  class GeometryGrid
  /** \cond */
    : public GridDefaultImplementation
      < HostGrid::dimension, CoordFunction::dimRange, typename HostGrid::ctype,
          GeoGrid::GridFamily< HostGrid, CoordFunction, Allocator > >,
      public GeoGrid::ExportParams< HostGrid, CoordFunction >,
      public GeoGrid::BackupRestoreFacilities< GeometryGrid< HostGrid, CoordFunction, Allocator > >
      /** \endcond */
  {
    typedef GeometryGrid< HostGrid, CoordFunction, Allocator > Grid;

    typedef GridDefaultImplementation
    < HostGrid::dimension, CoordFunction::dimRange, typename HostGrid::ctype,
        GeoGrid::GridFamily< HostGrid, CoordFunction, Allocator > >
    Base;

    friend class GeoGrid::HierarchicIterator< const Grid >;

    template< int, class, bool > friend class GeoGrid::EntityBase;
    template< int, int, class > friend class GeoGrid::Geometry;
    template< class, class, class, PartitionIteratorType > friend class GeoGrid::GridView;
    template< class, class > friend class GeoGrid::Intersection;
    template< class, class > friend class GeoGrid::IntersectionIterator;
    template< class, class > friend class GeoGrid::IdSet;
    template< class, class > friend class GeoGrid::IndexSet;
    template< class > friend struct HostGridAccess;

    template< class, class > friend class GeoGrid::CommDataHandle;

  public:
    /** \cond */
    typedef GeoGrid::GridFamily< HostGrid, CoordFunction, Allocator > GridFamily;
    /** \endcond */

    /** \name Traits
     *  \{ */

    //! type of the grid traits
    typedef typename GridFamily::Traits Traits;

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
    typedef typename Traits::HierarchicIterator HierarchicIterator;
    //! iterator over intersections with other entities on the leaf level
    typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
    //! iterator over intersections with other entities on the same level
    typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

    /** \} */

    /** \name Grid View Types
     *  \{ */

    /** \brief Types for GridView */
    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename GridFamily::Traits::template Partition< pitype >::LevelGridView
      LevelGridView;
      typedef typename GridFamily::Traits::template Partition< pitype >::LeafGridView
      LeafGridView;
    };

    /** \brief View types for All_Partition */
    typedef typename Partition< All_Partition >::LevelGridView LevelGridView;
    typedef typename Partition< All_Partition >::LeafGridView LeafGridView;

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
    typedef typename Traits::LeafIndexSet LeafIndexSet;

    /** \brief type of level index set
     *
     *  The index set assigns consecutive indices to the entities of a grid
     *  level. The indices are of integral type and can be used to access
     *  arrays.
     *
     *  The level index set is a model of Dune::IndexSet.
     */
    typedef typename Traits::LevelIndexSet LevelIndexSet;

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
    typedef typename Traits::GlobalIdSet GlobalIdSet;

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
    typedef typename Traits::LocalIdSet LocalIdSet;

    /** \} */

    /** \name Miscellaneous Types
     * \{ */

    //! type of vector coordinates (e.g., double)
    typedef typename Traits::ctype ctype;

    //! communicator with all other processes having some part of the grid
    typedef typename Traits::CollectiveCommunication CollectiveCommunication;

    /** \} */

    /** \name Construction and Destruction
     *  \{ */

    /** \brief constructor
     *
     *  The references to host grid and coordinate function are stored in the
     *  grid. Therefore, they must remain valid until the grid is destroyed.
     *
     *  \param[in]  hostGrid       reference to the grid to wrap
     *  \param[in]  coordFunction  reference to the coordinate function
     *  \param[in]  allocator      storage allocator
     */
    GeometryGrid ( HostGrid &hostGrid, CoordFunction &coordFunction, const Allocator &allocator = Allocator() )
      : hostGrid_( &hostGrid ),
        coordFunction_( coordFunction ),
        removeHostGrid_( false ),
        levelIndexSets_( hostGrid_->maxLevel()+1, nullptr, allocator ),
        storageAllocator_( allocator )
    {}

    /** \brief constructor
     *
     *  The grid takes ownership of the pointers to host grid and coordinate
     *  function. They will be deleted when the grid is destroyed.
     *
     *  \param[in]  hostGrid       pointer to the grid to wrap
     *  \param[in]  coordFunction  pointer to the coordinate function
     *  \param[in]  allocator      storage allocator
     */
    GeometryGrid ( HostGrid *hostGrid, CoordFunction *coordFunction, const Allocator &allocator = Allocator() )
      : hostGrid_( hostGrid ),
        coordFunction_( *coordFunction ),
        removeHostGrid_( true ),
        levelIndexSets_( hostGrid_->maxLevel()+1, nullptr, allocator ),
        storageAllocator_( allocator )
    {}

    /** \brief destructor
     */
    ~GeometryGrid ()
    {
      for( unsigned int i = 0; i < levelIndexSets_.size(); ++i )
      {
        if( levelIndexSets_[ i ] )
          delete( levelIndexSets_[ i ] );
      }

      if( removeHostGrid_ )
      {
        delete &coordFunction_;
        delete hostGrid_;
      }
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
      return levelGridView( level ).size( codim );
    }

    /** \brief obtain number of leaf entities
     *
     *  \param[in]  codim  codimension to consider
     *
     *  \returns number of leaf entities of codimension \em codim
     */
    int size ( int codim ) const
    {
      return leafGridView().size( codim );
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
      return levelGridView( level ).size( type );
    }

    /** \brief obtain number of leaf entities
     *
     *  \returns number of leaf entities with a geometry of type \em type
     */
    int size ( GeometryType type ) const
    {
      return leafGridView().size( type );
    }

    /** \brief returns the number of boundary segments within the macro grid
     *
     *  \returns number of boundary segments within the macro grid
     */
    size_t numBoundarySegments () const
    {
      return hostGrid().numBoundarySegments( );
    }
    /** \} */

    const GlobalIdSet &globalIdSet () const
    {
      if( !globalIdSet_ )
        globalIdSet_ = GlobalIdSet( hostGrid().globalIdSet() );
      assert( globalIdSet_ );
      return globalIdSet_;
    }

    const LocalIdSet &localIdSet () const
    {
      if( !localIdSet_ )
        localIdSet_ = LocalIdSet( hostGrid().localIdSet() );
      assert( localIdSet_ );
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

      LevelIndexSet *&levelIndexSet = levelIndexSets_[ level ];
      if( !levelIndexSet )
        levelIndexSet = new LevelIndexSet( hostGrid().levelIndexSet( level ) );
      assert( levelIndexSet );
      return *levelIndexSet;
    }

    const LeafIndexSet &leafIndexSet () const
    {
      if( !leafIndexSet_ )
        leafIndexSet_ = LeafIndexSet( hostGrid().leafIndexSet() );
      assert( leafIndexSet_ );
      return leafIndexSet_;
    }

    void globalRefine ( int refCount )
    {
      hostGrid().globalRefine( refCount );
      update();
    }

    bool mark ( int refCount, const typename Codim< 0 >::Entity &entity )
    {
      return hostGrid().mark( refCount, getHostEntity< 0 >( entity ) );
    }

    int getMark ( const typename Codim< 0 >::Entity &entity ) const
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
      hostGrid().postAdapt();
    }

    /** \name Parallel Data Distribution and Communication Methods
     *  \{ */

    /** \brief obtain size of overlap region for the leaf grid
     *
     *  \param[in]  codim  codimension for with the information is desired
     */
    int overlapSize ( int codim ) const
    {
      return leafGridView().overlapSize( codim );
    }

    /** \brief obtain size of ghost region for the leaf grid
     *
     *  \param[in]  codim  codimension for with the information is desired
     */
    int ghostSize( int codim ) const
    {
      return leafGridView().ghostSize( codim );
    }

    /** \brief obtain size of overlap region for a grid level
     *
     *  \param[in]  level  grid level (0, ..., maxLevel())
     *  \param[in]  codim  codimension (0, ..., dimension)
     */
    int overlapSize ( int level, int codim ) const
    {
      return levelGridView( level ).overlapSize( codim );
    }

    /** \brief obtain size of ghost region for a grid level
     *
     *  \param[in]  level  grid level (0, ..., maxLevel())
     *  \param[in]  codim  codimension (0, ..., dimension)
     */
    int ghostSize ( int level, int codim ) const
    {
      return levelGridView( level ).ghostSize( codim );
    }

    /** \brief communicate information on a grid level
     *
     *  \param      dataHandle  communication data handle (user defined)
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
    void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                       InterfaceType interface,
                       CommunicationDirection direction,
                       int level ) const
    {
      levelGridView( level ).communicate( dataHandle, interface, direction );
    }

    /** \brief communicate information on leaf entities
     *
     *  \param      dataHandle  communication data handle (user defined)
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
    void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                       InterfaceType interface,
                       CommunicationDirection direction ) const
    {
      leafGridView().communicate( dataHandle, interface, direction );
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

#if 0
    // data handle interface different between geo and interface

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
#endif

    /** \brief obtain EntityPointer from EntitySeed. */
    template< class EntitySeed >
    typename Traits::template Codim< EntitySeed::codimension >::EntityPointer
    DUNE_DEPRECATED_MSG("entityPointer() is deprecated and will be removed after the release of dune-grid 2.4. Use entity() instead to directly obtain an Entity object.")
    entityPointer ( const EntitySeed &seed ) const
    {
      typedef typename Traits::template Codim< EntitySeed::codimension >::Entity Entity;
      return DefaultEntityPointer< Entity >( entity( seed ) );
    }

    /** \brief obtain Entity from EntitySeed. */
    template< class EntitySeed >
    typename Traits::template Codim< EntitySeed::codimension >::Entity
    entity ( const EntitySeed &seed ) const
    {
      typedef typename Traits::template Codim< EntitySeed::codimension >::EntityImpl EntityImpl;
      return EntityImpl( *this, seed );
    }

    /** \} */

    /** \name Grid Views
     *  \{ */

    /** \brief View for a grid level */
    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LevelGridView levelGridView ( int level ) const
    {
      typedef typename Partition< pitype >::LevelGridView View;
      typedef typename View::GridViewImp ViewImp;
      return View( ViewImp( *this, hostGrid().levelGridView( level ) ) );
    }

    /** \brief View for the leaf grid */
    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LeafGridView leafGridView () const
    {
      typedef typename Traits::template Partition< pitype >::LeafGridView View;
      typedef typename View::GridViewImp ViewImp;
      return View( ViewImp( *this, hostGrid().leafGridView() ) );
    }

    /** \brief View for a grid level for All_Partition */
    LevelGridView levelGridView ( int level ) const
    {
      typedef typename LevelGridView::GridViewImp ViewImp;
      return LevelGridView( ViewImp( *this, hostGrid().levelGridView( level ) ) );
    }

    /** \brief View for the leaf grid for All_Partition*/
    LeafGridView leafGridView () const
    {
      typedef typename LeafGridView::GridViewImp ViewImp;
      return LeafGridView( ViewImp( *this, hostGrid().leafGridView() ) );
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
      GeoGrid::AdaptCoordFunction< typename CoordFunction::Interface >::adapt( coordFunction_ );

      const int newNumLevels = maxLevel()+1;
      const int oldNumLevels = levelIndexSets_.size();

      for( int i = newNumLevels; i < oldNumLevels; ++i )
      {
        if( levelIndexSets_[ i ] )
          delete levelIndexSets_[ i ];
      }
      levelIndexSets_.resize( newNumLevels, nullptr );
    }

    /** \} */

    using Base::getRealImplementation;

    const CoordFunction &coordFunction () const { return coordFunction_; }
    CoordFunction &coordFunction () { return coordFunction_; }

  protected:
    template< int codim >
    static const typename HostGrid::template Codim< codim >::Entity &
    getHostEntity( const typename Codim< codim >::Entity &entity )
    {
      return getRealImplementation( entity ).hostEntity();
    }

    void *allocateStorage ( std::size_t size ) const
    {
      return storageAllocator_.allocate( size );
    }

    void deallocateStorage ( void *p, std::size_t size ) const
    {
      storageAllocator_.deallocate( (char *)p, size );
    }

  private:
    HostGrid *const hostGrid_;
    CoordFunction &coordFunction_;
    bool removeHostGrid_;
    mutable std::vector< LevelIndexSet *, typename Allocator::template rebind< LevelIndexSet * >::other > levelIndexSets_;
    mutable LeafIndexSet leafIndexSet_;
    mutable GlobalIdSet globalIdSet_;
    mutable LocalIdSet localIdSet_;
    mutable typename Allocator::template rebind< char >::other storageAllocator_;
  };



  // GeometryGrid::Codim
  // -------------------

  template< class HostGrid, class CoordFunction, class Allocator >
  template< int codim >
  struct GeometryGrid< HostGrid, CoordFunction, Allocator >::Codim
    : public Base::template Codim< codim >
  {
    /** \name Entity and Entity Pointer Types
     *  \{ */

    /** \brief type of entity
     *
     *  The entity is a model of Dune::Entity.
     */
    typedef typename Traits::template Codim< codim >::Entity Entity;

    /** \brief type of entity pointer
     *
     *  The entity pointer is a model of Dune::EntityPointer.
     */
    typedef typename Traits::template Codim< codim >::EntityPointer EntityPointer;

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
    typedef typename Traits::template Codim< codim >::Geometry Geometry;

    /** \brief type of local geometry
     *
     *  Models the geomtry mapping into the reference element of dimension
     *  \em dimension.
     *
     *  The local geometry is a model of Dune::Geometry, implemented through
     *  the generic geometries provided by dune-grid.
     */
    typedef typename Traits::template Codim< codim >::LocalGeometry LocalGeometry;

    /** \} */

    /** \name Iterator Types
     *  \{ */

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename Traits::template Codim< codim >
      ::template Partition< pitype >::LeafIterator
      LeafIterator;
      typedef typename Traits::template Codim< codim >
      ::template Partition< pitype >::LevelIterator
      LevelIterator;
    };

    /** \brief type of level iterator
     *
     *  This iterator enumerates the entites of codimension \em codim of a
     *  grid level.
     *
     *  The level iterator is a model of Dune::LevelIterator.
     */
    typedef typename Partition< All_Partition >::LeafIterator LeafIterator;

    /** \brief type of leaf iterator
     *
     *  This iterator enumerates the entites of codimension \em codim of the
     *  leaf grid.
     *
     *  The leaf iterator is a model of Dune::LeafIterator.
     */
    typedef typename Partition< All_Partition >::LevelIterator LevelIterator;

    /** \} */
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_GRID_HH
