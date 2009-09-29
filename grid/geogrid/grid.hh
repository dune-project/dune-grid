// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_GRID_HH
#define DUNE_GEOGRID_GRID_HH

#include <string>
#include <map>

#include <dune/common/collectivecommunication.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/common/bitfield.hh>
#include <dune/grid/common/grid.hh>
#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>

struct Helix {
  enum {dimW = 3};
  enum {dimG = 2};
  typedef Dune::FieldVector<double,dimW> WCoord;
  typedef Dune::FieldVector<double,dimG> GCoord;
  Helix(const GCoord& x,WCoord& ret) {
    ret[0] = (x[0]+0.2)*cos(2.*M_PI*x[1]);
    ret[1] = (x[0]+0.2)*sin(2.*M_PI*x[1]);
    ret[2] = x[1];
  }
};
struct Circle {
  enum {dimW = 2};
  enum {dimG = 2};
  typedef Dune::FieldVector<double,dimW> WCoord;
  typedef Dune::FieldVector<double,dimG> GCoord;
  Circle(const GCoord& x,WCoord& ret) {
    ret[0] = (x[0]+0.2)*cos(x[1]);
    ret[1] = (x[0]+0.2)*sin(x[1]);
  }
};
typedef Circle func;

#include <dune/grid/geogrid/entity.hh>
#include <dune/grid/geogrid/entitypointer.hh>
#include <dune/grid/geogrid/intersectioniterator.hh>
#include <dune/grid/geogrid/leveliterator.hh>
#include <dune/grid/geogrid/leafiterator.hh>
#include <dune/grid/geogrid/hierarchiciterator.hh>
#include <dune/grid/geogrid/indexsets.hh>

#include <dune/grid/genericgeometry/geometry.hh>

namespace Dune
{

  // Forward declaration
  template< class HostGrid >
  class GeometryGrid;


  namespace GenericGeometry
  {
    template< class HostGrid >
    struct GeometryTraits< GeometryGrid< HostGrid > >
      : public DefaultGeometryTraits
        < typename HostGrid :: ctype, HostGrid :: dimension, func :: dimW >
    {};
  }

  template< class HostGrid >
  struct GeometryGridFamily
  {
    struct Traits
    {
      typedef GeometryGrid< HostGrid > Grid;

      dune_static_assert( (int)HostGrid :: dimensionworld == (int)func :: dimG,
                          "HostGrid and GridFunction are incompatible." );
      enum { dimension = HostGrid :: dimension };
      enum { dimensionworld = func :: dimW };

      typedef Intersection< const Grid, GeometryGridLeafIntersectionIterator >
      LeafIntersection;
      typedef Intersection< const Grid, GeometryGridLevelIntersectionIterator >
      LevelIntersection;

      typedef IntersectionIterator
      < const Grid, GeometryGridLeafIntersectionIterator,
          GeometryGridLeafIntersectionIterator >
      LeafIntersectionIterator;
      typedef IntersectionIterator
      < const Grid, GeometryGridLevelIntersectionIterator,
          GeometryGridLevelIntersectionIterator >
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

      typedef GeometryGridGlobalIdSet< const Grid > GlobalIdSet;
      typedef GeometryGridLocalIdSet< const Grid > LocalIdSet;

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



  //**********************************************************************
  //
  // --GeometryGrid
  //
  //**********************************************************************

  /** \brief [<em> provides \ref Dune::Grid </em>]
   *
   */
  template< class HostGrid >
  class GeometryGrid
    : public GridDefaultImplementation
      < HostGrid :: dimension, func :: dimW, double,
          GeometryGridFamily< HostGrid > >
  {
    friend class GeometryGridLevelIndexSet<const GeometryGrid<HostGrid> >;
    friend class GeometryGridLeafIndexSet<const GeometryGrid<HostGrid> >;
    friend class GeometryGridGlobalIdSet<const GeometryGrid<HostGrid> >;
    friend class GeometryGridLocalIdSet<const GeometryGrid<HostGrid> >;
    friend class GeometryGridHierarchicIterator<const GeometryGrid<HostGrid> >;
    friend class GeometryGridLevelIntersectionIterator<const GeometryGrid<HostGrid> >;
    friend class GeometryGridLeafIntersectionIterator<const GeometryGrid<HostGrid> >;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class GeometryGridLevelIterator;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class GeometryGridLeafIterator;


    template<int codim_, int dim_, class GridImp_>
    friend class GeometryGridEntity;

  public:
    /** \todo Should not be public */
    typedef HostGrid HostGridType;

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    //! type of the used GridFamily for this grid
    typedef GeometryGridFamily< HostGrid > GridFamily;

    //! the Traits
    typedef typename GridFamily :: Traits Traits;

    //! The type used to store coordinates, inherited from the HostGrid
    typedef typename HostGrid :: ctype ctype;

  private:
    HostGrid* hostgrid_;
    std :: vector
    < GeometryGridLevelIndexSet< const GeometryGrid< HostGrid > >* >
    levelIndexSets_;
    GeometryGridLeafIndexSet< const GeometryGrid< HostGrid > > leafIndexSet_;
    GeometryGridGlobalIdSet< const GeometryGrid< HostGrid > > globalIdSet_;
    GeometryGridLocalIdSet< const GeometryGrid< HostGrid > > localIdSet_;

  public:
    /** \brief Constructor
     */
    explicit GeometryGrid(HostGrid& hostgrid) :
      hostgrid_(&hostgrid),
      leafIndexSet_(*this),
      globalIdSet_(*this),
      localIdSet_(*this)
    {
      setIndices();
    }


    //! Desctructor
    ~GeometryGrid()
    {
      // Delete level index sets
      for (size_t i=0; i<levelIndexSets_.size(); i++)
        if (levelIndexSets_[i])
          delete (levelIndexSets_[i]);
    }


    //! return grid name
    std :: string name () const
    {
      return std :: string( "GeometryGrid< " )
             + hostGrid().name() + std :: string( " >" );
    }


    //! Return maximum level defined in this grid. Levels are numbered
    //! 0 ... maxlevel with 0 the coarsest level.
    int maxLevel() const {
      return hostgrid_->maxLevel();
    }


    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
      return GeometryGridLevelIterator<codim,All_Partition, const GeometryGrid<HostGrid> >(this, level);
    }


    //! one past the end on this level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
      return GeometryGridLevelIterator<codim,All_Partition, const GeometryGrid<HostGrid> >(this, level, true);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
      return GeometryGridLevelIterator<codim,PiType, const GeometryGrid<HostGrid> >(this, level);
    }


    //! one past the end on this level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
      return GeometryGridLevelIterator<codim,PiType, const GeometryGrid<HostGrid> >(this, level, true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
      return GeometryGridLeafIterator<codim,All_Partition, const GeometryGrid<HostGrid> >(this);
    }


    //! one past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend() const {
      return GeometryGridLeafIterator<codim,All_Partition, const GeometryGrid<HostGrid> >(this, true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
      return GeometryGridLeafIterator<codim,PiType, const GeometryGrid<HostGrid> >(this);
    }


    //! one past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
      return GeometryGridLeafIterator<codim,PiType, const GeometryGrid<HostGrid> >(this, true);
    }


    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const {
      return hostgrid_->size(level,codim);
    }


    //! number of leaf entities per codim in this process
    int size (int codim) const {
      return leafIndexSet().size(codim);
    }


    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const {
      return levelIndexSets_[level]->size(type);
    }


    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const
    {
      return leafIndexSet().size(type);
    }


    /** \brief Access to the GlobalIdSet */
    const typename Traits::GlobalIdSet& globalIdSet() const {
      return globalIdSet_;
    }


    /** \brief Access to the LocalIdSet */
    const typename Traits::LocalIdSet& localIdSet() const {
      return localIdSet_;
    }


    /** \brief Access to the LevelIndexSets */
    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
      return *levelIndexSets_[level];
    }


    /** \brief Access to the LeafIndexSet */
    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return leafIndexSet_;
    }


    /** @name Grid Refinement Methods */
    /*@{*/


    /** global refinement
     * \todo optimize implementation
     */
    void globalRefine (int refCount)
    {
      hostgrid_->globalRefine(refCount);
    }

    /** \brief Mark entity for refinement
     *
     * This only works for entities of codim 0.
     * The parameter is currently ignored
     *
     * \return <ul>
     * <li> true, if marking was succesfull </li>
     * <li> false, if marking was not possible </li>
     * </ul>
     */
    bool mark(int refCount, const typename Traits::template Codim<0>::EntityPointer & e)
    {
      return hostgrid_->mark(refCount, getHostEntity<0>(*e));
    }

    /** \brief Return refinement mark for entity
     *
     * \return refinement mark (1,0,-1)
     */
    int getMark(const typename Traits::template Codim<0>::EntityPointer & e) const
    {
      return hostgrid_->getMark(getHostEntity<0>(*e));
    }

    //! \todo Please doc me !
    bool preAdapt() {
      return hostgrid_->preAdapt();
    }


    //! Triggers the grid refinement process
    bool adapt()
    {
      return hostgrid_->adapt();
    }

    /** \brief Clean up refinement markers */
    void postAdapt() {
      return hostgrid_->postAdapt();
    }

    /*@}*/

    /** \brief Size of the overlap on the leaf level */
    unsigned int overlapSize(int codim) const {
      return hostgrid_->overlapSize(codim);
    }


    /** \brief Size of the ghost cell layer on the leaf level */
    unsigned int ghostSize(int codim) const {
      return hostgrid_->ghostSize(codim);
    }


    /** \brief Size of the overlap on a given level */
    unsigned int overlapSize(int level, int codim) const {
      return hostgrid_->overlapSize(level,codim);
    }


    /** \brief Size of the ghost cell layer on a given level */
    unsigned int ghostSize(int level, int codim) const {
      return hostgrid_->ghostSize(level,codim);
    }


#if 0
    /** \brief Distributes this grid over the available nodes in a distributed machine
     *
     * \param minlevel The coarsest grid level that gets distributed
     * \param maxlevel does currently get ignored
     */
    void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
      DUNE_THROW(NotImplemented, "GeometryGrid::loadBalance()");
    }

    /** \brief The communication interface
     *  @param T: array class holding data associated with the entities
     *  @param P: type used to gather/scatter data in and out of the message buffer
     *  @param codim: communicate entites of given codim
     *  @param if: one of the predifined interface types, throws error if it is not implemented
     *  @param level: communicate for entities on the given level
     *
     *  Implements a generic communication function sending an object of type P for each entity
     *  in the intersection of two processors. P has two methods gather and scatter that implement
     *  the protocol. Therefore P is called the "protocol class".
     */
    template<class T, template<class> class P, int codim>
    void communicate (T& t, InterfaceType iftype, CommunicationDirection dir, int level);

    /*! The new communication interface

       communicate objects for all codims on a given level
     */
    template<class DataHandle>
    void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {}

    template<class DataHandle>
    void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
    {}
#endif


    //! Obtain CollectiveCommunication object (taken from host grid)
    const CollectiveCommunication<GeometryGrid>& comm () const
    {
      return hostGrid().comm();
    }


    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    //! Returns the hostgrid this GeometryGrid lives in
    HostGridType& getHostGrid() const
    {
      return *hostgrid_;
    }

    //! Obtain the host grid wrapped this this GeometryGrid
    const HostGridType &hostGrid () const
    {
      return *hostgrid_;
    }

    //! Returns the hostgrid entity encapsulated in given subgrid entity
    template <int codim>
    typename HostGrid::Traits::template Codim<codim>::EntityPointer getHostEntity(const typename Traits::template Codim<codim>::Entity& e) const
    {
      return getRealImplementation(e).hostEntity_;
    }

  private:
    //! compute the grid indices and ids
    void setIndices ()
    {
      localIdSet_.update();
      globalIdSet_.update();

      // //////////////////////////////////////////
      //   Create the index sets
      // //////////////////////////////////////////
      for( int i = levelIndexSets_.size(); i <= maxLevel(); ++i )
      {
        GeometryGridLevelIndexSet< const GeometryGrid< HostGrid > > *p
          = new GeometryGridLevelIndexSet< const GeometryGrid< HostGrid > >();
        levelIndexSets_.push_back( p );
      }

      for( int i = 0; i <= maxLevel(); ++i )
      {
        if( levelIndexSets_[ i ] )
          levelIndexSets_[ i ]->update( *this, i );
      }
      leafIndexSet_.update(*this);
    }
  }; // end Class GeometryGrid



  // Capabilities
  // ------------

  namespace Capabilities
  {
    //! \todo Please doc me !
    template<class HostGrid, int codim>
    struct hasEntity< GeometryGrid<HostGrid>, codim>
    {
      static const bool v = hasEntity<HostGrid,codim>::v;
    };


    //! \todo Please doc me !
    template<class HostGrid>
    struct isParallel< GeometryGrid<HostGrid> >
    {
      static const bool v = isParallel<HostGrid>::v;
    };


    //! \todo Please doc me !
    template<class HostGrid>
    struct hasHangingNodes< GeometryGrid<HostGrid> >
    {
      static const bool v = hasHangingNodes<HostGrid>::v;
    };

    //! \todo Please doc me !
    template<class HostGrid>
    struct isLevelwiseConforming< GeometryGrid<HostGrid> >
    {
      static const bool v = isLevelwiseConforming<HostGrid>::v;
    };

    //! \todo Please doc me !
    template<class HostGrid>
    struct isLeafwiseConforming< GeometryGrid<HostGrid> >
    {
      static const bool v = isLeafwiseConforming<HostGrid>::v;
    };
  }

} // namespace Dune

#endif
