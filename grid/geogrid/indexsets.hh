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


  // GeometryGridLevelIndexSetTypes
  // ------------------------------

  template <class GridImp>
  struct GeometryGridLevelIndexSetTypes
  {
    //! The types
    template<int cd>
    struct Codim
    {
      template<PartitionIteratorType pitype>
      struct Partition
      {
        typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LevelIterator Iterator;
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




  /** \todo Take the index types from the host grid */
  template<class GridImp>
  class GeometryGridLevelIndexSet :
    public IndexSetDefaultImplementation<GridImp,GeometryGridLevelIndexSet<GridImp>,GeometryGridLevelIndexSetTypes<GridImp> >
  {
  public:

    typedef typename remove_const<GridImp>::type::HostGridType HostGrid;

    enum {dim = GridImp::dimension};

    typedef IndexSet<GridImp,GeometryGridLevelIndexSet<GridImp>,GeometryGridLevelIndexSetTypes<GridImp> > Base;

    //! get index of an entity
    template<int codim>
    int index (const typename GridImp::Traits::template Codim<codim>::Entity& e) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).template index<codim>(*grid_->template getHostEntity<codim>(e));
    }


    //! get index of subEntity of a codim 0 entity
    template<int codim>
    int subIndex (const typename GridImp::Traits::template Codim<0>::Entity& e, int i) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).template subIndex<codim>(*grid_->template getHostEntity<0>(e), i);
    }


    //! get number of entities of given codim, type and on this level
    int size (int codim) const {
      return grid_->hostgrid_->levelIndexSet(level_).size(codim);
    }


    //! get number of entities of given codim, type and on this level
    int size (GeometryType type) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).size(type);
    }


    /** \brief Deliver all geometry types used in this grid */
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).geomTypes(codim);
    }


    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Base::template Codim<cd>::template Partition<pitype>::Iterator begin () const
    {
      return grid_->template lbegin<cd,pitype>(level_);
    }


    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    typename Base::template Codim<cd>::template Partition<pitype>::Iterator end () const
    {
      return grid_->template lend<cd,pitype>(level_);
    }


    /** \brief Set up the index set */
    void update(const GridImp& grid, int level)
    {
      grid_ = &grid;
      level_ = level;
    }


    GridImp* grid_;

    int level_;
  };




  template <class GridImp>
  struct GeometryGridLeafIndexSetTypes
  {
    //! The types
    template<int cd>
    struct Codim
    {
      template<PartitionIteratorType pitype>
      struct Partition
      {
        typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LeafIterator Iterator;
      };
    };
  };




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




  template <class GridImp>
  class GeometryGridGlobalIdSet :
    public IdSet<GridImp,GeometryGridGlobalIdSet<GridImp>,
        typename remove_const<GridImp>::type::HostGridType::Traits::GlobalIdSet::IdType>
  {

    typedef typename remove_const<GridImp>::type::HostGridType HostGrid;


  public:
    //! constructor stores reference to a grid
    GeometryGridGlobalIdSet (const GridImp& g) : grid_(&g) {}

    //! define the type used for persistent indices
    typedef typename HostGrid::Traits::GlobalIdSet::IdType GlobalIdType;


    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instatiated yet.
     */
    template<int cd>
    GlobalIdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      // Return id of the host entity
      return grid_->hostgrid_->globalIdSet().id(*grid_->getRealImplementation(e).hostEntity_);
    }


    //! get id of subEntity
    /*
        We use the remove_const to extract the Type from the mutable class,
        because the const class is not instatiated yet.
     */
    template<int cc>
    GlobalIdType subId (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i) const
    {
      // Return sub id of the host entity
      return grid_->hostgrid_->globalIdSet().template subId<cc>(*grid_->getRealImplementation(e).hostEntity_,i);
    }


    /** \todo Should be private */
    void update() {}


    const GridImp* grid_;
  };




  template<class GridImp>
  class GeometryGridLocalIdSet :
    public IdSet<GridImp,GeometryGridLocalIdSet<GridImp>,
        typename remove_const<GridImp>::type::HostGridType::Traits::LocalIdSet::IdType>
  {
  private:

    typedef typename remove_const<GridImp>::type::HostGridType HostGrid;


  public:
    //! define the type used for persistent local ids
    typedef typename HostGrid::Traits::LocalIdSet::IdType LocalIdType;


    //! constructor stores reference to a grid
    GeometryGridLocalIdSet (const GridImp& g) : grid_(&g) {}


    //! get id of an entity
    /*
        We use the remove_const to extract the Type from the mutable class,
        because the const class is not instatiated yet.
     */
    template<int cd>
    LocalIdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      // Return id of the host entity
      return grid_->hostgrid_->localIdSet().id(*grid_->getRealImplementation(e).hostEntity_);
    }


    //! get id of subEntity
    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instatiated yet.
     */
    template<int cc>
    LocalIdType subId (const typename remove_const<GridImp>::type::template Codim<0>::Entity& e, int i) const
    {
      // Return sub id of the host entity
      return grid_->hostgrid_->localIdSet().template subId<cc>(*grid_->getRealImplementation(e).hostEntity_,i);
    }


    /** \todo Should be private */
    void update() {}


    const GridImp* grid_;
  };


}  // namespace Dune


#endif
