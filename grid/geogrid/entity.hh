// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITY_HH
#define DUNE_GEOGRID_ENTITY_HH

#include <dune/grid/common/referenceelements.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;

  template<int codim, class GridImp>
  class GeometryGridEntityPointer;

  template<int codim, PartitionIteratorType pitype, class GridImp>
  class GeometryGridLevelIterator;

  template<class GridImp>
  class GeometryGridLevelIntersectionIterator;

  template<class GridImp>
  class GeometryGridLeafIntersectionIterator;

  template<class GridImp>
  class GeometryGridHierarchicIterator;



  // Internal Forward Declarations
  // -----------------------------

  template< int codim, int dim, class Grid >
  class GeometryGridEntity;




  // GeometryGridMakeableEntity
  // --------------------------

  template<int codim, int dim, class GridImp>
  class GeometryGridMakeableEntity :
    public GridImp::template Codim<codim>::Entity
  {
  public:

    // The codimension of this entitypointer wrt the host grid
    enum {CodimInHostGrid = GridImp::HostGridType::dimension - GridImp::dimension + codim};

    // EntityPointer to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::EntityPointer HostGridEntityPointer;


    //! \todo Please doc me !
    GeometryGridMakeableEntity(const GridImp* identityGrid, const HostGridEntityPointer& hostEntity) :
      GridImp::template Codim<codim>::Entity (GeometryGridEntity<codim, dim, const GridImp>(identityGrid,hostEntity)),
      identityGrid_(identityGrid)
    {}


    //! \todo Please doc me !
    void setToTarget(const HostGridEntityPointer& hostEntity) {
      this->realEntity.setToTarget(hostEntity);
    }


    //! \todo Please doc me !
    const HostGridEntityPointer& getTarget() {
      return this->realEntity.hostEntity_;
    }


  private:

    const GridImp* identityGrid_;
  };



  /** \brief The implementation of entities in a GeometryGrid
   *   \ingroup GeometryGrid
   *
   *  A Grid is a container of grid entities. An entity is parametrized by the codimension.
   *  An entity of codimension c in dimension d is a d-c dimensional object.
   */
  template< int codim, int dim, class HostGrid, class CoordFunction >
  class GeometryGridEntity< codim, dim, const GeometryGrid< HostGrid, CoordFunction > >
    : public EntityDefaultImplementation
      < codim, dim, const GeometryGrid< HostGrid, CoordFunction >, GeometryGridEntity >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    friend class GeometryGridMakeableEntity< codim, dim, const Grid >;

    template< class > friend class GeometryGridLevelIndexSet;
    template< class > friend class GeometryGridLeafIndexSet;
    template< class > friend class GeometryGridLocalIdSet;
    template< class > friend class GeometryGridGlobalIdSet;
    template< class, int > friend class IndexSetter;

    typedef typename Grid :: ctype ctype;

    typedef typename HostGrid :: Traits :: template Codim< codim > :: Entity
    HostEntity;
    typedef typename HostGrid :: Traits :: template Codim< codim > :: EntityPointer
    HostEntityPointer;
    typedef typename HostGrid :: Traits :: template Codim< codim > :: Geometry
    HostGeometry;

  public:
    typedef typename Grid :: template Codim< codim > :: Geometry Geometry;

  private:
    typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;
    typedef typename GeometryImpl :: GlobalCoordinate GlobalCoordinate;

  public:
    HostEntityPointer hostEntity_;

  private:
    const Grid *grid_;
    mutable std :: vector< GlobalCoordinate > corners_;
    mutable Geometry *geo_;

  public:
    //! Constructor for an entity in a given grid level
    GeometryGridEntity( const Grid *grid, const HostEntityPointer &hostEntity )
      : hostEntity_( hostEntity ),
        grid_( grid ),
        geo_( 0 )
    {}

    //! \todo Please doc me !
    GeometryGridEntity( const GeometryGridEntity &other )
      : hostEntity_( other.hostEntity_ ),
        grid_( other.grid_ ),
        geo_( 0 )
    {}

    //! Destructor
    ~GeometryGridEntity ()
    {
      if( geo_ != 0 )
        delete geo_;
    }

    GeometryGridEntity &operator= ( const GeometryGridEntity &other )
    {
      if( this == &other )
        return *this;

      if( geo_ != 0 )
      {
        delete geo_;
        geo_ = 0;
      }

      grid_ = other.grid_;
      hostEntity_ = other.hostEntity_;
      return *this;
    }

    GeometryType type () const
    {
      return hostEntity().type();
    }

    //! level of this element
    int level () const
    {
      return hostEntity().level();
    }

    //! The partition type for parallel computing
    PartitionType partitionType () const
    {
      return hostEntity().partitionType();
    }

#if 0
    /** Intra-element access to entities of codimension cc > codim. Return number of entities
     * with codimension cc.
     */
    template<int cc> int count () const {
      return hostEntity_->template count<cc>();
    }
#endif

    //! Geometry of this entity
    const Geometry &geometry () const
    {
      if( geo_ ==0 )
      {
        const HostGeometry &hostGeo = hostEntity().geometry();
        corners_.resize( hostGeo.corners() );
        for( unsigned int i = 0; i < corners_.size(); ++i )
          coordFunction().evaluate( hostGeo[ i ], corners_[ i ] );
        geo_ = new MakeableGeometry( GeometryImpl( type(), corners_ ) );
      }
      return *geo_;
    }

    const HostEntity &hostEntity () const
    {
      return *hostEntity_;
    }

  private:
    const CoordFunction &coordFunction () const
    {
      return grid_->coordFunction();
    }

    void setToTarget( const HostEntityPointer &target )
    {
      if( geo_ != 0 )
      {
        delete geo_;
        geo_ = 0;
      }
      hostEntity_ = target;
    }
  };




  /** \brief Specialization for codim-0-entities.
   * \ingroup GeometryGrid
   *
   * This class embodies the topological parts of elements of the grid.
   * It has an extended interface compared to the general entity class.
   * For example, Entities of codimension 0  allow to visit all neighbors.
   */
  template< int dim, class HostGrid, class CoordFunction >
  class GeometryGridEntity< 0, dim, const GeometryGrid< HostGrid, CoordFunction > >
    : public EntityDefaultImplementation
      < 0, dim, const GeometryGrid< HostGrid, CoordFunction >, GeometryGridEntity >
  {
    enum { codim = 0 };

    typedef GeometryGrid< HostGrid, CoordFunction > Grid;

    friend class GeometryGridMakeableEntity< codim, dim, const Grid >;

    template< class > friend class GeometryGridLevelIndexSet;
    template< class > friend class GeometryGridLeafIndexSet;
    template< class > friend class GeometryGridLocalIdSet;
    template< class > friend class GeometryGridGlobalIdSet;
    template< class, int > friend class IndexSetter;

    typedef typename Grid :: ctype ctype;

    typedef typename HostGrid :: Traits :: template Codim< codim > :: Entity
    HostEntity;
    typedef typename HostGrid :: Traits :: template Codim< codim > :: EntityPointer
    HostEntityPointer;
    typedef typename HostGrid :: Traits :: template Codim< codim > :: Geometry
    HostGeometry;

  public:
    typedef typename Grid :: template Codim< codim > :: EntityPointer EntityPointer;
    typedef typename Grid :: template Codim< codim > :: Geometry Geometry;
    typedef typename Grid :: template Codim< codim > :: LocalGeometry LocalGeometry;

    typedef typename Grid :: HierarchicIterator HierarchicIterator;
    typedef typename Grid :: LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename Grid :: LevelIntersectionIterator LevelIntersectionIterator;

  private:
    typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
    typedef typename MakeableGeometry :: ImplementationType GeometryImpl;
    typedef typename GeometryImpl :: GlobalCoordinate GlobalCoordinate;

  public:
    HostEntityPointer hostEntity_;

  private:
    const Grid *grid_;
    mutable std :: vector< GlobalCoordinate > corners_;
    mutable Geometry *geo_;
    mutable LocalGeometry *geoInFather_;

  public:
    GeometryGridEntity ( const Grid *grid, const HostEntityPointer &hostEntity )
      : hostEntity_( hostEntity ),
        grid_( grid ),
        geo_( 0 ),
        geoInFather_( 0 )
    {}

    GeometryGridEntity( const GeometryGridEntity &other )
      : hostEntity_( other.hostEntity_ ),
        grid_( other.grid_ ),
        geo_( 0 ),
        geoInFather_( 0 )
    {}

    ~GeometryGridEntity ()
    {
      if( geo_ != 0 )
        delete geo_;
      if( geoInFather_ != 0 )
        delete geoInFather_;
    }


    GeometryGridEntity &operator= ( const GeometryGridEntity &other )
    {
      if( this == &other )
        return *this;

      if( geo_ != 0 )
      {
        delete geo_;
        geo_ = 0;
      }

      if( geoInFather_ != 0 )
      {
        delete geoInFather_;
        geoInFather_ = 0;
      }

      grid_ = other.grid_;
      hostEntity_ = other.hostEntity_;
      return *this;
    }

    GeometryType type () const
    {
      return hostEntity().type();
    }

    //! level of this element
    int level () const
    {
      return hostEntity().level();
    }

    //! The partition type for parallel computing
    PartitionType partitionType () const
    {
      return hostEntity().partitionType();
    }

    //! Geometry of this entity
    const Geometry &geometry () const
    {
      if( geo_ ==0 )
      {
        const HostGeometry &hostGeo = hostEntity().geometry();
        corners_.resize( hostGeo.corners() );
        for( unsigned int i = 0; i < corners_.size(); ++i )
          coordFunction().evaluate( hostGeo[ i ], corners_[ i ] );
        geo_ = new MakeableGeometry( GeometryImpl( type(), corners_ ) );
      }
      return *geo_;
    }

    template< int cc >
    int count () const
    {
      return hostEntity().template count< cc >();
    }

    template< int cc >
    typename Grid :: template Codim< cc > :: EntityPointer
    entity ( int i ) const
    {
      typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
      EntityPointerImpl impl( grid_, hostEntity().template entity< cc >( i ) );
      return MakeableEntityPointer( impl );
    }

    LevelIntersectionIterator ilevelbegin () const
    {
      typedef MakeableInterfaceObject< LevelIntersectionIterator >
      MakeableLevelIntersectionIterator;
      typedef typename MakeableLevelIntersectionIterator :: ImplementationType
      LevelIntersectionIteratorImpl;
      LevelIntersectionIteratorImpl impl( grid_, hostEntity().ilevelbegin() );
      return MakeableLevelIntersectionIterator( impl );
    }

    LevelIntersectionIterator ilevelend () const
    {
      typedef MakeableInterfaceObject< LevelIntersectionIterator >
      MakeableLevelIntersectionIterator;
      typedef typename MakeableLevelIntersectionIterator :: ImplementationType
      LevelIntersectionIteratorImpl;
      LevelIntersectionIteratorImpl impl( grid_, hostEntity().ilevelend() );
      return MakeableLevelIntersectionIterator( impl );
    }

    LeafIntersectionIterator ileafbegin () const
    {
      typedef MakeableInterfaceObject< LeafIntersectionIterator >
      MakeableLeafIntersectionIterator;
      typedef typename MakeableLeafIntersectionIterator :: ImplementationType
      LeafIntersectionIteratorImpl;
      LeafIntersectionIteratorImpl impl( grid_, hostEntity().ileafbegin() );
      return MakeableLeafIntersectionIterator( impl );
    }

    LeafIntersectionIterator ileafend () const
    {
      typedef MakeableInterfaceObject< LeafIntersectionIterator >
      MakeableLeafIntersectionIterator;
      typedef typename MakeableLeafIntersectionIterator :: ImplementationType
      LeafIntersectionIteratorImpl;
      LeafIntersectionIteratorImpl impl( grid_, hostEntity().ileafend() );
      return MakeableLeafIntersectionIterator( impl );
    }

    bool isLeaf () const
    {
      return hostEntity().isLeaf();
    }

    EntityPointer father () const
    {
      typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
      typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
      return MakeableEntityPointer( EntityPointerImpl( grid_, hostEntity().father() ) );
    }

    const LocalGeometry &geometryInFather () const
    {
      typedef MakeableInterfaceObject< LocalGeometry > MakeableLocalGeometry;
      typedef typename MakeableLocalGeometry :: ImplementationType LocalGeometryImpl;

      if( geoInFather_ == 0 )
        geoInFather_ = new MakeableLocalGeometry( LocalGeometryImpl( hostEntity().geometryInFather() ) );
      return *geoInFather_;
    }

    HierarchicIterator hbegin ( int maxLevel ) const
    {
      typedef MakeableInterfaceObject< HierarchicIterator > MakeableHierarchicIterator;
      typedef typename MakeableHierarchicIterator :: ImplementationType HierarchicIteratorImpl;
      return MakeableHierarchicIterator( HierarchicIteratorImpl( grid_, *this, maxLevel ) );
    }

    HierarchicIterator hend ( int maxLevel ) const
    {
      typedef MakeableInterfaceObject< HierarchicIterator > MakeableHierarchicIterator;
      typedef typename MakeableHierarchicIterator :: ImplementationType HierarchicIteratorImpl;
      return MakeableHierarchicIterator( HierarchicIteratorImpl( grid_, *this, maxLevel, true ) );
    }

    bool wasRefined () const
    {
      if( grid_->adaptationStep != Grid :: adaptDone )
        return false;

      int index = grid_->levelIndexSet( level() ).index( *this );
      return grid_->refinementMark_[ level ][ index ];
    }

    bool mightBeCoarsened () const
    {
      return true;
    }

    const HostEntity &hostEntity () const
    {
      return *hostEntity_;
    }

  private:
    const CoordFunction &coordFunction () const
    {
      return grid_->coordFunction();
    }

    void setToTarget( const HostEntityPointer &target )
    {
      if( geo_ != 0 )
      {
        delete geo_;
        geo_ = 0;
      }
      if( geoInFather_ != 0 )
      {
        delete geoInFather_;
        geoInFather_ = 0;
      }
      hostEntity_ = target;
    }
  }; // end of GeometryGridEntity codim = 0

} // namespace Dune

#endif
