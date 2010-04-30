// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITY_HH
#define DUNE_GEOGRID_ENTITY_HH

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid/geometrygrid/capabilities.hh>
#include <dune/grid/geometrygrid/cornerstorage.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    /** \class EntityBase
     *  \brief actual implementation of the entity
     *  \ingroup GeoGrid
     *
     *  \tparam  codim  codimension of the entity
     *  \tparam  Grid   GeometryGrid, this entity belongs to
     *  \tparam  fake   \b true, if the host grid does not provide this entity
     *                  (do not specify, the defualt value is already the
     *                  intended use)
     */
    template< int codim, class Grid,
        bool fake = !(Capabilities :: hasHostEntity< Grid, codim > :: v) >
    class EntityBase;

    /** \class Entity
     *  \brief DUNE-conform implementation of the entity
     *  \ingroup GeoGrid
     *
     *  This class merely changes the template parameters of the entity to make
     *  DUNE happy. The actual implementation of the entity can be found in
     *  EntityBase.
     *
     *  \tparam  codim  codimension of the entity
     *  \tparam  dim    dimension of the Grid (redundant information)
     *  \tparam  Grid   GeometryGrid, this entity belongs to
     */
    template< int codim, int dim, class Grid >
    class Entity;



    // EntityBase (real)
    // -----------------

    /** \copydoc EntityBase
     *
     *  This specialization implements the case, where the host grid provides
     *  the entity for this codimension, i.e., \em fake = \b false.
     *
     *  \nosubgrouping
     */
    template< int codim, class Grid >
    class EntityBase< codim, Grid, false >
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

    public:
      /** \name Attributes
       *  \{ */

      //! codimensioon of the entity
      static const int codimension = codim;
      //! dimension of the grid
      static const int dimension = Traits::dimension;
      //! dimension of the entity
      static const int mydimension = dimension - codimension;
      //! dimension of the world
      static const int dimensionworld = Traits::dimensionworld;

      //! \b true, if the entity is faked, i.e., if there is no corresponding host entity
      static const bool fake = false;

      /** \} */

      /** \name Types Required by DUNE
       *  \{ */

      //! coordinate type of the grid
      typedef typename Traits::ctype ctype;

      //! type of corresponding geometry
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;
      /** \} */

    private:
      typedef typename Traits::HostGrid HostGrid;
      typedef typename Traits::CoordFunction CoordFunction;

      typedef typename Traits::template Codim< codimension >::EntityNumbering EntityNumbering;

    public:
      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;
      //! type of corresponding host entity pointer
      typedef typename HostGrid::template Codim< codimension >::EntityPointer HostEntityPointer;

      //! type of host elements, i.e., of host entities of codimension 0
      typedef typename HostGrid::template Codim< 0 >::Entity HostElement;
      /** \} */

    private:
      typedef typename HostGrid::template Codim< codimension >::Geometry HostGeometry;

      typedef typename GenericGeometry::GlobalGeometryTraits< Grid >::template Codim<codimension>::CoordVector CoordVector;

      typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
      typedef typename MakeableGeometry::ImplementationType GeometryImpl;

    public:
      /** \name Construction, Initialization and Destruction
       *  \{ */

      /** \brief construct an initialized entity
       *
       *  \param[in]  grid        GeometryGrid this entity belongs to
       *  \param[in]  hostEntity  corresponding entity in the host grid
       *
       *  \note Both references must remain valid as long as this entity is in
       *        use.
       */
      EntityBase ( const Grid &grid, const HostEntity &hostEntity )
        : grid_( &grid ),
          hostEntity_( &hostEntity ),
          numbering_( grid.numbering()[ hostEntity ] ),
          geo_( GeometryImpl( grid.allocator() ) )
      {}

      EntityBase ( const EntityBase &other )
        : grid_( other.grid_ ),
          hostEntity_( other.hostEntity_ ),
          numbering_( other.numbering_ ),
          geo_( GeometryImpl( grid().allocator() ) )
      {}

      /** \} */

    private:
      EntityBase &operator= ( const EntityBase & );

    public:
      /** \name Methods Shared by Entities of All Codimensions
       *  \{ */

      /** \brief obtain the name of the corresponding reference element
       *
       *  This type can be used to access the DUNE reference element.
       */
      GeometryType type () const
      {
        return hostEntity().type();
      }

      unsigned int topologyId () const
      {
        return GenericGeometry::topologyId( type() );
      }

      /** \brief obtain the level of this entity */
      int level () const
      {
        return hostEntity().level();
      }

      /** \brief obtain the partition type of this entity */
      PartitionType partitionType () const
      {
        return hostEntity().partitionType();
      }

      /** obtain the geometry of this entity
       *
       *  Each DUNE entity encapsulates a geometry object, representing the map
       *  from the reference element to world coordinates. Wrapping the geometry
       *  is the main objective of the GeometryGrid.
       *
       *  The GeometryGrid provides geometries of order 1, obtained by
       *  interpolation of its corners \f$y_i\f$. There corners are calculated
       *  from the corners \f$x_i\f$ of the host geometry through the
       *  GeometryGrid's coordinate function \f$c\f$, i.e.,
       *  \f$y_i = c( x_i )\f$.
       *
       *  \returns a const reference to the geometry
       */
      const Geometry &geometry () const
      {
        GeometryImpl &geo = Grid::getRealImplementation( geo_ );
        if( !geo )
        {
          CoordVector coords( hostEntity(), grid().coordFunction() );
          geo = GeometryImpl( topologyId(), coords, grid().allocator() );
        }
        return geo_;
      }
      /** \} */


      /** \name Methods Supporting the Grid Implementation
       *  \{ */

      const Grid &grid () const
      {
        return *grid_;
      }

      const HostEntity &hostEntity () const
      {
        return *hostEntity_;
      }

      /** \brief obtain the entity's index from a host IndexSet
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param[in]  indexSet  host IndexSet to use
       */
      template< class HostIndexSet >
      typename HostIndexSet::IndexType
      index ( const HostIndexSet &indexSet ) const
      {
        return indexSet.template index< codimension >( hostEntity() );
      }

      /** \brief obtain the index of a subentity from a host IndexSet
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param[in]  indexSet  host IndexSet to use
       *  \param[in]  i         number of the subentity
       *  \param[in]  cd        codimension of the subentity
       */
      template< class HostIndexSet >
      typename HostIndexSet::IndexType
      subIndex ( const HostIndexSet &indexSet, int i, unsigned int cd ) const
      {
        return indexSet.subIndex( hostEntity(), i, cd );
      }

      /** \brief check whether the entity is contained in a host index set
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param  indexSet  host IndexSet to use
       */
      template< class HostIndexSet >
      bool isContained ( const HostIndexSet &indexSet ) const
      {
        return indexSet.contains( hostEntity() );
      }

      /** \brief obtain the entity's id from a host IdSet
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param  idSet  host IdSet to use
       */
      template< class HostIdSet >
      typename HostIdSet::IdType id ( const HostIdSet &idSet ) const
      {
        return idSet.template id< codimension >( hostEntity() );
      }
      /** \} */

    private:
      const Grid *grid_;
      const HostEntity *hostEntity_;
      EntityNumbering numbering_;
      mutable MakeableGeometry geo_;
    };



    // EntityBase (fake)
    // -----------------

    /** \copydoc EntityBase
     *
     *  This specialization implements the case, where the host grid does not
     *  provide the entity for this codimension, i.e., \em fake = \b true.
     *
     *  \nosubgrouping
     */
    template< int codim, class Grid >
    class EntityBase< codim, Grid, true >
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

    public:
      /** \name Attributes
       *  \{ */

      //! codimensioon of the entity
      static const int codimension = codim;
      //! dimension of the grid
      static const int dimension = Traits::dimension;
      //! dimension of the entity
      static const int mydimension = dimension - codimension;
      //! dimension of the world
      static const int dimensionworld = Traits::dimensionworld;

      //! \b true, if the entity is faked, i.e., if there is no corresponding host entity
      static const bool fake = true;
      /** \} */

      /** \name Types Required by DUNE
       *  \{ */

      //! coordinate type of the grid
      typedef typename Traits::ctype ctype;

      //! type of corresponding geometry
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;
      /** \} */

    private:
      typedef typename Traits::HostGrid HostGrid;
      typedef typename Traits::CoordFunction CoordFunction;

      typedef typename Traits::template Codim< 0 >::EntityNumbering ElementNumbering;

    public:
      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;
      //! type of corresponding host entity pointer
      typedef typename HostGrid::template Codim< codimension >::EntityPointer HostEntityPointer;

      //! type of host elements, i.e., of host entities of codimension 0
      typedef typename HostGrid::template Codim< 0 >::Entity HostElement;
      /** \} */

    private:
      typedef typename HostGrid::template Codim< 0 >::Geometry HostGeometry;
      typedef typename HostGrid::template Codim< dimension >::EntityPointer HostVertexPointer;

      typedef typename GenericGeometry::GlobalGeometryTraits<Grid> :: template Codim<codimension> :: CoordVector CoordVector;

      typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
      typedef typename MakeableGeometry::ImplementationType GeometryImpl;

    public:
      /** \name Construction, Initialization and Destruction
       *  \{ */

      /** \brief construct an initialized entity
       *
       *  \param[in]  grid        GeometryGrid this entity belongs to
       *  \param[in]  hostElement any host element containing the corresponding
       *                          host entity
       *  \param[in]  subEntity   number of this entity within the host element
       *
       *  \note Both references must remain valid as long as this entity is in
       *        use.
       */
      EntityBase ( const Grid &grid, const HostElement &hostElement, int subEntity )
        : grid_( &grid ),
          hostElement_( &hostElement ),
          numbering_( grid.numbering()[ hostElement ] ),
          subEntity_( subEntity ),
          geo_( GeometryImpl( grid.allocator() ) )
      {}

      EntityBase ( const EntityBase &other )
        : grid_( other.grid_ ),
          hostElement_( other.hostElement_ ),
          numbering_( other.numbering_ ),
          subEntity_( other.subEntity_ ),
          geo_( GeometryImpl( grid().allocator() ) )
      {}

      /** \} */

    private:
      EntityBase &operator= ( const EntityBase & );

    public:
      /** \name Methods Shared by Entities of All Codimensions
       *  \{ */

      /** \brief obtain the name of the corresponding reference element
       *
       *  This type can be used to access the DUNE reference element.
       */
      GeometryType type () const
      {
        const GenericReferenceElement< ctype, dimension > &refElement
          = GenericReferenceElements< ctype, dimension >::general( hostElement().type() );
        return refElement.type( subEntity_, codimension );
      }

      unsigned int topologyId () const
      {
        const GenericReferenceElement< ctype, dimension > &refElement
          = GenericReferenceElements< ctype, dimension >::general( hostElement().type() );
        return refElement.topologyId( subEntity_, codimension );
      }

      /** \brief obtain the level of this entity */
      int level () const
      {
        return hostElement().level();
      }

      /** \brief obtain the partition type of this entity */
      PartitionType partitionType () const
      {
        if( !(Capabilities::isParallel< HostGrid >::v) )
          return InteriorEntity;

        const GenericReferenceElement< ctype, dimension > &refElement
          = GenericReferenceElements< ctype, dimension >::general( hostElement().type() );

        PartitionType type = vertexPartitionType( refElement, 0 );
        if( (type != BorderEntity) && (type != FrontEntity) )
          return type;

        const int numVertices = refElement.size( subEntity_, codimension, dimension );
        for( int i = 1; i < numVertices; ++i )
        {
          PartitionType vtxType = vertexPartitionType( refElement, i );
          if( (vtxType != BorderEntity) && (vtxType != FrontEntity) )
            return vtxType;
          if( type != vtxType )
            return OverlapEntity;
        }
        assert( (type == BorderEntity) || (type == FrontEntity) );
        return type;
      }

      /** obtain the geometry of this entity
       *
       *  Each DUNE entity encapsulates a geometry object, representing the map
       *  from the reference element to world coordinates. Wrapping the geometry
       *  is the main objective of the GeometryGrid.
       *
       *  The GeometryGrid provides geometries of order 1, obtained by
       *  interpolation of its corners \f$y_i\f$. There corners are calculated
       *  from the corners \f$x_i\f$ of the host geometry through the
       *  GeometryGrid's coordinate function \f$c\f$, i.e.,
       *  \f$y_i = c( x_i )\f$.
       *
       *  \returns a const reference to the geometry
       */
      const Geometry &geometry () const
      {
        GeometryImpl &geo = Grid::getRealImplementation( geo_ );
        if( !geo )
        {
          CoordVector coords( hostElement(), subEntity_, grid().coordFunction() );
          geo = GeometryImpl( topologyId(), coords, grid().allocator() );
        }
        return geo_;
      }
      /** \} */

      /** \name Methods Supporting the Grid Implementation
       *  \{ */

      const Grid &grid () const
      {
        return *grid_;
      }

      const HostElement &hostElement () const
      {
        return *hostElement_;
      }

      int subEntity () const
      {
        return subEntity_;
      }

      const HostEntity &hostEntity () const
      {
        DUNE_THROW( NotImplemented, "HostGrid has no entities of codimension "
                    << codimension << "." );
      }

      /** \brief obtain the entity's index from a host IndexSet
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param[in]  indexSet  host IndexSet to use
       */
      template< class HostIndexSet >
      typename HostIndexSet::IndexType index ( const HostIndexSet &indexSet ) const
      {
        return indexSet.subIndex( hostElement(), subEntity_, codimension );
      }

      /** \brief obtain the index of a subentity from a host IndexSet
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param[in]  indexSet  host IndexSet to use
       *  \param[in]  i         number of the subentity
       *  \param[in]  cd        codimension of the subentity
       */
      template< class HostIndexSet >
      typename HostIndexSet::IndexType
      subIndex ( const HostIndexSet &indexSet, int i, unsigned int cd ) const
      {
        const GenericReferenceElement< ctype, dimension > &refElement
          = GenericReferenceElements< ctype, dimension >::general( hostElement().type() );
        const int j = refElement.subEntity( subEntity_, codimension, i, codimension+cd );
        return indexSet.subIndex( hostElement(), j, codimension+cd );
      }

      /** \brief check whether the entity is contained in a host index set
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param  indexSet  host IndexSet to use
       */
      template< class HostIndexSet >
      bool isContained ( const HostIndexSet &indexSet ) const
      {
        return indexSet.contains( hostElement() );
      }

      /** \brief obtain the entity's id from a host IdSet
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param  idSet  host IdSet to use
       */
      template< class HostIdSet >
      typename HostIdSet::IdType id ( const HostIdSet &idSet ) const
      {
        return idSet.subId( hostElement(), subEntity_, codimension );
      }
      /** \} */

    private:
      PartitionType
      vertexPartitionType ( const GenericReferenceElement< ctype, dimension > &refElement, int i ) const
      {
        const int j = refElement.subEntity( subEntity_, codimension, 0, dimension );
        return hostElement().template subEntity< dimension >( j )->partitionType();
      }

    private:
      const Grid *grid_;
      const HostElement *hostElement_;
      ElementNumbering numbering_;
      unsigned int subEntity_;
      mutable Geometry geo_;
    };



    // EntityBase for codimension 0
    // ----------------------------

    /** \copydoc EntityBase
     *
     *  This specialization implements the extended interface for \em codim = 0.
     *  As all host grids provide this entity, we need only specialize the case
     *  \em fake = \b false.
     *
     *  \nosubgrouping
     */
    template< class Grid >
    class EntityBase< 0, Grid, false >
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

    public:
      /** \name Attributes
       *  \{ */

      //! codimensioon of the entity
      static const int codimension = 0;
      //! dimension of the grid
      static const int dimension = Traits::dimension;
      //! dimension of the entity
      static const int mydimension = dimension - codimension;
      //! dimension of the world
      static const int dimensionworld = Traits::dimensionworld;

      //! \b true, if the entity is faked, i.e., if there is no corresponding host entity
      static const bool fake = false;
      /** \} */

      /** \name Types Required by DUNE
       *  \{ */

      //! coordinate type of the grid
      typedef typename Traits::ctype ctype;

      //! type of corresponding geometry
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;
      //! type of corresponding local geometry
      typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;
      //! type of corresponding entity pointer
      typedef typename Traits::template Codim< codimension >::EntityPointer EntityPointer;

      //! type of hierarchic iterator
      typedef typename Traits::HierarchicIterator HierarchicIterator;
      //! type of leaf intersection iterator
      typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
      //! type of level intersection iterator
      typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

      /** \} */

    private:
      typedef typename Traits::HostGrid HostGrid;
      typedef typename Traits::CoordFunction CoordFunction;

      typedef typename Traits::template Codim< codimension >::EntityNumbering EntityNumbering;

    public:
      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGrid::template Codim< codimension > :: Entity HostEntity;
      //! type of corresponding host entity pointer
      typedef typename HostGrid::template Codim< codimension > :: EntityPointer HostEntityPointer;

      //! type of host elements, i.e., of host entities of codimension 0
      typedef typename HostGrid::template Codim< 0 >::Entity HostElement;
      /** \} */

    private:
      typedef typename HostGrid::template Codim< codimension >::Geometry HostGeometry;

      typedef typename GenericGeometry::GlobalGeometryTraits< Grid >::template Codim< codimension >::CoordVector
      CoordVector;

      typedef MakeableInterfaceObject< Geometry > MakeableGeometry;
      typedef typename MakeableGeometry::ImplementationType GeometryImpl;
      typedef typename GeometryImpl::GlobalCoordinate GlobalCoordinate;

    public:
      /** \name Construction, Initialization and Destruction
       *  \{ */

      /** \brief construct an initialized entity
       *
       *  \param[in]  grid        GeometryGrid this entity belongs to
       *  \param[in]  hostEntity  corresponding entity in the host grid
       *
       *  \note Both references must remain valid as long as this entity is in
       *        use.
       */
      EntityBase ( const Grid &grid, const HostEntity &hostEntity )
        : grid_( &grid ),
          hostEntity_( &hostEntity ),
          numbering_( grid.numbering()[ hostEntity ] ),
          geo_( GeometryImpl( grid.allocator() ) )
      {}

      EntityBase ( const EntityBase &other )
        : grid_( other.grid_ ),
          hostEntity_( other.hostEntity_ ),
          numbering_( other.numbering_ ),
          geo_( GeometryImpl( grid().allocator() ) )
      {}

      /** \} */

    private:
      EntityBase &operator= ( const EntityBase & );

    public:
      /** \name Methods Shared by Entities of All Codimensions
       *  \{ */

      /** \brief obtain the name of the corresponding reference element
       *
       *  This type can be used to access the DUNE reference element.
       */
      GeometryType type () const
      {
        return hostEntity().type();
      }

      unsigned int topologyId () const
      {
        return GenericGeometry::topologyId( type() );
      }

      /** \brief obtain the level of this entity */
      int level () const
      {
        return hostEntity().level();
      }

      /** \brief obtain the partition type of this entity */
      PartitionType partitionType () const
      {
        return hostEntity().partitionType();
      }

      /** obtain the geometry of this entity
       *
       *  Each DUNE entity encapsulates a geometry object, representing the map
       *  from the reference element to world coordinates. Wrapping the geometry
       *  is the main objective of the GeometryGrid.
       *
       *  The GeometryGrid provides geometries of order 1, obtained by
       *  interpolation of its corners \f$y_i\f$. There corners are calculated
       *  from the corners \f$x_i\f$ of the host geometry through the
       *  GeometryGrid's coordinate function \f$c\f$, i.e.,
       *  \f$y_i = c( x_i )\f$.
       *
       *  \returns a const reference to the geometry
       */
      const Geometry &geometry () const
      {
        GeometryImpl &geo = Grid::getRealImplementation( geo_ );

        if( !geo )
        {
          CoordVector coords( hostEntity(), grid().coordFunction() );
          geo = GeometryImpl( topologyId(), coords, grid().allocator() );
        }
        return geo_;
      }
      /** \} */

      template< int codim >
      int count () const
      {
        return hostEntity().template count< codim >();
      }

      template< int codim >
      typename Grid::template Codim< codim >::EntityPointer
      subEntity ( int i ) const
      {
        typedef typename Grid::template Codim< codim >::EntityPointer EntityPointer;
        typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
        typedef typename MakeableEntityPointer::ImplementationType EntityPointerImpl;

        EntityPointerImpl impl( *grid_, hostEntity(), i );
        return MakeableEntityPointer( impl );
      }

      LevelIntersectionIterator ilevelbegin () const
      {
        typedef MakeableInterfaceObject< LevelIntersectionIterator >
        MakeableLevelIntersectionIterator;
        typedef typename MakeableLevelIntersectionIterator :: ImplementationType
        LevelIntersectionIteratorImpl;
        LevelIntersectionIteratorImpl impl( *this, hostEntity().ilevelbegin() );
        return MakeableLevelIntersectionIterator( impl );
      }

      LevelIntersectionIterator ilevelend () const
      {
        typedef MakeableInterfaceObject< LevelIntersectionIterator >
        MakeableLevelIntersectionIterator;
        typedef typename MakeableLevelIntersectionIterator :: ImplementationType
        LevelIntersectionIteratorImpl;
        LevelIntersectionIteratorImpl impl( *this, hostEntity().ilevelend() );
        return MakeableLevelIntersectionIterator( impl );
      }

      LeafIntersectionIterator ileafbegin () const
      {
        typedef MakeableInterfaceObject< LeafIntersectionIterator >
        MakeableLeafIntersectionIterator;
        typedef typename MakeableLeafIntersectionIterator :: ImplementationType
        LeafIntersectionIteratorImpl;
        LeafIntersectionIteratorImpl impl( *this, hostEntity().ileafbegin() );
        return MakeableLeafIntersectionIterator( impl );
      }

      LeafIntersectionIterator ileafend () const
      {
        typedef MakeableInterfaceObject< LeafIntersectionIterator >
        MakeableLeafIntersectionIterator;
        typedef typename MakeableLeafIntersectionIterator :: ImplementationType
        LeafIntersectionIteratorImpl;
        LeafIntersectionIteratorImpl impl( *this, hostEntity().ileafend() );
        return MakeableLeafIntersectionIterator( impl );
      }

      bool hasBoundaryIntersections () const
      {
        return hostEntity().hasBoundaryIntersections();
      }

      bool isLeaf () const
      {
        return hostEntity().isLeaf();
      }

      EntityPointer father () const
      {
        typedef MakeableInterfaceObject< EntityPointer > MakeableEntityPointer;
        typedef typename MakeableEntityPointer :: ImplementationType EntityPointerImpl;
        return MakeableEntityPointer( EntityPointerImpl( *grid_, hostEntity().father() ) );
      }

      bool hasFather () const
      {
        return hostEntity().hasFather();
      }

      const LocalGeometry &geometryInFather () const
      {
        return hostEntity().geometryInFather();
      }

      HierarchicIterator hbegin ( int maxLevel ) const
      {
        typedef MakeableInterfaceObject< HierarchicIterator > MakeableIterator;
        typedef typename MakeableIterator :: ImplementationType Impl;
        Impl impl( *grid_, hostEntity().hbegin( maxLevel ) );
        return MakeableIterator( impl );
      }

      HierarchicIterator hend ( int maxLevel ) const
      {
        typedef MakeableInterfaceObject< HierarchicIterator > MakeableIterator;
        typedef typename MakeableIterator :: ImplementationType Impl;
        Impl impl( *grid_, hostEntity().hend( maxLevel ) );
        return MakeableIterator( impl );
      }

      bool isRegular () const
      {
        return hostEntity().isRegular();
      }

      bool isNew () const
      {
        return hostEntity().isNew();
      }

      bool mightVanish () const
      {
        return hostEntity().mightVanish();
      }

      /** \name Methods Supporting the Grid Implementation
       *  \{ */
      const Grid &grid () const
      {
        return *grid_;
      }

      const HostEntity &hostEntity () const
      {
        return *hostEntity_;
      }

      /** \brief obtain the entity's index from a host IndexSet
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param[in]  indexSet  host IndexSet to use
       */
      template< class HostIndexSet >
      typename HostIndexSet::IndexType index ( const HostIndexSet &indexSet ) const
      {
        return indexSet.template index< codimension >( hostEntity() );
      }

      /** \brief obtain the index of a subentity from a host IndexSet
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param[in]  indexSet  host IndexSet to use
       *  \param[in]  i         number of the subentity
       *  \param[in]  cd        codimension of the subentity
       */
      template< class HostIndexSet >
      typename HostIndexSet::IndexType
      subIndex ( const HostIndexSet &indexSet, int i, unsigned int cd ) const
      {
        return indexSet.subIndex( hostEntity(), i, cd );
      }

      /** \brief check whether the entity is contained in a host index set
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param  indexSet  host IndexSet to use
       */
      template< class HostIndexSet >
      bool isContained ( const HostIndexSet &indexSet ) const
      {
        return indexSet.contains( hostEntity() );
      }

      /** \brief obtain the entity's id from a host IdSet
       *
       *  \internal This method is provided by the entity, because its
       *  implementation is different for fake and non-fake entities.
       *
       *  \param  idSet  host IdSet to use
       */
      template< class HostIdSet >
      typename HostIdSet::IdType id ( const HostIdSet &idSet ) const
      {
        return idSet.template id< codimension >( hostEntity() );
      }

      /** \} */

    private:
      const Grid *grid_;
      const HostEntity *hostEntity_;
      EntityNumbering numbering_;
      mutable MakeableGeometry geo_;
    };



    template< int codim, int dim, class Grid >
    class Entity
      : public EntityBase< codim, Grid >
    {
      typedef EntityBase< codim, Grid > Base;

    public:
      typedef typename Base::HostEntity HostEntity;
      typedef typename Base::HostElement HostElement;

      Entity ( const Grid &grid, const HostEntity &hostEntity )
        : Base( grid, hostEntity )
      {}

      Entity ( const Grid &grid, const HostElement &hostElement, int subEntity )
        : Base( grid, hostElement, subEntity )
      {}
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_ENTITY_HH
