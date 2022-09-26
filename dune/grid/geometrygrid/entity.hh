// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITY_HH
#define DUNE_GEOGRID_ENTITY_HH

#include <dune/geometry/referenceelements.hh>

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
     *                  (do not specify, the default value is already the
     *                  intended use)
     */
    template< int codim, class Grid, bool fake = !(Capabilities::hasHostEntity< Grid, codim >::v) >
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



    // External Forward Declarations
    // -----------------------------

    template< class Grid >
    class HierarchicIterator;

    template< class Grid, class HostIntersectionIterator >
    class IntersectionIterator;



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
      typedef typename std::remove_const< Grid >::type::Traits Traits;

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

    public:
      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;

      //! type of corresponding entity seed
      typedef typename Traits::template Codim< codimension >::EntitySeed EntitySeed;

      //! type of host elements, i.e., of host entities of codimension 0
      typedef typename HostGrid::template Codim< 0 >::Entity HostElement;
      /** \} */

      typedef typename Traits::template Codim< codim >::GeometryImpl GeometryImpl;

    private:
      typedef typename HostGrid::template Codim< codimension >::Geometry HostGeometry;

      typedef GeoGrid::CoordVector< mydimension, Grid, fake > CoordVector;

    public:
      /** \name Construction, Initialization and Destruction
       *  \{ */

      EntityBase ()
        : hostEntity_()
        , grid_( nullptr )
        , geo_()
      {}

      EntityBase ( const Grid &grid, const EntitySeed &seed )
        : hostEntity_( grid.hostGrid().entity( seed.impl().hostEntitySeed() ) )
        , grid_( &grid )
      {}

      EntityBase ( const Grid &grid, const HostElement &hostElement, int i )
        : hostEntity_( hostElement.template subEntity<codim>(i) )
        , grid_( &grid )
      {}


      EntityBase ( const GeometryImpl &geo, const HostEntity &hostEntity )
        : hostEntity_( hostEntity )
        , grid_( &geo.grid() )
        , geo_( geo )
      {}

      EntityBase ( const GeometryImpl &geo, HostEntity&& hostEntity )
        : hostEntity_( std::move( hostEntity ) )
        , grid_( &geo.grid() )
        , geo_( geo )
      {}

      EntityBase ( const Grid &grid, const HostEntity& hostEntity )
        : hostEntity_( hostEntity )
        , grid_( &grid )
      {}

      EntityBase ( const Grid &grid, HostEntity&& hostEntity )
        : hostEntity_( std::move( hostEntity ) )
        , grid_( &grid )
      {}


      EntityBase ( const EntityBase &other )
        : hostEntity_( other.hostEntity_ )
        , grid_( other.grid_ )
        , geo_( other.geo_ )
      {}

      EntityBase ( EntityBase&& other )
        : hostEntity_( std::move( other.hostEntity_ ) )
        , grid_( other.grid_ )
        , geo_( std::move( other.geo_ ) )
      {}

      /** \} */

      const EntityBase &operator= ( const EntityBase &other )
      {
        hostEntity_ = other.hostEntity_;
        grid_ = other.grid_;
        geo_ = other.geo_;
        return *this;
      }

      const EntityBase &operator= ( EntityBase&& other )
      {
        hostEntity_ = std::move( other.hostEntity_ );
        grid_ = std::move( other.grid_ );
        geo_ = std::move( other.geo_ );
        return *this;
      }

      /** \brief compare two entities */
      bool equals ( const EntityBase &other) const
      {
        return hostEntity_ == other.hostEntity_;
      }

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
      Geometry geometry () const
      {
        if( !geo_ )
        {
          CoordVector coords( hostEntity(), grid().coordFunction() );
          geo_ = GeometryImpl( grid(), type(), coords );
        }
        return Geometry( geo_ );
      }

      unsigned int subEntities ( unsigned int cc ) const
      {
        return hostEntity().subEntities( cc );
      }

      /** \brief return EntitySeed of host grid entity */
      EntitySeed seed () const { return typename EntitySeed::Implementation( hostEntity().seed() ); }
      /** \} */


      /** \name Methods Supporting the Grid Implementation
       *  \{ */

      const Grid &grid () const { assert( grid_ ); return *grid_; }

      const HostEntity &hostEntity () const
      {
        return hostEntity_;
      }

      /** \brief initiliaze an entity
       *
       *  \param[in]  hostEntity  reference to the host entity
       *
       */
      void initialize ( const HostEntity &hostEntity ) { hostEntity_ = hostEntity; }

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
      HostEntity hostEntity_;
      const Grid *grid_;
      mutable GeometryImpl geo_;
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
      typedef typename std::remove_const< Grid >::type::Traits Traits;

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

    public:
      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;

      //! type of corresponding entity seed
      typedef typename Traits::template Codim< codimension >::EntitySeed EntitySeed;

      //! type of host elements, i.e., of host entities of codimension 0
      typedef typename HostGrid::template Codim< 0 >::Entity HostElement;
      /** \} */

      typedef typename Traits::template Codim< codimension >::GeometryImpl GeometryImpl;

    private:
      typedef typename HostGrid::template Codim< 0 >::Geometry HostGeometry;

      typedef GeoGrid::CoordVector< mydimension, Grid, fake > CoordVector;

    public:
      /** \name Construction, Initialization and Destruction
       *  \{ */

      EntityBase ()
        : hostElement_()
        , subEntity_(-1)
        , grid_(nullptr)
        , geo_()
      {}

      EntityBase(const Grid& grid, const HostElement& hostElement, unsigned int subEntity)
        : hostElement_(hostElement)
        , subEntity_(subEntity)
        , grid_(&grid)
      {}

      EntityBase ( const Grid &grid, const EntitySeed &seed )
        : hostElement_( grid.hostGrid().entity( seed.impl().hostElementSeed() ) )
        , subEntity_( seed.impl().subEntity() )
        , grid_( &grid )
      {}

      EntityBase ( const EntityBase &other )
        : hostElement_( other.hostElement_ )
        , subEntity_( other.subEntity_ )
        , grid_(other.grid_)
        , geo_( other.geo_ )
      {}

      EntityBase ( EntityBase &&other )
        : hostElement_( std::move( other.hostElement_ ) )
        , subEntity_( std::move( other.subEntity_ ) )
        , grid_( std::move( other.grid_ ) )
        , geo_( std::move( other.geo_ ) )
      {}

      /*
       * This method is required by constructors in the `Entity` class
       * below, however it cannot do anything useful for fake
       * entities.
       */
      EntityBase(const Grid& grid, const HostEntity& hostEntity)
      {
        DUNE_THROW(Dune::Exception, "GeometryGrid: Cannot create fake entity of codim " << codimension << " from real host entity.");
      }

      /** \} */

      const EntityBase &operator= ( const EntityBase &other )
      {
        hostElement_ = other.hostElement_;
        subEntity_ = other.subEntity_;
        grid_ = other.grid_;
        geo_ = other.geo_;
        return *this;
      }

      const EntityBase &operator= ( EntityBase&& other )
      {
        hostElement_ = std::move( other.hostElement_ );
        subEntity_ = std::move( other.subEntity_ );
        grid_ = std::move( other.grid_ );
        geo_ = std::move( other.geo_ );
        return *this;
      }

      /** \brief compare two entities */
      bool equals ( const EntityBase &other) const
      {
        const bool thisEnd = (subEntity() < 0);
        const bool otherEnd = (other.subEntity() < 0);
        if( thisEnd || otherEnd )
          return thisEnd && otherEnd;

        const int lvl = level();
        if( lvl != other.level() )
          return false;

        const typename Traits::HostGrid::Traits::LevelIndexSet &indexSet
          = grid().hostGrid().levelIndexSet( lvl );

        const HostElement &thisElement = hostElement();
        assert( indexSet.contains( thisElement ) );
        const HostElement &otherElement = other.hostElement();
        assert( indexSet.contains( otherElement ) );

        const int thisIndex = indexSet.subIndex( thisElement, subEntity(), codimension );
        const int otherIndex = indexSet.subIndex( otherElement, other.subEntity(), codimension );
        return (thisIndex == otherIndex);
      }

      /** \name Methods Shared by Entities of All Codimensions
       *  \{ */

      /** \brief obtain the name of the corresponding reference element
       *
       *  This type can be used to access the DUNE reference element.
       */
      GeometryType type () const
      {
        auto refElement = referenceElement< ctype, dimension >( hostElement().type() );
        return refElement.type( subEntity_, codimension );
      }

      /** \brief obtain the level of this entity */
      int level () const
      {
        return hostElement().level();
      }

      /** \brief obtain the partition type of this entity */
      PartitionType partitionType () const
      {
        auto refElement = referenceElement< ctype, dimension >( hostElement().type() );

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
       *  \returns a geometry object
       */
      Geometry geometry () const
      {
        if( !geo_ )
        {
          CoordVector coords( hostElement(), subEntity_, grid().coordFunction() );
          geo_ = GeometryImpl( grid(), type(), coords );
        }
        return Geometry( geo_ );
      }

      unsigned int subEntities ( unsigned int cc ) const
      {
        auto refElement = referenceElement< ctype, dimension >( hostElement().type() );
        return refElement.size( subEntity_, codimension, cc );
      }

      /** \brief return EntitySeed of host grid entity */
      EntitySeed seed () const { return typename EntitySeed::Implementation( hostElement().seed(), subEntity_ ); }
      /** \} */

      /** \name Methods Supporting the Grid Implementation
       *  \{ */

      const Grid &grid () const { assert( grid_ ); return *grid_; }

      const HostEntity &hostEntity () const
      {
        DUNE_THROW( NotImplemented, "HostGrid has no entities of codimension " << codimension << "." );
      }

      const HostElement &hostElement () const
      {
        return hostElement_;
      }

      int subEntity () const { return subEntity_; }

      /** \brief initiliaze an entity
       *
       *  \param[in]  hostElement  reference to the host element
       *
       *  \note The reference must remain valid as long as this entity is in
       *        use.
       */
      void initialize ( const HostElement &hostElement ) { hostElement_ = hostElement; }

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
        auto refElement = referenceElement< ctype, dimension >( hostElement().type() );
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
      vertexPartitionType ( Dune::Transitional::ReferenceElement< ctype, Dim<dimension> > refElement, int i ) const
      {
        const int j = refElement.subEntity( subEntity_, codimension, i, dimension );
        return hostElement().template subEntity< dimension >( j ).partitionType();
      }

    private:
      HostElement hostElement_;
      unsigned int subEntity_;
      const Grid *grid_;
      mutable GeometryImpl geo_;
    };



    // Entity
    // ------

    template< int codim, int dim, class Grid >
    class Entity
      : public EntityBase< codim, Grid >
    {
      typedef EntityBase< codim, Grid > Base;

    public:
      typedef typename Base::HostEntity HostEntity;
      typedef typename Base::HostElement HostElement;
      typedef typename Base::GeometryImpl GeometryImpl;
      typedef typename Base::EntitySeed EntitySeed;

      Entity () : Base() {}

      Entity ( const Grid &grid, const EntitySeed &seed ) : Base( grid, seed ) {}

      Entity ( const Grid &grid, const HostEntity &hostEntity ) : Base( grid, hostEntity ) {}
      Entity ( const Grid &grid, HostEntity&& hostEntity ) : Base( grid, std::move( hostEntity ) ) {}

      Entity ( const Grid &grid, const HostElement &hostEntity, int i ) : Base( grid, hostEntity, i ) {}

    };



    // Entity for codimension 0
    // ------------------------

    template< int dim, class Grid >
    class Entity< 0, dim, Grid >
      : public EntityBase< 0, Grid >
    {
      typedef EntityBase< 0, Grid > Base;

      typedef typename std::remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::HostGrid HostGrid;

    public:
      /** \name Attributes
       *  \{ */

      //! codimensioon of the entity
      static const int codimension = Base::codimension;
      //! dimension of the grid
      static const int dimension = Base::dimension;
      //! dimension of the entity
      static const int mydimension = Base::mydimension;
      //! dimension of the world
      static const int dimensionworld = Base::dimensionworld;

      //! \b true, if the entity is faked, i.e., if there is no corresponding host entity
      static const bool fake = Base::fake;
      /** \} */

      /** \name Types Required by DUNE
       *  \{ */

      //! type of corresponding local geometry
      typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;

      typedef Dune::Entity< 0, dim, Grid, Dune::GeoGrid::Entity > EntityFacade;

      //! type of hierarchic iterator
      typedef typename Traits::HierarchicIterator HierarchicIterator;
      //! type of leaf intersection iterator
      typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
      //! type of level intersection iterator
      typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

      /** \} */

      typedef typename Base::HostEntity HostEntity;
      typedef typename Base::HostElement HostElement;
      typedef typename Base::GeometryImpl GeometryImpl;
      typedef typename Base::EntitySeed EntitySeed;

      using Base::grid;
      using Base::hostEntity;

      Entity () : Base() {}

      Entity ( const Grid &g, const HostEntity &hostE ) : Base( g, hostE ) {}
      Entity ( const Grid &g, HostEntity&& hostE ) : Base( g, std::move( hostE ) ) {}
      Entity ( const GeometryImpl &geo, const HostEntity& hostE ) : Base( geo, hostE ) {}
      Entity ( const GeometryImpl &geo, HostEntity &&hostE ) : Base( geo, std::move( hostE ) ) {}

      Entity ( const Grid &g, const EntitySeed &seed ) : Base( g, seed ) {}

      Entity ( const Grid &g, const HostEntity &hostE, int i ) : Base( g, hostE )
      {
        assert( i == 0 );
      }

      template< int codim >
      typename Grid::template Codim< codim >::Entity
      subEntity ( int i ) const
      {
        typedef typename Traits::template Codim< codim >::EntityImpl EntityImpl;
        return EntityImpl( grid(), hostEntity(), i );
      }

      LevelIntersectionIterator ilevelbegin () const
      {
        typedef GeoGrid::IntersectionIterator< Grid, typename HostGrid::LevelIntersectionIterator > LevelIntersectionIteratorImpl;
        return LevelIntersectionIteratorImpl( *this, hostEntity().ilevelbegin() );
      }

      LevelIntersectionIterator ilevelend () const
      {
        typedef GeoGrid::IntersectionIterator< Grid, typename HostGrid::LevelIntersectionIterator > LevelIntersectionIteratorImpl;
        return LevelIntersectionIteratorImpl( *this, hostEntity().ilevelend() );
      }

      LeafIntersectionIterator ileafbegin () const
      {
        typedef GeoGrid::IntersectionIterator< Grid, typename HostGrid::LeafIntersectionIterator > LeafIntersectionIteratorImpl;
        return LeafIntersectionIteratorImpl( *this, hostEntity().ileafbegin() );
      }

      LeafIntersectionIterator ileafend () const
      {
        typedef GeoGrid::IntersectionIterator< Grid, typename HostGrid::LeafIntersectionIterator > LeafIntersectionIteratorImpl;
        return LeafIntersectionIteratorImpl( *this, hostEntity().ileafend() );
      }

      bool hasBoundaryIntersections () const
      {
        return hostEntity().hasBoundaryIntersections();
      }

      bool isLeaf () const
      {
        return hostEntity().isLeaf();
      }

      EntityFacade father () const
      {
        return Entity( grid(), hostEntity().father() );
      }

      bool hasFather () const
      {
        return hostEntity().hasFather();
      }

      LocalGeometry geometryInFather () const
      {
        return hostEntity().geometryInFather();
      }

      HierarchicIterator hbegin ( int maxLevel ) const
      {
        typedef GeoGrid::HierarchicIterator< Grid > HierarchicIteratorImpl;
        return HierarchicIteratorImpl( grid(), hostEntity().hbegin( maxLevel ) );
      }

      HierarchicIterator hend ( int maxLevel ) const
      {
        typedef GeoGrid::HierarchicIterator< Grid > HierarchicIteratorImpl;
        return HierarchicIteratorImpl( grid(), hostEntity().hend( maxLevel ) );
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
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_ENTITY_HH
