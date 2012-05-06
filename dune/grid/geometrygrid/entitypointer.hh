// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITYPOINTER_HH
#define DUNE_GEOGRID_ENTITYPOINTER_HH

#include <dune/grid/common/grid.hh>

#include <dune/grid/geometrygrid/declaration.hh>
#include <dune/grid/geometrygrid/capabilities.hh>
#include <dune/grid/geometrygrid/entityseed.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // External Forward Declarations
    // -----------------------------

    template< int, int, class >
    class Entity;

    template< class, class >
    struct ExportParams;




    // Internal Forward Declarations
    // -----------------------------

    template< int codim, class Grid >
    struct EntityPointerTraits;

    template< class Traits, bool fake = Traits::fake >
    class EntityPointer;



    // EntityPointerTraits
    // -------------------

    template< int codim, class Grid >
    struct EntityPointerTraits;

    /** \cond */
    template< int codim, class Grid >
    struct EntityPointerTraits< codim, const Grid >
      : public EntityPointerTraits< codim, Grid >
    {};
    /** \endcond */

    template< int codim, class HostGrid, class CoordFunction, class Allocator >
    struct EntityPointerTraits< codim, GeometryGrid< HostGrid, CoordFunction, Allocator > >
      : public ExportParams< HostGrid, CoordFunction >
    {
      typedef Dune::GeometryGrid< HostGrid, CoordFunction, Allocator > Grid;

      static const bool fake = !Capabilities::hasHostEntity< Grid, codim >::v;

      typedef typename HostGrid::ctype ctype;

      static const int dimension = HostGrid::dimension;
      static const int codimension = codim;

      typedef Dune::Entity< codimension, dimension, const Grid, GeoGrid::Entity > Entity;
      typedef GeoGrid::EntitySeed< codimension, const Grid > EntitySeed;

      typedef typename HostGrid::template Codim< codim >::Entity HostEntity;
      typedef typename HostGrid::template Codim< codim >::EntityPointer HostEntityPointer;
      typedef HostEntityPointer HostEntityIterator;

      typedef typename HostGrid::template Codim< 0 >::Entity HostElement;
      typedef typename HostGrid::template Codim< 0 >::EntityPointer HostElementIterator;
    };



    // EntityPointer (real)
    // --------------------

    template< class Traits >
    class EntityPointer< Traits, false >
    {
      typedef EntityPointer< Traits, false > This;

      typedef typename Traits::Grid Grid;

      typedef EntityPointerTraits< Traits::codimension, const Grid > BaseTraits;
      friend class EntityPointer< BaseTraits, false >;

    public:
      static const int dimension = Traits::dimension;
      static const int codimension = Traits::codimension;

      typedef typename Traits::Entity Entity;

      static const bool fake = Traits::fake;

      typedef EntityPointer< BaseTraits, fake > EntityPointerImp;

    protected:
      typedef typename Traits::HostEntityPointer HostEntityPointer;
      typedef typename Traits::HostEntityIterator HostEntityIterator;
      typedef typename Traits::HostElement HostElement;

      typedef typename Traits::EntitySeed EntitySeed;

      typedef GeoGrid::Entity< codimension, dimension, const Grid > EntityImpl;
      typedef typename EntityImpl::GeometryImpl GeometryImpl;

    public:
      EntityPointer ( const GeometryImpl &geo, const HostEntityIterator &hostEntityIterator )
        : entity_( EntityImpl( geo ) ),
          hostEntityIterator_( hostEntityIterator )
      {}

      EntityPointer ( const Grid &grid, const HostEntityIterator &hostEntityIterator )
        : entity_( EntityImpl( grid ) ),
          hostEntityIterator_( hostEntityIterator )
      {}

      EntityPointer ( const Grid &grid, const HostElement &hostElement, int subEntity )
        : entity_( EntityImpl( grid ) ),
          hostEntityIterator_( hostElement.template subEntity< codimension >( subEntity ) )
      {}

      EntityPointer ( const Grid &grid, const EntitySeed &seed )
        : entity_( EntityImpl( grid ) ),
          hostEntityIterator_( grid.hostGrid().entityPointer( seed.hostEntitySeed() ) )
      {}

      explicit EntityPointer ( const EntityImpl &entity )
        : entity_( entity ),
          hostEntityIterator_( entity.hostEntity() )
      {}

      EntityPointer ( const This &other )
        : entity_( other.entityImpl() ),
          hostEntityIterator_( other.hostEntityIterator_ )
      {}

      template< class T >
      explicit EntityPointer ( const EntityPointer< T, fake > &other )
        : entity_( other.entityImpl() ),
          hostEntityIterator_( other.hostEntityIterator_ )
      {}

      const This &operator= ( const This &other )
      {
        entityImpl() = other.entityImpl();
        hostEntityIterator_ = other.hostEntityIterator_;
        return *this;
      }

      template< class T >
      bool equals ( const EntityPointer< T, fake > &other ) const
      {
        return (hostIterator() == other.hostIterator());
      }

      Entity &dereference () const
      {
        if( !entityImpl() )
          entityImpl().initialize( *hostIterator() );
        return entity_;
      }

      int level () const { return hostIterator().level(); }

      const HostEntityIterator &hostIterator() const { return hostEntityIterator_; }

      const Grid &grid () const { return entityImpl().grid(); }

    protected:
      EntityImpl &entityImpl () const
      {
        return Grid::getRealImplementation( entity_ );
      }

    private:
      mutable Entity entity_;

    protected:
      HostEntityIterator hostEntityIterator_;
    };



    // EntityPointer (fake)
    // --------------------

    template< class Traits >
    class EntityPointer< Traits, true >
    {
      typedef EntityPointer< Traits, true > This;

      typedef typename Traits::Grid Grid;

      typedef EntityPointerTraits< Traits::codimension, const Grid > BaseTraits;
      friend class EntityPointer< BaseTraits, true >;

    public:
      static const int dimension = Traits::dimension;
      static const int codimension = Traits::codimension;

      typedef typename Traits::Entity Entity;

      static const bool fake = Traits::fake;

      typedef EntityPointer< BaseTraits, fake > EntityPointerImp;

    protected:
      typedef typename Traits::HostEntityPointer HostEntityPointer;
      typedef typename Traits::HostElementIterator HostElementIterator;
      typedef typename Traits::HostElement HostElement;

      typedef typename Traits::EntitySeed EntitySeed;

      typedef GeoGrid::Entity< codimension, dimension, const Grid > EntityImpl;
      typedef typename EntityImpl::GeometryImpl GeometryImpl;

    public:
      EntityPointer ( const GeometryImpl &geo, const HostElementIterator &hostElementIterator, int subEntity )
        : entity_( EntityImpl( geo, subEntity ) ),
          hostElementIterator_( hostElementIterator )
      {}

      EntityPointer ( const Grid &grid, const HostElementIterator &hostElementIterator, int subEntity )
        : entity_( EntityImpl( grid, subEntity ) ),
          hostElementIterator_( hostElementIterator )
      {}

      EntityPointer ( const Grid &grid, const HostElement &hostElement, int subEntity )
        : entity_( EntityImpl( grid, subEntity ) ),
          hostElementIterator_( hostElement )
      {}

      EntityPointer ( const Grid &grid, const EntitySeed &seed )
        : entity_( EntityImpl( grid, seed.subEntity() ) ),
          hostElementIterator_( grid.hostGrid().entityPointer( seed.hostElementSeed() ) )
      {}

      explicit EntityPointer ( const EntityImpl &entity )
        : entity_( entity ),
          hostElementIterator_( entity.hostElement() )
      {}

      EntityPointer ( const This &other )
        : entity_( other.entityImpl() ),
          hostElementIterator_( other.hostElementIterator_ )
      {}

      template< class T >
      explicit EntityPointer ( const EntityPointer< T, fake > &other )
        : entity_( other.entityImpl() ),
          hostElementIterator_( other.hostElementIterator_ )
      {}

      const This &operator= ( const This &other )
      {
        entityImpl() = other.entityImpl();
        hostElementIterator_ = other.hostElementIterator_;
        return *this;
      }

      template< class T >
      bool equals ( const EntityPointer< T, fake > &other ) const
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

        const HostElement &thisElement = *hostElementIterator();
        assert( indexSet.contains( thisElement ) );
        const HostElement &otherElement = *(other.hostElementIterator());
        assert( indexSet.contains( otherElement ) );

        const int thisIndex = indexSet.subIndex( thisElement, subEntity(), codimension );
        const int otherIndex = indexSet.subIndex( otherElement, other.subEntity(), codimension );
        return (thisIndex == otherIndex);
      }

      Entity &dereference () const
      {
        if( !entityImpl() )
          entityImpl().initialize( *hostElementIterator() );
        return entity_;
      }

      int level () const { return hostElementIterator().level(); }

      const Grid &grid () const { return entityImpl().grid(); }
      int subEntity () const { return entityImpl().subEntity(); }

    protected:
      EntityImpl &entityImpl () const
      {
        return Grid::getRealImplementation( entity_ );
      }

      const HostElementIterator &hostElementIterator () const
      {
        return hostElementIterator_;
      }

    private:
      mutable Entity entity_;

    protected:
      HostElementIterator hostElementIterator_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_ENTITYPOINTER_HH
