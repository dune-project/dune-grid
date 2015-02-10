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
    class ExportParams;




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
      typedef Dune::EntitySeed< const Grid, GeoGrid::EntitySeed< codimension, const Grid > > EntitySeed;

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

      EntityPointer ()
        : hostEntityIterator_()
        , grid_( nullptr )
      {}

      EntityPointer ( const Grid &grid, const HostEntityIterator &hostEntityIterator )
        : hostEntityIterator_( hostEntityIterator )
        , grid_( &grid )
      {}

      EntityPointer ( const Grid &grid, HostEntityIterator&& hostEntityIterator )
        : hostEntityIterator_( std::move( hostEntityIterator ) )
        , grid_( &grid )
      {}

      EntityPointer ( const Grid &grid, const HostElement &hostElement, int subEntity )
        : hostEntityIterator_( hostElement.template subEntity< codimension >( subEntity ) )
        , grid_( &grid )
      {}

      EntityPointer ( const Grid &grid, const EntitySeed &seed )
        : hostEntityIterator_( grid.hostGrid().entityPointer( grid.getRealImplementation(seed).hostEntitySeed() ) )
        , grid_( &grid )
      {}

      EntityPointer ( const This &other )
        : hostEntityIterator_( other.hostEntityIterator_ )
        , grid_( other.grid_ )
      {}

      EntityPointer ( This&& other )
        : hostEntityIterator_( std::move( other.hostEntityIterator_ ) )
        , grid_( other.grid_ )
      {}

      template< class T >
      explicit EntityPointer ( const EntityPointer< T, fake > &other )
        : hostEntityIterator_( other.hostEntityIterator_ )
        , grid_( other.grid_ )
      {}

      EntityPointer ( const EntityImpl &entity )
        : hostEntityIterator_( entity.hostEntity() )
        , grid_( &entity.grid() )
      {}

      This &operator= ( const This &other )
      {
        hostEntityIterator_ = other.hostEntityIterator_;
        grid_ = other.grid_;
        return *this;
      }

      This &operator= ( This&& other )
      {
        hostEntityIterator_ = std::move( other.hostEntityIterator_ );
        grid_ = other.grid_;
        return *this;
      }

      template< class T >
      This &operator= ( const EntityPointer< T, fake > &other )
      {
        hostEntityIterator_ = other.hostEntityIterator_;
        grid_ = other.grid_;
        return *this;
      }

      template< class T >
      bool equals ( const EntityPointer< T, fake > &other ) const
      {
        return (hostIterator() == other.hostIterator());
      }

      Entity dereference () const
      {
        return EntityImpl( grid(), *hostIterator() );
      }

      int level () const { return hostIterator().level(); }

      const HostEntityIterator &hostIterator() const { return hostEntityIterator_; }

      const Grid &grid () const { return *grid_; }

    protected:
      HostEntityIterator hostEntityIterator_;
      const Grid *grid_;
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

      EntityPointer ()
        : hostElementIterator_()
        , grid_( nullptr )
        , subEntity_ ( 0 )
      {}

      EntityPointer ( const Grid &grid, const HostElementIterator &hostElementIterator, int subEntity )
        : hostElementIterator_( hostElementIterator )
        , grid_( &grid )
        , subEntity_( subEntity )
      {}

      EntityPointer ( const Grid &grid, HostElementIterator&& hostElementIterator, int subEntity )
        : hostElementIterator_( std::move( hostElementIterator ) )
        , grid_( &grid )
        , subEntity_( subEntity )
      {}

      EntityPointer ( const Grid &grid, const HostElement &hostElement, int subEntity )
        : hostElementIterator_( hostElement )
        , grid_( &grid )
        , subEntity_( subEntity )
      {}

      EntityPointer ( const Grid &grid, const EntitySeed &seed )
        : hostElementIterator_( grid.hostGrid().entityPointer( grid.getRealImplementation(seed).hostElementSeed() ) )
        , grid_( &grid )
        , subEntity_( grid.getRealImplementation(seed).subEntity() )
      {}

      explicit EntityPointer ( const EntityImpl &entity )
        : hostElementIterator_( entity.hostElement() )
        , grid_( entity.grid() )
        , subEntity_( entity.subEntity() )
      {}

      EntityPointer ( const This &other )
        : hostElementIterator_( other.hostElementIterator_ )
        , grid_( other.grid_ )
        , subEntity_( other.subEntity_ )
      {}

      EntityPointer ( This&& other )
        : hostElementIterator_( std::move( other.hostElementIterator_ ) )
        , grid_( other.grid_ )
        , subEntity_( other.subEntity_ )
      {}

      template< class T >
      explicit EntityPointer ( const EntityPointer< T, fake > &other )
        : hostElementIterator_( other.hostElementIterator_ )
        , grid_( other.grid_ )
        , subEntity_( other.subEntity_ )
      {}

      This &operator= ( const This &other )
      {
        hostElementIterator_ = other.hostElementIterator_;
        grid_ = other.grid_;
        subEntity_ = other.subEntity_;
        return *this;
      }

      This &operator= ( This&& other )
      {
        hostElementIterator_ = std::move( other.hostElementIterator_ );
        grid_ = other.grid_;
        subEntity_ = other.subEntity_;
        return *this;
      }

      template< class T >
      This &operator= ( const EntityPointer< T, fake > &other )
      {
        hostElementIterator_ = other.hostElementIterator_;
        grid_ = other.grid_;
        subEntity_ = other.subEntity_;
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

      Entity dereference () const
      {
        return EntityImpl( grid(), *hostElementIterator(), subEntity() );
      }

      int level () const { return hostElementIterator()->level(); }

      const Grid &grid () const { return grid_; }
      int subEntity () const { return subEntity_; }

    protected:

      const HostElementIterator &hostElementIterator () const
      {
        return hostElementIterator_;
      }

    protected:
      HostElementIterator hostElementIterator_;
      const Grid* grid_;
      int subEntity_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_ENTITYPOINTER_HH
