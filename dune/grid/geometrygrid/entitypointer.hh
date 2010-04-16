// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITYPOINTER_HH
#define DUNE_GEOGRID_ENTITYPOINTER_HH

#include <dune/grid/common/grid.hh>

#include <dune/grid/geometrygrid/capabilities.hh>
#include <dune/grid/geometrygrid/storage.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CordFunction >
  class GeometryGrid;



  namespace GeoGrid
  {

    // External Forward Declarations
    // -----------------------------

    template< int, int, class >
    class Entity;

    template< class >
    class EntityWrapper;

    template< class HostGrid, class CoordFunction >
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

    template< int codim, class HostGrid, class CoordFunction >
    struct EntityPointerTraits< codim, GeometryGrid< HostGrid, CoordFunction > >
      : public ExportParams< HostGrid, CoordFunction >
    {
      typedef Dune::GeometryGrid< HostGrid, CoordFunction > Grid;

      static const bool fake = !Capabilities::hasHostEntity< Grid, codim >::v;

      typedef typename HostGrid::ctype ctype;

      static const int dimension = HostGrid::dimension;
      static const int codimension = codim;

      typedef Dune::Entity< codimension, dimension, const Grid, GeoGrid::Entity > Entity;

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

      typedef EntityPointerTraits< Traits :: codimension, const Grid > BaseTraits;
      friend class EntityPointer< BaseTraits, false >;

    public:
      static const int dimension = Traits::dimension;
      static const int codimension = Traits::codimension;

      typedef typename Traits::Entity Entity;

      static const bool fake = Traits::fake;

      typedef EntityPointer< BaseTraits, fake > EntityPointerImp;

    private:
      typedef GeoGrid::EntityWrapper< Entity > EntityWrapper;
      typedef GeoGrid::Storage< EntityWrapper > EntityStorage;

      const Grid *grid_;
      mutable EntityWrapper *entity_;

    protected:
      typedef typename Traits::HostEntityPointer HostEntityPointer;
      typedef typename Traits::HostEntityIterator HostEntityIterator;
      typedef typename Traits::HostElement HostElement;

      HostEntityIterator hostEntityIterator_;

    public:
      EntityPointer ( const Grid &grid, const HostEntityIterator &hostEntityIterator )
        : grid_( &grid ),
          entity_( 0 ),
          hostEntityIterator_( hostEntityIterator )
      {}

      EntityPointer ( const Grid &grid, const HostElement &hostElement, int subEntity )
        : grid_( &grid ),
          entity_( 0 ),
          hostEntityIterator_( hostElement.template subEntity< codimension >( subEntity ) )
      {}

      EntityPointer ( const typename EntityWrapper::Implementation &entity )
        : grid_( &entity.grid() ),
          entity_( 0 ),
          hostEntityIterator_( entity.hostEntity() )
      {}

      EntityPointer ( const This &other )
        : grid_( other.grid_ ),
          entity_( 0 ),
          hostEntityIterator_( other.hostEntityIterator_ )
      {}

      template< class T >
      explicit EntityPointer ( const EntityPointer< T, fake > &other )
        : grid_( other.grid_ ),
          entity_( 0 ),
          hostEntityIterator_( other.hostEntityIterator_ )
      {}

      ~EntityPointer ()
      {
        EntityStorage::free( entity_ );
      }

      This &operator= ( const This &other )
      {
        grid_ = other.grid_;
        hostEntityIterator_ = other.hostEntityIterator_;
        releaseEntity();
        return *this;
      }

      operator const EntityPointerImp & () const
      {
        return reinterpret_cast< const EntityPointerImp & >( *this );
      }

      template< class T >
      bool equals ( const EntityPointer< T, fake > &other ) const
      {
        return (hostIterator() == other.hostIterator());
      }

      Entity &dereference () const
      {
        if( entity_ == 0 )
        {
          entity_ = EntityStorage::alloc();
          entity_->initialize( grid(), *hostIterator() );
        }
        return *entity_;
      }

      int level () const
      {
        return hostIterator().level();
      }

      void compactify ()
      {
        hostEntityIterator_.compactify();
        releaseEntity();
      }

      const HostEntityIterator &hostIterator() const
      {
        return hostEntityIterator_;
      }

      const Grid &grid () const
      {
        return *grid_;
      }

    protected:
      void releaseEntity ()
      {
        EntityStorage::free( entity_ );
        entity_ = 0;
      }
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

    private:
      typedef GeoGrid::EntityWrapper< Entity > EntityWrapper;
      typedef GeoGrid::Storage< EntityWrapper > EntityStorage;

      const Grid *grid_;
      mutable EntityWrapper *entity_;

    protected:
      typedef typename Traits::HostEntityPointer HostEntityPointer;
      typedef typename Traits::HostElementIterator HostElementIterator;
      typedef typename Traits::HostElement HostElement;

      int subEntity_;
      HostElementIterator hostElementIterator_;

    public:
      EntityPointer ( const Grid &grid,
                      const HostElementIterator &hostElementIterator, int subEntity )
        : grid_( &grid ),
          entity_( 0 ),
          subEntity_( subEntity ),
          hostElementIterator_( hostElementIterator )
      {}

      EntityPointer ( const Grid &grid, const HostElement &hostElement, int subEntity )
        : grid_( &grid ),
          entity_( 0 ),
          subEntity_( subEntity ),
          hostElementIterator_( hostElement )
      {}

      EntityPointer ( const typename EntityWrapper::Implementation &entity )
        : grid_( &entity.grid() ),
          entity_( 0 ),
          subEntity_( entity.subEntity() ),
          hostElementIterator_( entity.hostElement() )
      {}

      EntityPointer ( const This &other )
        : grid_( other.grid_ ),
          entity_( 0 ),
          subEntity_( other.subEntity_ ),
          hostElementIterator_( other.hostElementIterator_ )
      {}

      template< class T >
      explicit EntityPointer ( const EntityPointer< T, fake > &other )
        : grid_( other.grid_ ),
          entity_( 0 ),
          subEntity_( other.subEntity_ ),
          hostElementIterator_( other.hostElementIterator_ )
      {}

      ~EntityPointer ()
      {
        EntityStorage::free( entity_ );
      }

      This &operator= ( const This &other )
      {
        grid_ = other.grid_;
        subEntity_ = other.subEntity_;
        hostElementIterator_ = other.hostElementIterator_;
        releaseEntity();
        return *this;
      }

      operator const EntityPointerImp & () const
      {
        return reinterpret_cast< const EntityPointerImp & >( *this );
      }

      operator EntityPointerImp & ()
      {
        return reinterpret_cast< EntityPointerImp & >( *this );
      }

      template< class T >
      bool equals ( const EntityPointer< T, fake > &other ) const
      {
        const bool thisEnd = (subEntity_ < 0);
        const bool otherEnd = (other.subEntity_ < 0);
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

        const int thisIndex = indexSet.subIndex( thisElement, subEntity_, codimension );
        const int otherIndex = indexSet.subIndex( otherElement, other.subEntity_, codimension );
        return (thisIndex == otherIndex);
      }

      Entity &dereference () const
      {
        if( entity_ == 0 )
        {
          entity_ = EntityStorage::alloc();
          entity_->initialize( grid(), *hostElementIterator(), subEntity_ );
        }
        return *entity_;
      }

      int level () const
      {
        return hostElementIterator().level();
      }

      void compactify ()
      {
        hostElementIterator_.compactify();
        releaseEntity();
      }

      const Grid &grid () const
      {
        return *grid_;
      }

    protected:
      void releaseEntity ()
      {
        EntityStorage::free( entity_ );
        entity_ = 0;
      }

      const HostElementIterator &hostElementIterator () const
      {
        return hostElementIterator_;
      }
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_ENTITYPOINTER_HH
