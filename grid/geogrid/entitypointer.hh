// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITYPOINTER_HH
#define DUNE_GEOGRID_ENTITYPOINTER_HH

#include <dune/grid/common/grid.hh>

#include <dune/grid/geogrid/capabilities.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int, int, class >
  class GeometryGridEntityAdapter;

  template< class HostGrid, class CoordFunction >
  struct GeometryGridExportParams;

  template< class HostGrid, class CordFunction >
  class GeometryGrid;



  // Internal Forward Declarations
  // -----------------------------

  template< int codim, class Grid >
  struct GeometryGridEntityPointerTraits;

  template< class Traits, bool fake = Traits :: fake >
  class GeometryGridEntityPointer;



  /************************************************************************
  * Warning:
  *
  * The wrapped iterators are not directly derived from the entity
  * pointer. They differ in the type of host iterator they store
  * internally.
  * This means that Iterator::Base is not binary compatible with
  * EntityPointer. Since the DUNE entity pointer interface contains a
  * dirty reinterpret_cast, this may cause trouble.
  * In this case, the fields in the entity pointer can be reordered such
  * that the wrapped iterator is the last field. Then, the copy
  * constructor for arbitrary entity pointer has to be marked explicit and
  * a two cast operators of the form
  *
  * operator Base & ()
  * {
  *   return reinterpret_cast< Base & >( *this );
  * }
  *
  * have to be added.
  *
  * This should only be done in the emergency that this code causes a
  * segmentation fault!
  ************************************************************************/



  // GeometryGridEntityPointerTraits
  // -------------------------------

  template< int codim, class Grid >
  struct GeometryGridEntityPointerTraits;

  /** \cond */
  template< int codim, class Grid >
  struct GeometryGridEntityPointerTraits< codim, const Grid >
    : public GeometryGridEntityPointerTraits< codim, Grid >
  {};
  /** \endcond */

  template< int codim, class HostGrid, class CoordFunction >
  struct GeometryGridEntityPointerTraits
  < codim, GeometryGrid< HostGrid, CoordFunction > >
    : public GeometryGridExportParams< HostGrid, CoordFunction >
  {
    typedef Dune :: GeometryGrid< HostGrid, CoordFunction > Grid;

    static const bool fake = !Capabilities :: hasHostEntity< Grid, codim > :: v;

    typedef typename HostGrid :: ctype ctype;

    static const int dimension = HostGrid :: dimension;
    static const int codimension = codim;

    typedef Dune :: Entity
    < codimension, dimension, const Grid, GeometryGridEntityAdapter >
    Entity;

    typedef typename HostGrid :: template Codim< codim > :: Entity HostEntity;
    typedef typename HostGrid :: template Codim< codim > :: EntityPointer
    HostEntityPointer;
    typedef HostEntityPointer HostEntityIterator;

    typedef typename HostGrid :: template Codim< 0 > :: Entity HostElement;
    typedef typename HostGrid :: template Codim< 0 > :: EntityPointer
    HostElementPointer;
    typedef HostElementPointer HostElementIterator;
  };



  // GeometryGridEntityPointer (real)
  // --------------------------------

  template< class Traits >
  class GeometryGridEntityPointer< Traits, false >
  {
    typedef GeometryGridEntityPointer< Traits, false > This;

    typedef typename Traits :: Grid Grid;

    typedef GeometryGridEntityPointerTraits< Traits :: codimension, const Grid >
    BaseTraits;
    friend class GeometryGridEntityPointer< BaseTraits, false >;

  public:
    static const int dimension = Traits :: dimension;
    static const int codimension = Traits :: codimension;

    typedef typename Traits :: Entity Entity;

    static const bool fake = Traits :: fake;

    typedef GeometryGridEntityPointer< BaseTraits, fake > Base;
    typedef GeometryGridEntityPointer< BaseTraits, fake > base;

  protected:
    typedef typename Traits :: HostEntityPointer HostEntityPointer;
    typedef typename Traits :: HostEntityIterator HostEntityIterator;
    typedef typename Traits :: HostElement HostElement;

    HostEntityIterator hostEntityIterator_;

  private:
    typedef MakeableInterfaceObject< Entity > MakeableEntity;
    typedef typename MakeableEntity :: ImplementationType EntityImpl;

    mutable MakeableEntity virtualEntity_;

  public:
    GeometryGridEntityPointer ( const Grid &grid,
                                const HostEntityIterator &hostEntityIterator )
      : hostEntityIterator_( hostEntityIterator ),
        virtualEntity_( EntityImpl( grid ) )
    {}

    GeometryGridEntityPointer ( const Grid &grid,
                                const HostElement &hostElement,
                                int subEntity )
      : hostEntityIterator_( hostElement.template entity< codimension >( subEntity ) ),
        virtualEntity_( EntityImpl( grid ) )
    {}

    GeometryGridEntityPointer ( const This &other )
      : hostEntityIterator_( other.hostEntityIterator_ ),
        virtualEntity_( other.grid() )
    {}

    template< class T >
    GeometryGridEntityPointer ( const GeometryGridEntityPointer< T, fake > &other )
      : hostEntityIterator_( other.hostEntityIterator_ ),
        virtualEntity_( other.grid() )
    {}

    This &operator= ( const This &other )
    {
      hostEntityIterator_ = other.hostEntityIterator_;
      update();
      return *this;
    }

    template< class T >
    bool equals ( const GeometryGridEntityPointer< T, fake > &other ) const
    {
      return (hostEntityPointer() == other.hostEntityPointer());
    }

    Entity &dereference () const
    {
      EntityImpl &impl = Grid :: getRealImplementation( virtualEntity_ );
      if( !impl.isValid() )
        impl.setToTarget( *hostEntityPointer() );
      return virtualEntity_;
    }

    int level () const
    {
      return hostEntityPointer().level();
    }

    const HostEntityPointer &hostEntityPointer () const
    {
      return hostEntityIterator_;
    }

  protected:
    const Grid &grid () const
    {
      return Grid :: getRealImplementation( virtualEntity_ ).grid();
    }

    void update ()
    {
      Grid :: getRealImplementation( virtualEntity_ ).invalidate();
    }
  };



  // GeometryGridEntityPointer (fake)
  // --------------------------------

  template< class Traits >
  class GeometryGridEntityPointer< Traits, true >
  {
    typedef GeometryGridEntityPointer< Traits, true > This;

    typedef typename Traits :: Grid Grid;

    typedef GeometryGridEntityPointerTraits< Traits :: codimension, const Grid >
    BaseTraits;
    friend class GeometryGridEntityPointer< BaseTraits, true >;

  public:
    static const int dimension = Traits :: dimension;
    static const int codimension = Traits :: codimension;

    typedef typename Traits :: Entity Entity;

    static const bool fake = Traits :: fake;

    typedef GeometryGridEntityPointer< BaseTraits, fake > Base;
    typedef GeometryGridEntityPointer< BaseTraits, fake > base;

  protected:
    typedef typename Traits :: HostEntityPointer HostEntityPointer;
    typedef typename Traits :: HostElementPointer HostElementPointer;
    typedef typename Traits :: HostElementIterator HostElementIterator;
    typedef typename Traits :: HostElement HostElement;

    typedef MakeableInterfaceObject< Entity > MakeableEntity;
    typedef typename MakeableEntity :: ImplementationType EntityImpl;

    HostElementIterator hostElementIterator_;
    int subEntity_;
    mutable MakeableEntity virtualEntity_;

  public:
    GeometryGridEntityPointer ( const Grid &grid,
                                const HostElementIterator &hostElementIterator,
                                int subEntity )
      : hostElementIterator_( hostElementIterator ),
        subEntity_( subEntity ),
        virtualEntity_( EntityImpl( grid ) )
    {}

    GeometryGridEntityPointer ( const Grid &grid,
                                const HostElement &hostElement,
                                int subEntity )
      : hostElementIterator_( hostElement.template entity< 0 >( 0 ) ),
        subEntity_( subEntity ),
        virtualEntity_( EntityImpl( grid ) )
    {}

    GeometryGridEntityPointer ( const This &other )
      : hostElementIterator_( other.hostElementIterator_ ),
        subEntity_( other.subEntity_ ),
        virtualEntity_( other.grid() )
    {}

    template< class T >
    GeometryGridEntityPointer ( const GeometryGridEntityPointer< T, fake > &other )
      : hostElementIterator_( other.hostElementIterator_ ),
        subEntity_( other.subEntity_ ),
        virtualEntity_( other.grid() )
    {}

    This &operator= ( const This &other )
    {
      hostElementIterator_ = other.hostElementIterator_;
      subEntity_ = other.subEntity_;
      update();
      return *this;
    }

    template< class T >
    bool equals ( const GeometryGridEntityPointer< T, fake > &other ) const
    {
      const int thisSub = subEntity_;
      const int otherSub = other.subEntity_;

      if( (thisSub < 0) || (otherSub < 0) )
        return (thisSub * otherSub >= 0);

      const int lvl = level();
      if( lvl != other.level() )
        return false;

      const typename Traits :: HostGrid :: Traits :: LevelIndexSet &indexSet
        = grid().hostGrid().levelIndexSet( lvl );

      const HostElement &thisElement = *hostElementPointer();
      assert( indexSet.contains( thisElement ) );
      const HostElement &otherElement = *(other.hostElementPointer());
      assert( indexSet.contains( otherElement ) );

      const int thisIndex
        = indexSet.template subIndex< codimension >( thisElement, thisSub );
      const int otherIndex
        = indexSet.template subIndex< codimension >( otherElement, otherSub );
      return thisIndex == otherIndex;
    }

    Entity &dereference () const
    {
      EntityImpl &impl = Grid :: getRealImplementation( virtualEntity_ );
      if( !impl.isValid() )
        impl.setToTarget( *hostElementPointer(), subEntity_ );
      return virtualEntity_;
    }

    int level () const
    {
      return hostElementPointer().level();
    }

    const HostEntityPointer &hostEntityPointer () const
    {
      DUNE_THROW( NotImplemented, "HostGrid has no entities of codimension "
                  << codimension << "." );
    }

  protected:
    const Grid &grid () const
    {
      return Grid :: getRealImplementation( virtualEntity_ ).grid();
    }

    void update ()
    {
      Grid :: getRealImplementation( virtualEntity_ ).invalidate();
    }

    const HostElementPointer &hostElementPointer () const
    {
      return hostElementIterator_;
    }
  };

}

#endif
