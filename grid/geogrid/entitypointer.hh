// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITYPOINTER_HH
#define DUNE_GEOGRID_ENTITYPOINTER_HH

#include <dune/grid/common/grid.hh>

#include <dune/grid/geogrid/capabilities.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int codim, class Grid,
      bool fake = !Capabilities :: hasHostEntity< Grid, codim > :: v >
  class GeometryGridEntityPointer;



  // GeometryGridEntityPointer (real)
  // --------------------------------

  template< int codim, class Grid >
  class GeometryGridEntityPointer< codim, Grid, false >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

  public:
    enum { dimension = Traits :: dimension };
    enum { codimension = codim };

    typedef typename Traits :: template Codim< codim > :: Entity Entity;

    typedef GeometryGridEntityPointer< codim, Grid > Base;
    typedef GeometryGridEntityPointer< codim, Grid > base;

    static const bool fake = false;

  protected:
    typedef typename Traits :: HostGrid HostGrid;

    typedef typename HostGrid :: template Codim< codimension > :: EntityPointer
    HostEntityPointer;
    typedef typename HostGrid :: template Codim< 0 > :: Entity HostElement;

    typedef MakeableInterfaceObject< Entity > MakeableEntity;
    typedef typename MakeableEntity :: ImplementationType EntityImpl;

    HostEntityPointer hostEntityPointer_;
    mutable MakeableEntity virtualEntity_;

  public:
    GeometryGridEntityPointer ( const Grid &grid,
                                const HostEntityPointer &hostEntityPointer )
      : hostEntityPointer_( hostEntityPointer ),
        virtualEntity_( EntityImpl( grid ) )
    {}

    GeometryGridEntityPointer ( const Grid &grid,
                                const HostElement &hostElement,
                                int subEntity )
      : hostEntityPointer_( hostElement.template entity< codimension >( subEntity ) ),
        virtualEntity_( EntityImpl( grid ) )
    {}

    bool equals ( const GeometryGridEntityPointer &other ) const
    {
      return (hostEntityPointer() == other.hostEntityPointer());
    }

    Entity &dereference () const
    {
      EntityImpl &impl = Grid :: getRealImplementation( virtualEntity_ );
      if( !impl.isValid() )
        impl.setToTarget( *hostEntityPointer_ );
      return virtualEntity_;
    }

    //! ask for level of entity
    int level () const
    {
      return hostEntityPointer_.level();
    }

    const HostEntityPointer &hostEntityPointer () const
    {
      return hostEntityPointer_;
    }

  protected:
    const Grid &grid () const
    {
      EntityImpl &impl = Grid :: getRealImplementation( virtualEntity_ );
      return impl.grid();
    }

    void setToTarget ( const HostEntityPointer &target )
    {
      hostEntityPointer_ = target;
      Grid :: getRealImplementation( virtualEntity_ ).invalidate();
    }
  };



  // GeometryGridEntityPointer (fake)
  // --------------------------------

  template< int codim, class Grid >
  class GeometryGridEntityPointer< codim, Grid, true >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

  public:
    enum { dimension = Traits :: dimension };
    enum { codimension = codim };

    typedef typename Traits :: template Codim< codim > :: Entity Entity;

    typedef GeometryGridEntityPointer< codim, Grid > Base;
    typedef GeometryGridEntityPointer< codim, Grid > base;

    static const bool fake = true;

  protected:
    typedef typename Traits :: HostGrid HostGrid;

    typedef typename HostGrid :: template Codim< codim > :: EntityPointer
    HostEntityPointer;
    typedef typename HostGrid :: template Codim< 0 > :: Entity HostElement;
    typedef typename HostGrid :: template Codim< 0 > :: EntityPointer
    HostElementPointer;

    typedef MakeableInterfaceObject< Entity > MakeableEntity;
    typedef typename MakeableEntity :: ImplementationType EntityImpl;

    HostElementPointer hostElementPointer_;
    int subEntity_;
    mutable MakeableEntity virtualEntity_;

  public:
    GeometryGridEntityPointer ( const Grid &grid,
                                const HostElementPointer &hostElementPointer,
                                int subEntity )
      : hostElementPointer_( hostElementPointer ),
        subEntity_( subEntity ),
        virtualEntity_( EntityImpl( grid ) )
    {}

    GeometryGridEntityPointer ( const Grid &grid,
                                const HostElement &hostElement,
                                int subEntity )
      : hostElementPointer_( hostElement.template entity< 0 >( 0 ) ),
        subEntity_( subEntity ),
        virtualEntity_( EntityImpl( grid ) )
    {}

    bool equals ( const GeometryGridEntityPointer &other ) const
    {
      const int thisSub = subEntity_;
      const int otherSub = other.subEntity_;

      if( (thisSub < 0) || (otherSub < 0) )
        return (thisSub * otherSub >= 0);

      const int level = hostElementPointer_.level();
      if( level != other.hostElementPointer_.level() )
        return false;

      const typename HostGrid :: Traits :: LevelIndexSet &indexSet
        = grid().hostGrid().levelIndexSet( level );

      const HostElement &thisElement = *hostElementPointer_;
      assert( indexSet.contains( thisElement ) );
      const HostElement &otherElement = *(other.hostElementPointer_);
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
        impl.setToTarget( *hostElementPointer_, subEntity_ );
      return virtualEntity_;
    }

    //! ask for level of entity
    int level () const
    {
      return hostElementPointer_.level();
    }

    const HostEntityPointer &hostEntityPointer () const
    {
      DUNE_THROW( NotImplemented, "HostGrid has no entities of codimension "
                  << codimension << "." );
    }

  protected:
    const Grid &grid () const
    {
      EntityImpl &impl = Grid :: getRealImplementation( virtualEntity_ );
      return impl.grid();
    }

    const HostElementPointer &hostElementPointer () const
    {
      return hostElementPointer_;
    }

    int subEntity () const
    {
      return subEntity_;
    }

    void setToTarget ( const HostEntityPointer &target )
    {
      DUNE_THROW( NotImplemented, "HostGrid has no entities of codimension "
                  << codimension << "." );
    }

    void setToTarget ( const HostElementPointer &target, int subEntity )
    {
      hostElementPointer_ = target;
      subEntity_ = subEntity;
      Grid :: getRealImplementation( virtualEntity_ ).invalidate();
    }
  };

} // end namespace Dune

#endif
