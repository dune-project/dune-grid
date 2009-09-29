// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITYPOINTER_HH
#define DUNE_GEOGRID_ENTITYPOINTER_HH

#include <dune/grid/common/grid.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGridFamily;

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  // Internal Forward Declarations
  // -----------------------------

  template< int codim, class Grid >
  class GeometryGridEntityPointer;



  // GeometryGridEntityPointer
  // -------------------------

  template< int codim, class HostGrid, class CoordFunction >
  class GeometryGridEntityPointer< codim, const GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef typename GeometryGridFamily< HostGrid, CoordFunction > :: Traits Traits;

    typedef typename Traits :: Grid Grid;

  public:
    enum { dimension = Traits :: dimension };
    enum { codimension = codim };

    typedef typename Grid :: template Codim< codim > :: Entity Entity;

    typedef GeometryGridEntityPointer< codim, const Grid > Base;
    typedef GeometryGridEntityPointer< codim, const Grid > base;

  protected:
    typedef typename HostGrid :: template Codim< codim > :: EntityPointer HostEntityPointer;

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
      return dereference().level();
    }

    const HostEntityPointer &hostEntityPointer () const
    {
      return hostEntityPointer_;
    }

  protected:
    void setToTarget ( const HostEntityPointer &target )
    {
      hostEntityPointer_ = target;
      Grid :: getRealImplementation( virtualEntity_ ).invalidate();
    }
  };

} // end namespace Dune

#endif
