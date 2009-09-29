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

    mutable MakeableEntity virtualEntity_;

  public:
    GeometryGridEntityPointer ( const Grid &grid, const HostEntityPointer &hostEntity )
      : virtualEntity_( EntityImpl( &grid, hostEntity ) )
    {}

    bool equals ( const GeometryGridEntityPointer &other ) const
    {
      const HostEntityPointer &thisHostEntity
        = Grid :: template getHostEntity< codim >( virtualEntity_ );
      const HostEntityPointer &otherHostEntity
        = Grid :: template getHostEntity< codim >( other.virtualEntity_ );
      return (thisHostEntity == otherHostEntity);
    }

    Entity &dereference () const
    {
      return virtualEntity_;
    }

    //! ask for level of entity
    int level () const
    {
      return dereference().level();
    }

  protected:
    void setToTarget ( const HostEntityPointer &target )
    {
      Grid :: getRealImplementation( virtualEntity_ ).setToTarget( target );
    }
  };

} // end namespace Dune

#endif
