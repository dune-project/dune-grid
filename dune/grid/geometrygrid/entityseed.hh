// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITYSEED_HH
#define DUNE_GEOGRID_ENTITYSEED_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/geometrygrid/capabilities.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< int codim, class Grid, bool fake = !(Capabilities::hasHostEntity< Grid, codim >::v) >
    class EntitySeed;



    // EntitySeed (real)
    // -----------------

    template< int codim, class Grd >
    class EntitySeed< codim, Grd, false >
    {
      typedef typename remove_const< Grd >::type::Traits Traits;

    public:
      static const int codimension = codim;
      static const int dimension = Traits::dimension;
      static const int mydimension = dimension - codimension;
      static const int dimensionworld = Traits::dimensionworld;

      static const bool fake = false;

      typedef typename Traits::Grid Grid;
      typedef typename Traits::template Codim< codim >::Entity Entity;

      typedef typename Traits::HostGrid HostGrid;
      typedef typename HostGrid::template Codim< codim >::EntitySeed HostEntitySeed;

      explicit EntitySeed ( const HostEntitySeed &hostEntitySeed )
        : hostEntitySeed_( hostEntitySeed )
      {}

      bool operator== ( const EntitySeed &other )
      {
        return hostEntitySeed() == other.hostEntitySeed();
      }

      bool operator!= ( const EntitySeed &other )
      {
        return hostEntitySeed() != other.hostEntitySeed();
      }

      const HostEntitySeed &hostEntitySeed () const { return hostEntitySeed_; }

    private:
      HostEntitySeed hostEntitySeed_;
    };



    // EntitySeed (fake)
    // -----------------

    template< int codim, class Grd >
    class EntitySeed< codim, Grd, true >
    {
      typedef typename remove_const< Grd >::type::Traits Traits;

    public:
      static const int codimension = codim;
      static const int dimension = Traits::dimension;
      static const int mydimension = dimension - codimension;
      static const int dimensionworld = Traits::dimensionworld;

      static const bool fake = true;

      typedef typename Traits::Grid Grid;
      typedef typename Traits::template Codim< codim >::Entity Entity;

      typedef typename Traits::HostGrid HostGrid;
      typedef typename HostGrid::template Codim< 0 >::EntitySeed HostElementSeed;

      explicit EntitySeed ( const Grid &grid, const HostElementSeed &hostElementSeed, unsigned int subEntity )
        : grid_( &grid ),
          hostElementSeed_( hostElementSeed ),
          subEntity_( subEntity )
      {}

      bool operator== ( const EntitySeed &other )
      {
        return id() == other.id();
      }

      bool operator!= ( const EntitySeed &other )
      {
        return id() != other.id();
      }

      const HostElementSeed &hostElementSeed () const { return hostElementSeed_; }
      unsigned int subEntity () const { return subEntity_; }

    private:
      // very slow; only required for comparison
      typename Traits::LocalIdSet::IdType id () const
      {
        const HostGrid &hostGrid = grid_->hostGrid();
        return hostGrid.localIdSet().subId( *hostGrid.entityPointer( hostElementSeed() ), subEntity(), codim );
      }

      const Grid *grid_; // only required for comparison
      HostElementSeed hostElementSeed_;
      unsigned int subEntity_;
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_ENTITYSEED_HH
