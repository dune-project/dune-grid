// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_ENTITYSEED_HH
#define DUNE_GEOGRID_ENTITYSEED_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/entityseed.hh>
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
      typedef typename std::remove_const< Grd >::type::Traits Traits;

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

      //! default construct an invalid entity seed
      EntitySeed ( )
      {}

      explicit EntitySeed ( const HostEntitySeed &hostEntitySeed )
        : hostEntitySeed_( hostEntitySeed )
      {}

      //! check whether the EntitySeed refers to a valid Entity
      bool isValid() const
      {
        return hostEntitySeed_.isValid();
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
      typedef typename std::remove_const< Grd >::type::Traits Traits;

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

      //! default construct an invalid entity seed
      EntitySeed ( )
      {}

      explicit EntitySeed ( const HostElementSeed &hostElementSeed, unsigned int subEntity )
        : hostElementSeed_( hostElementSeed ),
          subEntity_( subEntity )
      {}

      //! check whether the EntitySeed refers to a valid Entity
      bool isValid() const
      {
        return hostElementSeed_.isValid();
      }

      const HostElementSeed &hostElementSeed () const { return hostElementSeed_; }
      unsigned int subEntity () const { return subEntity_; }

    private:
      HostElementSeed hostElementSeed_;
      unsigned int subEntity_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_ENTITYSEED_HH
