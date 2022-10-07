// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_HOSTCORNERS_HH
#define DUNE_GEOGRID_HOSTCORNERS_HH

#include <dune/geometry/type.hh>

#include <dune/grid/common/entity.hh>

namespace Dune
{
  namespace GeoGrid
  {

    // HostCorners
    // -----------

    template< class HostEntity >
    class HostCorners
    {
      typedef typename HostEntity::Geometry HostGeometry;

    public:
      typedef typename HostGeometry::GlobalCoordinate Coordinate;

      explicit HostCorners ( const HostEntity &hostEntity )
      : hostGeometry_( hostEntity.geometry() )
      {}

      GeometryType type () const
      {
        return hostGeometry_.type();
      }

      Coordinate operator[] ( int i ) const
      {
        return hostGeometry_.corner( i );
      }

      std::size_t size () const
      {
        return hostGeometry_.corners();
      }

    private:
      HostGeometry hostGeometry_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_HOSTCORNERS_HH
