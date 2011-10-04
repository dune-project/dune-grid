// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_HOSTCORNERS_HH
#define DUNE_GEOGRID_HOSTCORNERS_HH

#include <dune/geometry/type.hh>

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

      HostCorners ( const HostEntity &hostEntity )
        : hostGeometry_( hostEntity.geometry() )
      {}

      GeometryType type () const
      {
        return hostGeometry_.type();
      }

      Coordinate corner ( const int i ) const
      {
        return hostGeometry_.corner( i );
      }

      unsigned int numCorners () const
      {
        return hostGeometry_.corners();
      }

    private:
      const HostGeometry &hostGeometry_;
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_HOSTCORNERS_HH
