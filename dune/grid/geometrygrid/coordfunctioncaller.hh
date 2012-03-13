// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_COORDFUNCTIONCALLER_HH
#define DUNE_GEOGRID_COORDFUNCTIONCALLER_HH

#include <dune/grid/geometrygrid/hostcorners.hh>
#include <dune/grid/geometrygrid/coordfunction.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // CoordFunctionCaller
    // -------------------

    template< class HostEntity, class CoordFunctionInterface >
    class CoordFunctionCaller;

    template< class HostEntity, class ct, unsigned int dimD, unsigned int dimR, class Impl >
    class CoordFunctionCaller< HostEntity, AnalyticalCoordFunctionInterface< ct, dimD, dimR, Impl > >
    {
      typedef AnalyticalCoordFunctionInterface< ct, dimD, dimR, Impl > CoordFunctionInterface;
      typedef CoordFunctionCaller< HostEntity, CoordFunctionInterface > This;

      static const int codimension = HostEntity::codimension;

    public:
      typedef typename CoordFunctionInterface::RangeVector RangeVector;

      CoordFunctionCaller ( const HostEntity &hostEntity,
                            const CoordFunctionInterface &coordFunction )
        : hostCorners_( hostEntity ),
          coordFunction_( coordFunction )
      {}

      void evaluate ( unsigned int i, RangeVector &y ) const
      {
        coordFunction_.evaluate( hostCorners_.corner( i ), y );
      }

      GeometryType type () const
      {
        return hostCorners_.type();
      }

      unsigned int numCorners () const
      {
        return hostCorners_.numCorners();
      }

    private:
      const HostCorners< HostEntity > hostCorners_;
      const CoordFunctionInterface &coordFunction_;
    };

    template< class HostEntity, class ct, unsigned int dimR, class Impl >
    class CoordFunctionCaller< HostEntity, DiscreteCoordFunctionInterface< ct, dimR, Impl > >
    {
      typedef DiscreteCoordFunctionInterface< ct, dimR, Impl > CoordFunctionInterface;
      typedef CoordFunctionCaller< HostEntity, CoordFunctionInterface > This;

      typedef typename CoordFunctionInterface::RangeVector RangeVector;

    public:
      CoordFunctionCaller ( const HostEntity &hostEntity,
                            const CoordFunctionInterface &coordFunction )
        : hostEntity_( hostEntity ),
          coordFunction_( coordFunction )
      {}

      void evaluate ( unsigned int i, RangeVector &y ) const
      {
        coordFunction_.evaluate( hostEntity_, i, y );
      }

      GeometryType type () const
      {
        return hostEntity_.type();
      }

      unsigned int numCorners () const
      {
        return hostEntity_.geometry().corners();
      }

    private:
      const HostEntity &hostEntity_;
      const CoordFunctionInterface &coordFunction_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_COORDFUNCTIONCALLER_HH
