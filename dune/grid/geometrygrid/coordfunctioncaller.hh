// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
        coordFunction_.evaluate( hostCorners_[ i ], y );
      }

      GeometryType type () const
      {
        return hostCorners_.type();
      }

      std::size_t size () const
      {
        return hostCorners_.size();
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

      std::size_t size () const
      {
        auto refElement = referenceElement< ct, HostEntity::mydimension >( type() );
        return refElement.size( HostEntity::mydimension );
      }

    private:
      const HostEntity &hostEntity_;
      const CoordFunctionInterface &coordFunction_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_COORDFUNCTIONCALLER_HH
