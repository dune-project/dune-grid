// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_TEST_FUNCTIONS_HH
#define DUNE_GEOGRID_TEST_FUNCTIONS_HH

#include <dune/grid/geometrygrid/coordfunction.hh>
#include <dune/grid/geometrygrid/identity.hh>

namespace Dune
{

  class Helix
    : public AnalyticalCoordFunction< double, 2, 3, Helix >
  {
    typedef Helix This;
    typedef AnalyticalCoordFunction< double, 2, 3, This > Base;

  public:
    typedef Base :: DomainVector DomainVector;
    typedef Base :: RangeVector RangeVector;

    void evaluate ( const DomainVector &x, RangeVector &y ) const
    {
      y[ 0 ] = (x[ 0 ] + 0.2) * cos( 2.0 * M_PI * x[ 1 ] );
      y[ 1 ] = (x[ 0 ] + 0.2) * sin( 2.0 * M_PI * x[ 1 ] );
      y[ 2 ] = x[ 1 ];
    }
  };


  class Circle
    : public AnalyticalCoordFunction< double, 2, 2, Circle >
  {
    typedef Circle This;
    typedef AnalyticalCoordFunction< double, 2, 2, This > Base;

  public:
    typedef Base :: DomainVector DomainVector;
    typedef Base :: RangeVector RangeVector;

    void evaluate ( const DomainVector &x, RangeVector &y ) const
    {
      y[ 0 ] = (x[ 0 ] + 0.2) * cos( 2.0 * M_PI * x[ 1 ] );
      y[ 1 ] = (x[ 0 ] + 0.2) * sin( 2.0 * M_PI * x[ 1 ] );
    }
  };


  class ThickHelix
    : public AnalyticalCoordFunction< double, 3, 3, ThickHelix >
  {
    typedef ThickHelix This;
    typedef AnalyticalCoordFunction< double, 3, 3, This > Base;

  public:
    typedef Base :: DomainVector DomainVector;
    typedef Base :: RangeVector RangeVector;

    void evaluate ( const DomainVector &x, RangeVector &y ) const
    {
      double angle = x[1]-3.;
      if (std::abs(angle)<1.) {
        y[0] = x[0]+0.2;
        y[1] = angle;
        y[2] = x[2];
      }
      else {
        double fac = 1.;
        if (angle>0) {
          y[1] = 1.;
          angle -= 1.;
          fac=3.;
        }
        else {
          y[1] = -1;
          angle += 1.;
          fac=5;
        }
        y[ 0 ] = (x[ 0 ] + 0.2) * cos( 2.0 * M_PI * angle );
        y[ 1 ] += (x[ 0 ] + 0.2) * sin( 2.0 * M_PI * angle );
        y[ 2 ] = x[ 2 ] + fac*std::abs(angle);
      }
    }
  };

}

#endif
