// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_IDENTITY_HH
#define DUNE_GEOGRID_IDENTITY_HH

#include <dune/grid/geometrygrid/coordfunction.hh>

namespace Dune
{

  template< class ctype, unsigned int dim >
  class IdenticalCoordFunction
    : public AnalyticalCoordFunction
      < ctype, dim, dim, IdenticalCoordFunction< ctype, dim > >
  {
    typedef IdenticalCoordFunction< ctype, dim > This;
    typedef AnalyticalCoordFunction< ctype, dim, dim, This > Base;

  public:
    typedef typename Base :: DomainVector DomainVector;
    typedef typename Base :: RangeVector RangeVector;

    template< typename... Args >
    IdenticalCoordFunction( Args&... )
    {}

    RangeVector operator()(const DomainVector& x) const
    {
      return x;
    }

  };

}

#endif
