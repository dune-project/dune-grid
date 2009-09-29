// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_IDENTITY_HH
#define DUNE_GEOGRID_IDENTITY_HH

#include <dune/common/fvector.hh>

namespace Dune
{

  template< class ctype, unsigned int dim >
  struct IdenticalCoordFunction
  {
    enum { dimRange = dim };
    enum { dimDomain = dim };

    typedef Dune :: FieldVector< ctype, dim > Vector;

    void evaluate ( const Vector &x, Vector &y )
    {
      y = x;
    }
  };

}

#endif
