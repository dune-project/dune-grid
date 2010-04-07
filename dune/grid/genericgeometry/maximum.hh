// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAXIMUM_HH
#define DUNE_GENERICGEOMETRY_MAXIMUM_HH

namespace Dune
{

  template< unsigned int v1 = 0, unsigned int v2 = 0, unsigned int v3 = 0, unsigned int v4 = 0,
      unsigned int v5 = 0, unsigned int v6 = 0, unsigned int v7 = 0, unsigned int v8 = 0 >
  struct Maximum
  {
    template< unsigned int w1, unsigned int w2 >
    struct M { static const unsigned int v = (w1 > w2 ? w1 : w2); };

    static const int v = M< M< M< v1, v2 >::v, M< v3, v4 >::v >::v, M< M< v5, v6 >::v, M< v7, v8 >::v >::v >::v;
  };

}

#endif // #ifndef DUNE_GENERICGEOMETRY_MAXIMUM_HH
