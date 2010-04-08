// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAXIMUM_HH
#define DUNE_GENERICGEOMETRY_MAXIMUM_HH

namespace Dune
{

  template< template< int > class Value, int first, int last >
  struct Maximum
  {
    template< int v1, int v2 >
    struct M { static const int v = (v1 > v2 ? v1 : v2); };

    static const int v = M< Value< first >::v, Maximum< Value, first+1, last >::v >::v;

  private:
    dune_static_assert( (first <= last), "Maximum: first > last" );
  };

  template< template< int > class Value, int last >
  struct Maximum< Value, last, last >
  {
    static const int v = Value< last >::v;
  };

}

#endif // #ifndef DUNE_GENERICGEOMETRY_MAXIMUM_HH
