// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAXIMUM_HH
#define DUNE_GENERICGEOMETRY_MAXIMUM_HH

#include <dune/common/forloop.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // StaticMaximum
    // -------------

    template< class A, class B >
    struct StaticMaximum
    {
      static const int v = (A::v > B::v ? A::v : B::v);
    };



    // Maximum
    // -------

    template< template< int > class Value, int first, int last >
    struct Maximum
      : public GenericForLoop< StaticMaximum, Value, first, last >
    {};

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_MAXIMUM_HH
