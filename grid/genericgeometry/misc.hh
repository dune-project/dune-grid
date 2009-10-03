// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MISC_HH
#define DUNE_GENERICGEOMETRY_MISC_HH

#include <dune/common/static_assert.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< bool condition, template< bool > class True, template< bool > class False >
    struct ProtectedIf;

    template< template< bool > class True, template< bool > class False >
    struct ProtectedIf< true, True, False >
      : public True< true >
    {};

    template< template< bool > class True, template< bool > class False >
    struct ProtectedIf< false, True, False >
      : public False< false >
    {};

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_MISC_HH
