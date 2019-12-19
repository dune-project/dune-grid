// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_ANY_TYPE_HH
#define DUNE_GRID_CONCEPTS_ANY_TYPE_HH

#include <dune/common/concept.hh>

namespace Dune {
  namespace Concept
  {

    struct AnyType {
      template<class E>
      auto require(E&& e) -> decltype(
        requireTrue<AlwaysTrue<E>::value>()
      );
    };

    struct NoType {
      template<class E>
      auto require(E&& e) -> decltype(
        requireTrue<AlwaysFalse<E>::value>()
      );
    };

  }
}  // end namespace Dune

#endif