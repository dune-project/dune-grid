#ifndef DUNE_GRID_CONCEPTS_ANY_TYPE_HH
#define DUNE_GRID_CONCEPTS_ANY_TYPE_HH

#include <dune/common/concept.hh>

namespace Dune {
  namespace Concept {
    namespace Fallback {
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
    } // nampespace Fallback
  } // nampespace Concept
} // end namespace Dune

#endif