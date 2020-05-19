#ifndef DUNE_GRID_CONCEPTS_INTERSECTION_ITERATOR_HH
#define DUNE_GRID_CONCEPTS_INTERSECTION_ITERATOR_HH

#include <dune/grid/concepts/intersection.hh>

#include <dune/common/concept.hh>

namespace Dune {
  namespace Concept
  {

    struct IntersectionIterator
    {
      template<class I>
      auto require(I&& i) -> decltype(
        requireType<typename I::Intersection>(),
        requireConcept<Dune::Concept::Intersection, typename I::Intersection>(),
        i++, // FIXME set type requirement
        ++i, // FIXME set type requirement
        *i, // FIXME set type requirement
        i.operator ->(), // FIXME set type requirement
        requireConvertible<bool>(i==i),
        requireConvertible<bool>(i!=i),
        I{}, // default constructible
        I{i} // copy constructible
      );
    };

  }

  template <class I>
  constexpr bool isIntersectionIterator()
  {
    return models<Concept::IntersectionIterator, I>();
  }

}  // end namespace Dune

#endif
