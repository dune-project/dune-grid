// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_ENTITY_ITERATOR_HH
#define DUNE_GRID_CONCEPTS_ENTITY_ITERATOR_HH

#include <dune/grid/concepts/entity.hh>

#include <dune/common/concept.hh>

namespace Dune {
  namespace Concept
  {

    struct EntityIterator
    {
      template<class I>
      auto require(I&& i) -> decltype(
        requireConcept<Dune::Concept::Entity<I::Entity::codimension>, typename I::Entity>(),
        requireType<typename I::Reference>(),
        i++, // FIXME set type requirement
        ++i, // FIXME set type requirement
        requireConvertible<typename I::Reference>(*i),
        i.operator ->(), // FIXME set type requirement
        requireConvertible<bool>(i==i),
        requireConvertible<bool>(i!=i),
        I{} // default constructible
      );
    };

  }

  template <class I>
  constexpr bool isEntityIterator()
  {
    return models<Concept::EntityIterator, I>();
  }

}  // end namespace Dune

#endif
