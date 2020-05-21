// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_ENTITY_ITERATOR_HH
#define DUNE_GRID_CONCEPTS_ENTITY_ITERATOR_HH

#include <dune/grid/concepts/entity.hh>

#include <dune/common/concept.hh>

#if DUNE_HAVE_CXX_CONCEPTS
#include <dune/common/std/concepts.hh>
#endif

namespace Dune {
  namespace Concept {

#if DUNE_HAVE_CXX_CONCEPTS

    template<class I>
    concept EntityIterator = requires(I i)
    {
      requires Entity<typename I::Entity>;
      i++; // FIXME set type requirement
      ++i; // FIXME set type requirement
      { *i } -> Std::convertible_to<typename I::Reference>;
      i.operator ->(); // FIXME set type requirement
      { i==i } -> Std::convertible_to<bool>;
      { i!=i } -> Std::convertible_to<bool>;
      requires Std::default_initializable<I>;
    };

#endif
    namespace Fallback {
      struct EntityIterator
      {
        template<class I>
        auto require(I&& i) -> decltype(
          requireConcept<Entity,typename I::Entity>(),
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
    } // nampespace Fallback
  } // nampespace Concept

  template <class I>
  constexpr void expectEntityIterator()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Dune::Concept::EntityIterator<I>);
#else
    static_assert(models<Concept::Fallback::EntityIterator, I>());
#endif
  }

} // end namespace Dune

#endif
