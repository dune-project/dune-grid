#ifndef DUNE_GRID_CONCEPTS_ENTITY_ITERATOR_HH
#define DUNE_GRID_CONCEPTS_ENTITY_ITERATOR_HH

#include <dune/grid/concepts/entity.hh>

namespace Dune::Concept {

/**
 * @brief Model of an entity iterator
 * @ingroup GridConcepts
 * @details Dune::EntityIterator is a template for this model
 */
template<class It>
concept EntityIterator = requires(It it)
{
  requires Entity<typename It::Entity>;
  requires std::forward_iterator<It>;
  requires std::default_initializable<It>;
};

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_ENTITY_ITERATOR_HH
