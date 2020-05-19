#ifndef DUNE_GRID_CONCEPTS_INTERSECTION_ITERATOR_HH
#define DUNE_GRID_CONCEPTS_INTERSECTION_ITERATOR_HH

#include <dune/grid/concepts/intersection.hh>

namespace Dune::Concept {

/**
 * @brief Model of an intersection iterator
 * @ingroup GridConcepts
 * @details Dune::IntersectionIterator is a template for this model
 */
template<class It>
concept IntersectionIterator = requires(It it)
{
  requires Intersection<typename It::Intersection>;
  requires std::forward_iterator<It>;
  requires std::default_initializable<It>;
};

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_INTERSECTION_ITERATOR_HH
