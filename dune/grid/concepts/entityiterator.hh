// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_ENTITY_ITERATOR_HH
#define DUNE_GRID_CONCEPTS_ENTITY_ITERATOR_HH

#include <concepts>
#include <iterator>

#include <dune/grid/concepts/entity.hh>

namespace Dune::Concept {

/**
 * @brief Model of an entity iterator
 * @ingroup GridConcepts
 * @details Dune::EntityIterator is a template for this model
 */
template<class It>
concept EntityIterator =
  std::forward_iterator<It> &&
  std::default_initializable<It> &&
  Entity<typename It::Entity>;

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_ENTITY_ITERATOR_HH
