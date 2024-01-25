// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_INDEX_SET_HH
#define DUNE_GRID_CONCEPTS_INDEX_SET_HH

#include <concepts>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/concepts/container.hh>
#include <dune/common/concepts/hashable.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/concepts/entity.hh>

namespace Dune::Concept {
namespace Impl {

  template<class IS, int codim>
  concept IndexSetEntityCodim = Entity<typename IS::template Codim<codim>::Entity> &&
    requires(const IS is, int i, unsigned int cc, const typename IS::template Codim<codim>::Entity& entity)
  {
    { is.template index<codim>(entity)         } -> std::same_as<typename IS::IndexType>;
    { is.index(entity)                         } -> std::same_as<typename IS::IndexType>;
    { is.template subIndex<codim>(entity,i,cc) } -> std::same_as<typename IS::IndexType>;
    { is.subIndex(entity,i,cc)                 } -> std::same_as<typename IS::IndexType>;
    { is.contains(entity)                      } -> std::convertible_to<bool>;
  };

  template<class IS, std::size_t... c>
  void indexSetEntityAllCodims(std::integer_sequence<std::size_t,c...>)
    requires (IndexSetEntityCodim<IS,int(c)> &&...);

} // end namespace Impl

/**
 * @brief Model of an index set
 * @ingroup GridConcepts
 * @details Dune::Grid::LevelIndexSet and Dune::Grid::LeafIndexSet are templates
 * for this model
 */
template<class IS>
concept IndexSet = requires(const IS is, Dune::GeometryType type, int codim)
{
  { IS::dimension   } -> std::convertible_to<int>;

  requires RandomAccessContainer<typename IS::Types>;
  { is.types(codim) } -> std::same_as<typename IS::Types>;

  requires std::integral<typename IS::IndexType>;
  { is.size(type)   } -> std::convertible_to<typename IS::IndexType>;
  { is.size(codim)  } -> std::convertible_to<typename IS::IndexType>;
} &&
Impl::IndexSetEntityCodim<IS,0> &&
requires (index_constant<1> from, index_constant<IS::dimension+1> to) {
  Impl::indexSetEntityAllCodims<IS>(range(from, to).to_integer_sequence());
};

namespace Impl {

  template<class IS, int codim>
  concept IdSetEntityCodim = Entity<typename IS::template Codim<codim>::Entity> &&
    requires(const IS is, const typename IS::template Codim<codim>::Entity& entity)
  {
    { is.template id<codim>(entity) } -> std::same_as<typename IS::IdType>;
    { is.id(entity)                 } -> std::same_as<typename IS::IdType>;
  };

  template<class IS, std::size_t... c>
  void idSetEntityAllCodims(std::integer_sequence<std::size_t,c...>)
    requires (IdSetEntityCodim<IS,int(c)> &&...);

} // end namespace Impl

/**
 * @brief Model of an id set
 * @ingroup GridConcepts
 * @details Dune::Grid::GlobalIdSet and Dune::Grid::LocalIdSet are templates
 * for this model
 */
template<class IS>
concept IdSet = requires(const IS is, const typename IS::template Codim<0>::Entity& entity, int i, unsigned int cc)
{
  requires Hashable<typename IS::IdType>;
  requires std::totally_ordered<typename IS::IdType>;
  { is.subId(entity,i,cc) } -> std::same_as<typename IS::IdType>;
} &&
Impl::IdSetEntityCodim<IS,0> &&
requires (index_constant<1> from, index_constant<IS::dimension+1> to) {
  Impl::idSetEntityAllCodims<IS>(range(from, to).to_integer_sequence());
};

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_INDEX_SET_HH
