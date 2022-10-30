// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_INDEX_SET_HH
#define DUNE_GRID_CONCEPTS_INDEX_SET_HH

#include <concepts>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/concepts/entity.hh>

namespace Dune::Concept {
namespace Impl {

  template<class IS, int codim,
    class E = typename IS::template Codim<codim>::Entity>
  concept IndexSetEntityCodim = requires(const IS is, int i, unsigned int sub_codim, const E& entity)
  {
    requires Entity<E>;
    { is.template index<codim>(entity)                  } -> std::convertible_to<typename IS::IndexType>;
    { is.index(entity)                                  } -> std::convertible_to<typename IS::IndexType>;
    { is.template subIndex<codim>(entity, i, sub_codim) } -> std::convertible_to<typename IS::IndexType>;
    { is.subIndex(entity, i, sub_codim)                 } -> std::convertible_to<typename IS::IndexType>;
    { is.contains(entity)                               } -> std::convertible_to<bool>;
  };

  template<class IS, int dim>
  concept AllIndexSetEntityCodims = requires(std::make_integer_sequence<int,dim+1> codims)
  {
    []<int... cc>(std::integer_sequence<int,cc...>)
        requires (IndexSetEntityCodim<IS,cc> &&...) {} (codims);
  };

} // end namespace Impl

/**
 * @brief Model of an index set
 * @ingroup GridConcepts
 * @details Dune::Grid::LevelIndexSet and Dune::Grid::LeafIndexSet are templates
 * for this model
 */
template<class IS>
concept IndexSet = requires(const IS is, Dune::GeometryType type, int sub_codim)
{
  { IS::dimension       } -> std::convertible_to<int>;
  { is.types(sub_codim) } -> std::convertible_to<typename IS::Types>;
  { is.size(type)       } -> std::convertible_to<typename IS::IndexType>;
  { is.size(sub_codim)  } -> std::convertible_to<typename IS::IndexType>;
  requires std::integral<typename IS::IndexType>;
  typename IS::Types;
  requires Impl::AllIndexSetEntityCodims<IS, IS::dimension>;
};


namespace Impl {

  template<class IS, int codim,
    class E = typename IS::template Codim<codim>::Entity>
  concept IdSetEntityCodim = requires(const IS is, const E& entity)
  {
    requires Entity<E>;
    { is.template id<codim>(entity) } -> std::convertible_to<typename IS::IdType>;
    { is.id(entity)                 } -> std::convertible_to<typename IS::IdType>;
  };

  template<class IS, int dim>
  concept AllIdSetEntityCodims = requires(std::make_integer_sequence<int,dim+1> codims)
  {
    []<int... cc>(std::integer_sequence<int,cc...>)
        requires (IdSetEntityCodim<IS,cc> &&...) {} (codims);
  };

} // end namespace Impl

/**
 * @brief Model of an id set
 * @ingroup GridConcepts
 * @details Dune::Grid::GlobalIdSet and Dune::Grid::LocalIdSet are templates
 * for this model
 */
template<class IS>
concept IdSet = requires(const IS is, const typename IS::template Codim<0>::Entity& entity, int i, unsigned int codim)
{
  { is.subId(entity,i,codim) } -> std::convertible_to<typename IS::IdType>;
  requires Impl::AllIdSetEntityCodims<IS, IS::dimension>;
};

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_INDEX_SET_HH