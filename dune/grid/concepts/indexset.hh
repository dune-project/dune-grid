#ifndef DUNE_GRID_CONCEPTS_INDEX_SET_HH
#define DUNE_GRID_CONCEPTS_INDEX_SET_HH

#include <dune/grid/concepts/entity.hh>

#include <dune/geometry/type.hh>

#include <dune/common/indices.hh>

#include <type_traits>

namespace Dune::Concept {

namespace Impl {

  template<class IS, int codim>
  concept IndexSetEntityCodim = requires(IS is, int i, unsigned int sub_codim, const typename IS::template Codim<codim>::Entity& entity)
  {
    requires Entity<typename IS::template Codim<codim>::Entity>;
    { is.template index<codim>(entity)                    } -> std::convertible_to<typename IS::IndexType   >;
    { is.index(entity)                                    } -> std::convertible_to<typename IS::IndexType   >;
    { is.template subIndex<codim>(entity, i, sub_codim )  } -> std::convertible_to<typename IS::IndexType   >;
    { is.subIndex(entity, i, sub_codim )                  } -> std::convertible_to<typename IS::IndexType   >;
    { is.contains(entity)                                 } -> std::convertible_to<bool                     >;
  };

  template<class IS, std::size_t dim>
  concept AllIndexSetEntityCodims = requires(std::make_index_sequence<dim+1> codims)
  {
    []<std::size_t... cc>(std::index_sequence<cc...>)
        requires (IndexSetEntityCodim<IS,cc> &&...) {} (codims);
  };

}

/**
 * @brief Model of an index set
 * @ingroup GridConcepts
 * @details Dune::Grid::LevelIndexSet and Dune::Grid::LeafIndexSet are templates
 * for this model
 */
template<class IS>
concept IndexSet = requires(IS is, Dune::GeometryType type, int sub_codim)
{
  { IS::dimension       } -> std::convertible_to< int                     >;
  { is.types(sub_codim) } -> std::convertible_to< typename IS::Types      >;
  { is.size(type)       } -> std::convertible_to< typename IS::IndexType  >;
  { is.size(sub_codim)  } -> std::convertible_to< typename IS::IndexType  >;
  requires std::is_integral<typename IS::IndexType>::value;
  requires (not std::is_copy_constructible<IS>::value);
  requires (not std::is_copy_assignable<IS>::value);
  typename IS::Types;
  requires Impl::AllIndexSetEntityCodims<IS, IS::dimension>;
};

namespace Impl {

  template<class IS, int codim>
  concept IdSetEntityCodim = requires(IS is, const typename IS::template Codim<codim>::Entity& entity)
  {
    requires Entity<typename IS::template Codim<codim>::Entity>;
    { is.template id<codim>(entity)                    } -> std::convertible_to<typename IS::IdType   >;
    { is.id(entity)                                    } -> std::convertible_to<typename IS::IdType   >;
  };

  template<class IS, std::size_t dim>
  concept AllIdSetEntityCodims = requires(std::make_index_sequence<dim+1> codims)
  {
    []<std::size_t... cc>(std::index_sequence<cc...>)
        requires (IdSetEntityCodim<IS,cc> &&...) {} (codims);
  };

}

/**
 * @brief Model of an id set
 * @ingroup GridConcepts
 * @details Dune::Grid::GlobalIdSet and Dune::Grid::LocalIdSet are templates
 * for this model
 */
template<class IS>
concept IdSet = requires(IS is, const typename IS::template Codim<0>::Entity& entity, int i, unsigned int codim)
{
  { is.subId(entity,i,codim) } -> std::convertible_to<typename IS::IdType>;
  requires Impl::AllIdSetEntityCodims<IS, IS::dimension>;
};

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_INDEX_SET_HH
