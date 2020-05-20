#ifndef DUNE_GRID_CONCEPTS_INDEX_SET_HH
#define DUNE_GRID_CONCEPTS_INDEX_SET_HH

#include <dune/grid/concepts/anytype.hh>
#include <dune/grid/concepts/entity.hh>

#include <dune/geometry/type.hh>

#include <dune/common/concept.hh>

namespace Dune {
  namespace Concept
  {

    template<int codim>
    struct IndexSetEntityCodim : public Refines<IndexSetEntityCodim<codim-1>>
    {
      template<class IS>
      auto require(IS&& is) -> decltype(
        requireConcept<Dune::Concept::Entity,typename IS::template Codim<codim>::Entity>(),
        requireConvertible<typename IS::IndexType>(is.template index<codim>(/*entity*/ std::declval<const typename IS::template Codim<codim>::Entity&>())),
        requireConvertible<typename IS::IndexType>(is.index(/*entity*/ std::declval<const typename IS::template Codim<codim>::Entity&>())),
        requireConvertible<typename IS::IndexType>(is.template subIndex<codim>(/*entity*/ std::declval<const typename IS::template Codim<codim>::Entity&>(), /*i*/ int{}, /*codim*/ (unsigned int){} )),
        requireConvertible<typename IS::IndexType>(is.subIndex(/*entity*/ std::declval<const typename IS::template Codim<codim>::Entity&>(), /*i*/ int{}, /*codim*/ (unsigned int){} )),
        requireConvertible<bool>(is.contains(/*entity*/ std::declval<const typename IS::template Codim<codim>::Entity&>()))
      );
    };

    // stop recursion
    template<>
    struct IndexSetEntityCodim<-1> : public AnyType {};

    struct IndexSet
    {
      template<class IS>
      auto require(IS&& is) -> decltype(
        requireConvertible<int>(IS::dimension),
        requireType<typename IS::IndexType>(),
        requireTrue<std::is_integral<typename IS::IndexType>::value>(),
        requireType<typename IS::Types>(),
        requireConvertible<typename IS::Types>(is.types(/*codim*/ int{} )),
        requireConvertible<typename IS::IndexType>(is.size(/*geometry_type*/ Dune::GeometryType{} )),
        requireConvertible<typename IS::IndexType>(is.size(/*codim*/ int{} )),
        requireConcept<IndexSetEntityCodim<IS::dimension>,IS>(),
        requireTrue<not std::is_copy_constructible<IS>::value>(),
        requireTrue<not std::is_copy_assignable<IS>::value>()
      );
    };

    template<int codim>
    struct IdSetEntityCodim : public Refines<IdSetEntityCodim<codim-1>>
    {
      template<class IS>
      auto require(IS&& is) -> decltype(
        requireType<typename IS::Grid::template Codim<codim>::Entity>(),
        requireConcept<Dune::Concept::Entity,typename IS::Grid::template Codim<codim>::Entity>(),
        requireConvertible<typename IS::IdType>(is.template id<codim>(/*entity*/ std::declval<const typename IS::Grid::template Codim<codim>::Entity&>())),
        requireConvertible<typename IS::IdType>(is.id(/*entity*/ std::declval<const typename IS::Grid::template Codim<codim>::Entity&>()))
      );
    };

    // stop recursion
    template<>
    struct IdSetEntityCodim<-1> : public AnyType {};


    struct IdSet
    {
      template<class IS>
      auto require(IS&& is) -> decltype(
        requireType<typename IS::IdType>(),
        requireConcept<IdSetEntityCodim<IS::Grid::dimension>,IS>(),
        requireConvertible<typename IS::IdType>(is.subId(/*entity*/ std::declval<const typename IS::Grid::template Codim<0>::Entity&>(), /*i*/ int{}, /*codim*/ (unsigned int){} ))
      );
    };
  }

  template <class IS>
  constexpr void expectIndexSet()
  {
    static_assert(models<Concept::IndexSet, IS>());
  }

  template <class IS>
  constexpr void expectIdSet()
  {
    static_assert(models<Concept::IdSet, IS>());
  }

}  // end namespace Dune

#endif
