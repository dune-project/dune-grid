#ifndef DUNE_GRID_CONCEPTS_INDEX_SET_HH
#define DUNE_GRID_CONCEPTS_INDEX_SET_HH

#include <dune/grid/concepts/anytype.hh>
#include <dune/grid/concepts/entity.hh>

#include <dune/geometry/type.hh>

#include <dune/common/concept.hh>

#if DUNE_HAVE_CXX_CONCEPTS
#include <dune/common/std/concepts.hh>
#endif

namespace Dune {
  namespace Concept {

#if DUNE_HAVE_CXX_CONCEPTS

    template<class IS, int codim>
    concept IndexSetEntityCodim = requires(IS is, int i, unsigned int sub_codim, const typename IS::template Codim<codim>::Entity& entity)
    {
      requires Entity<typename IS::template Codim<codim>::Entity>;
      { is.template index<codim>(entity)                    } -> Std::convertible_to<typename IS::IndexType   >;
      { is.index(entity)                                    } -> Std::convertible_to<typename IS::IndexType   >;
      { is.template subIndex<codim>(entity, i, sub_codim )  } -> Std::convertible_to<typename IS::IndexType   >;
      { is.subIndex(entity, i, sub_codim )                  } -> Std::convertible_to<typename IS::IndexType   >;
      { is.contains(entity)                                 } -> Std::convertible_to<bool                     >;
    };

    template<class IS, int codim = IS::dimension>
    struct is_index_set_entity_codim : std::conjunction<std::bool_constant<IndexSetEntityCodim<IS,codim>>,is_index_set_entity_codim<IS,codim-1>> {};

    // Stop recursion
    template<class IS>
    struct is_index_set_entity_codim<IS,0> : std::bool_constant<IndexSetEntityCodim<IS,0>> {};

#endif

    namespace Fallback {
      template<int codim>
      struct IndexSetEntityCodim : public Refines<IndexSetEntityCodim<codim-1>>
      {
        template<class IS>
        auto require(IS&& is) -> decltype(
          requireConcept<Entity,typename IS::template Codim<codim>::Entity>(),
          requireConvertible<typename IS::IndexType   >( is.template index<codim>(/*entity*/ std::declval<const typename IS::template Codim<codim>::Entity&>())                                                 ),
          requireConvertible<typename IS::IndexType   >( is.index(/*entity*/ std::declval<const typename IS::template Codim<codim>::Entity&>())                                                                 ),
          requireConvertible<typename IS::IndexType   >( is.template subIndex<codim>(/*entity*/ std::declval<const typename IS::template Codim<codim>::Entity&>(), /*i*/ int{}, /*sub_codim*/ (unsigned int){} )    ),
          requireConvertible<typename IS::IndexType   >( is.subIndex(/*entity*/ std::declval<const typename IS::template Codim<codim>::Entity&>(), /*i*/ int{}, /*sub_codim*/ (unsigned int){} )                    ),
          requireConvertible<bool                     >( is.contains(/*entity*/ std::declval<const typename IS::template Codim<codim>::Entity&>())                                                              )
        );
      };

      // Stop recursion
      template<>
      struct IndexSetEntityCodim<-1> : public AnyType {};

    } // nampespace Fallback

#if DUNE_HAVE_CXX_CONCEPTS

    template<class IS>
    concept IndexSet = requires(IS is, Dune::GeometryType type, int sub_codim)
    {
      { IS::dimension       } -> Std::convertible_to< int                     >;
      { is.types(sub_codim) } -> Std::convertible_to< typename IS::Types      >;
      { is.size(type)       } -> Std::convertible_to< typename IS::IndexType  >;
      { is.size(sub_codim)  } -> Std::convertible_to< typename IS::IndexType  >;
      requires std::is_integral<typename IS::IndexType>::value;
      requires (not std::is_copy_constructible<IS>::value);
      requires (not std::is_copy_assignable<IS>::value);
      typename IS::Types;
      requires IndexSetEntityCodim<IS,0>; // Force compiler to issue errors on codim 0
      requires is_index_set_entity_codim<IS>::value; // Start recursion on codim entities
    };
#endif

    namespace Fallback {
      struct IndexSet
      {
        template<class IS>
        auto require(IS&& is) -> decltype(
        requireConvertible< int                     >( IS::dimension                                    ),
        requireConvertible< typename IS::Types      >( is.types(/*codim*/ int{} )                       ),
        requireConvertible< typename IS::IndexType  >( is.size(/*geometry_type*/ Dune::GeometryType{} ) ),
        requireConvertible< typename IS::IndexType  >( is.size(/*codim*/ int{} )                        ),
        requireTrue< std::is_integral<typename IS::IndexType>::value  >(),
        requireTrue< not std::is_copy_constructible<IS>::value        >(),
        requireTrue< not std::is_copy_assignable<IS>::value           >(),
        requireType<typename IS::Types>(),
        requireConcept<IndexSetEntityCodim<IS::dimension>,IS>() // Start recursion codim entities
        );
      };
    } // nampespace Fallback


#if DUNE_HAVE_CXX_CONCEPTS

    template<class IS, int codim>
    concept IdSetEntityCodim = requires(IS is, const typename IS::Grid::template Codim<codim>::Entity& entity)
    {
      requires Entity<typename IS::Grid::template Codim<codim>::Entity>;
      { is.template id<codim>(entity)                    } -> Std::convertible_to<typename IS::IdType   >;
      { is.id(entity)                                    } -> Std::convertible_to<typename IS::IdType   >;
    };

    template<class IS, int codim = IS::Grid::dimension>
    struct is_id_set_entity_codim : std::conjunction<std::bool_constant<IdSetEntityCodim<IS,codim>>,is_id_set_entity_codim<IS,codim-1>> {};

    // Stop recursion
    template<class IS>
    struct is_id_set_entity_codim<IS,0> : std::bool_constant<IdSetEntityCodim<IS,0>> {};

#endif

    namespace Fallback {

      template<int codim>
      struct IdSetEntityCodim : public Refines<IdSetEntityCodim<codim-1>>
      {
        template<class IS>
        auto require(IS&& is) -> decltype(
          requireConcept<Entity,typename IS::Grid::template Codim<codim>::Entity>(),
          requireConvertible<typename IS::IdType>(is.template id<codim>(/*entity*/ std::declval<const typename IS::Grid::template Codim<codim>::Entity&>())),
          requireConvertible<typename IS::IdType>(is.id(/*entity*/ std::declval<const typename IS::Grid::template Codim<codim>::Entity&>()))
        );
      };

      // stop recursion
      template<>
      struct IdSetEntityCodim<-1> : public AnyType {};

    }  // nampespace Fallback

#if DUNE_HAVE_CXX_CONCEPTS

    template<class IS>
    concept IdSet = requires(IS is, const typename IS::Grid::template Codim<0>::Entity& entity, int i, unsigned int codim)
    {
      { is.subId(entity,i,codim) } -> Std::convertible_to<typename IS::IdType>;
      requires IdSetEntityCodim<IS,0>; // Force compiler to issue errors on codim 0
      requires is_id_set_entity_codim<IS>::value; // Start recursion on codim entities
    };

#endif

    namespace Fallback {
      struct IdSet
      {
        template<class IS>
        auto require(IS&& is) -> decltype(
          requireConcept<IdSetEntityCodim<IS::Grid::dimension>,IS>(), // Start recursion codim entities
          requireConvertible<typename IS::IdType>(is.subId(/*entity*/ std::declval<const typename IS::Grid::template Codim<0>::Entity&>(), /*i*/ int{}, /*codim*/ (unsigned int){} ))
        );
      };
    } // nampespace Fallback
  } // nampespace Concept

  template <class IS>
  constexpr void expectIndexSet()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::IndexSet<IS>);
#else
    static_assert(models<Concept::Fallback::IndexSet, IS>());
#endif
  }

  template <class IS>
  constexpr void expectIdSet()
  {

#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::IdSet<IS>);
#else
    static_assert(models<Concept::Fallback::IdSet, IS>());
#endif
  }

} // end namespace Dune

#endif
