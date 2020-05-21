#ifndef DUNE_GRID_CONCEPTS_ENTITY_HH
#define DUNE_GRID_CONCEPTS_ENTITY_HH

#include <dune/grid/concepts/anytype.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/geometry/type.hh>

#include <dune/common/concept.hh>

#if DUNE_HAVE_CXX_CONCEPTS
#include <dune/common/std/concepts.hh>
#endif

namespace Dune {
  namespace Concept {
#if DUNE_HAVE_CXX_CONCEPTS
namespace Concept {
    template<class S>
    concept EntitySeed = requires(S seed)
    {
      requires Std::default_initializable<S>;
      { S::codimension  } -> Std::convertible_to<int>;
      { seed.isValid()  } -> Std::convertible_to<bool>;
    };
}
#endif

    // namespace Fallback {
      struct EntitySeed
      {
        template<class S>
        auto require(S&& seed) -> decltype(
          requireConvertible<int  >( S::codimension ),
          requireConvertible<bool >( seed.isValid() ),
          S{}
        );
      };
    // }

#if DUNE_HAVE_CXX_CONCEPTS
namespace Concept {
    template<class E>
    concept EntityGeneral = requires(E e)
    {
      requires Geometry<typename E::Geometry>;
      requires EntitySeed<typename E::EntitySeed>;
      E::mydimension==(E::dimension-E::codimension);
      { e.level()                                   } -> Std::convertible_to<int>;
      { e.partitionType()                           } -> Std::convertible_to<Dune::PartitionType>;
      { e.geometry()                                } -> Std::convertible_to<typename E::Geometry>;
      { e.type()                                    } -> Std::convertible_to<Dune::GeometryType>;
      { e.subEntities(/*codim*/ (unsigned int){})   } -> Std::convertible_to<unsigned int>;
      { e.seed()                                    } -> Std::convertible_to<typename E::EntitySeed>;
      { e==e                                        } -> Std::convertible_to<bool>;
      { e!=e                                        } -> Std::convertible_to<bool>;
      requires Std::default_initializable<E>;
      requires Std::copy_constructible<E>;
      requires Std::move_constructible<E>;
      e = e;
      e = std::move(e);
    };
}
#endif

    struct EntityGeneral
    {
      template<class E>
      auto require(E&& e) -> decltype(
        requireConcept<Dune::Concept::Geometry,typename E::Geometry>(),
        requireConcept<Dune::Concept::EntitySeed,typename E::EntitySeed>(),
        requireTrue<(int)E::mydimension==((int)E::dimension-(int)E::codimension)>(),
        requireConvertible<int                      >( e.level()                                    ),
        requireConvertible<Dune::PartitionType      >( e.partitionType()                            ),
        requireConvertible<typename E::Geometry     >( e.geometry()                                 ),
        requireConvertible<Dune::GeometryType       >( e.type()                                     ),
        requireConvertible<unsigned int             >( e.subEntities(/*codim*/ (unsigned int){})    ),
        requireConvertible<typename E::EntitySeed   >( e.seed()                                     ),
        requireConvertible<bool                     >( e==e                                         ),
        requireConvertible<bool                     >( e!=e                                         ),
        E{},              // default constructible
        E{e},             // copy constructible
        E{std::move(e)},  // move constructible
        e = e,            // copy assignable
        e = std::move(e)  // move assignable
      );
    };

#if DUNE_HAVE_CXX_CONCEPTS
namespace Concept {
    template<class E, int codim>
    concept EntityCodimExtended = requires(E e)
    {
      { e.template subEntity<codim>(/*sub_entity*/ int{}) } -> EntityGeneral;
    };

    template<class E, int codim = E::dimension>
    struct is_entity_codim_extended : std::conjunction<std::bool_constant<EntityCodimExtended<E,codim>>,is_entity_codim_extended<E,codim-1>> {};

    // Stop recursion
    template<class E>
    struct is_entity_codim_extended<E,0> : std::bool_constant<EntityCodimExtended<E,0>> {};
}
#endif

    template<int codim>
    struct EntityCodimExtended : public Refines<EntityCodimExtended<codim-1>>
    {
      template<class E>
      auto require(E&& e) -> decltype(
        requireConcept<Dune::Concept::EntityGeneral,typename E::template Codim<codim>::Entity>(),
        requireConvertible<typename E::template Codim<codim>::Entity>(e.template subEntity<codim>(/*sub_entity*/ int{}))
      );
    };

    // stop recursion
    template<>
    struct EntityCodimExtended<-1> : public AnyType {};


#if DUNE_HAVE_CXX_CONCEPTS
namespace Concept {
    template<class E>
    concept EntityExtended = requires(E e)
    {
      requires (E::codimension == 0);
      requires EntityGeneral<E>;
      requires Geometry<typename E::LocalGeometry>;
      { e.father()                      } -> Std::convertible_to<E>;
      { e.hasFather()                   } -> Std::convertible_to<bool>;
      { e.isLeaf()                      } -> Std::convertible_to<bool>;
      { e.isRegular()                   } -> Std::convertible_to<bool>;
      { e.geometryInFather()            } -> Std::convertible_to<typename E::LocalGeometry>;
      { e.hbegin(/*maxLevel*/ int{})    } -> Std::convertible_to<typename E::HierarchicIterator>;
      { e.hend(/*maxLevel*/ int{})      } -> Std::convertible_to<typename E::HierarchicIterator>;
      { e.isNew()                       } -> Std::convertible_to<bool>;
      { e.mightVanish()                 } -> Std::convertible_to<bool>;
      { e.hasBoundaryIntersections()    } -> Std::convertible_to<bool>;
      requires is_entity_codim_extended<E>::value; // Start recursion on codim entities
      requires std::is_same<E,typename E::template Codim<0>::Entity>::value;
    };
}
#endif

    struct EntityExtended : public Refines<Dune::Concept::EntityGeneral>
    {
      template<class E>
      auto require(E&& e) -> decltype(
        requireTrue<E::codimension == 0>(),
        requireConcept<Dune::Concept::Geometry,typename E::LocalGeometry>(),
        requireConvertible< E                               >( e.father()                   ),
        requireConvertible< bool                            >( e.hasFather()                ),
        requireConvertible< bool                            >( e.isLeaf()                   ),
        requireConvertible< bool                            >( e.isRegular()                ),
        requireConvertible< typename E::LocalGeometry       >( e.geometryInFather()         ),
        requireConvertible< typename E::HierarchicIterator  >( e.hbegin(/*maxLevel*/ int{}) ),
        requireConvertible< typename E::HierarchicIterator  >( e.hend(/*maxLevel*/ int{})   ),
        requireConvertible< bool                            >( e.isNew()                    ),
        requireConvertible< bool                            >( e.mightVanish()              ),
        requireConvertible< bool                            >( e.hasBoundaryIntersections() ),
        requireConcept<EntityCodimExtended<E::dimension>,E>(), // Start recursion codim entities
        requireTrue<std::is_same<E,typename E::template Codim<0>::Entity>::value>()
      );
    };

#if DUNE_HAVE_CXX_CONCEPTS
namespace Concept
{
    template<class E>
    concept Entity = EntityExtended<E> || EntityGeneral<E>;
}
#endif

    struct Entity
    {
      template<class E>
      auto require(E&& e) -> decltype(
        requireConcept<std::conditional_t<E::codimension == 0,EntityExtended,EntityGeneral>,E>()
      );
    };
  }

  template <class S>
  constexpr void expectEntitySeed()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::Concept::EntitySeed<S>);
#else
    static_assert(models<Concept::EntitySeed, S>());
#endif
  }

  template <class E>
  constexpr void expectEntityGeneral()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::Concept::EntityGeneral<E>);
#else
    static_assert(models<Concept::EntityGeneral, E>());
#endif
  }

  template <class E>
  constexpr void expectEntityExtended()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::Concept::EntityExtended<E>);
#else
    static_assert(models<Concept::EntityExtended, E>());
#endif
  }

  template <class E>
  constexpr void expectEntity()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Dune::Concept::Concept::Entity<E>);
#else
    static_assert(models<Dune::Concept::Entity,E>());
#endif
  }
}  // end namespace Dune

#endif
