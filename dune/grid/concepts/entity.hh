#ifndef DUNE_GRID_CONCEPTS_ENTITY_HH
#define DUNE_GRID_CONCEPTS_ENTITY_HH

#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/geometry/type.hh>

#include <dune/common/concept.hh>

#if DUNE_HAVE_CXX_CONCEPTS
#include <dune/common/std/concepts.hh>
#endif

namespace Dune {
  namespace Concept {

/*!@defgroup ConceptEntitySeed General entity
 * @{
  * @ingroup Concepts
 *  @par Description
 *    This concept models how _any_ entity object should look like at compilation time.
 *    Dune::Entity is a template for this model.
 *  @snippet this entity-seed-concept
 * @}
 */

#if DUNE_HAVE_CXX_CONCEPTS

    //! [entity-seed-concept]
    template<class S>
    concept EntitySeed = requires(S seed)
    {
      requires Std::default_initializable<S>;
      { S::codimension  } -> Std::convertible_to<int>;
      { seed.isValid()  } -> Std::convertible_to<bool>;
    };
    //! [entity-seed-concept]

#endif

    namespace Fallback {
      struct EntitySeed
      {
        template<class S>
        auto require(S&& seed) -> decltype(
          requireConvertible<int  >( S::codimension ),
          requireConvertible<bool >( seed.isValid() ),
          S{}
        );
      };
    } // nampespace Fallback

/*!@defgroup ConceptEntityGeneral General entity
 * @{
  * @ingroup Concepts
 *  @par Description
 *    This concept models how _any_ entity object should look like at compilation time.
 *    Dune::Entity is a template for this model.
 *  @snippet this general-entity-concept
 * @}
 */

#if DUNE_HAVE_CXX_CONCEPTS

    //! [general-entity-concept]
    template<class E>
    concept EntityGeneral = requires(E e, unsigned int codim)
    {
      requires Geometry<typename E::Geometry>;
      requires EntitySeed<typename E::EntitySeed>;
      E::mydimension==(E::dimension-E::codimension);
      { e.level()             } -> Std::convertible_to<int>;
      { e.partitionType()     } -> Std::convertible_to<Dune::PartitionType>;
      { e.geometry()          } -> Std::convertible_to<typename E::Geometry>;
      { e.type()              } -> Std::convertible_to<Dune::GeometryType>;
      { e.subEntities(codim)  } -> Std::convertible_to<unsigned int>;
      { e.seed()              } -> Std::convertible_to<typename E::EntitySeed>;
      { e==e                  } -> Std::convertible_to<bool>;
      { e!=e                  } -> Std::convertible_to<bool>;
      requires Std::default_initializable<E>;
      requires Std::copy_constructible<E>;
      requires Std::move_constructible<E>;
      e = e;
      e = std::move(e);
    };
    //! [general-entity-concept]

#endif

    namespace Fallback {
      struct EntityGeneral
      {
        template<class E>
        auto require(E&& e) -> decltype(
          requireConcept<Geometry,typename E::Geometry>(),
          requireConcept<EntitySeed,typename E::EntitySeed>(),
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
    } // nampespace Fallback

#if DUNE_HAVE_CXX_CONCEPTS

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

#endif

    namespace Fallback {
      template<int codim>
      struct EntityCodimExtended : public Refines<EntityCodimExtended<codim-1>>
      {
        template<class E>
        auto require(E&& e) -> decltype(
          requireConcept<EntityGeneral,typename E::template Codim<codim>::Entity>(),
          requireConvertible<typename E::template Codim<codim>::Entity>(e.template subEntity<codim>(/*sub_entity*/ int{}))
        );
      };

      // stop recursion
      template<>
      struct EntityCodimExtended<-1> {
        template<class E>
        std::true_type require(E&& e);
      };

    } // nampespace Fallback


/*!@defgroup ConceptEntityExtended Extended entity
 * @{
 *  @ingroup Concepts
 *  @par Description
 *    This concept models how an entity object that lives in codimension 0 should
 *    look like at compilation time. The specialization of Dune::Entity with
 *    codimension 0 is a template for this model.
 *  @snippet this general-entity-concept
 *  @par Refines
 *    - @ref ConceptEntityGeneral
 *    - @ref ConceptGeometry
 * @}
 */
#if DUNE_HAVE_CXX_CONCEPTS

    //! [extended-entity-concept]
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
      requires EntityCodimExtended<E,0>; //! Force compiler to issue errors on codim 0
      requires is_entity_codim_extended<E>::value; //! Start recursion on codim entities
      requires std::is_same<E,typename E::template Codim<0>::Entity>::value;
    };
    //! [extended-entity-concept]

#endif

    namespace Fallback {

      struct EntityExtended : public Refines<EntityGeneral>
      {
        template<class E>
        auto require(E&& e) -> decltype(
          requireTrue<E::codimension == 0>(),
          requireConcept<Geometry,typename E::LocalGeometry>(),
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
    } // nampespace Fallback


/*!@defgroup ConceptEntity Entity
 * @{
 *  @ingroup Concepts
 *  @par Description
 *    This concept models how an entity object should look like at compilation time.
 *    Dune::Entity is a template for this model.
 *  @snippet this entity-concept
 *  @par Refines
 *    - @ref ConceptEntityGeneral
 *  @par Uses
 *    - @ref ConceptEntityExtended
 * @}
 */

#if DUNE_HAVE_CXX_CONCEPTS

    //! [entity-concept]
    template<class E>
    concept Entity = EntityExtended<E> || EntityGeneral<E>;
    //! [entity-concept]
#endif

    namespace Fallback {

      struct Entity
      {
        template<class E>
        auto require(E&& e) -> decltype(
          requireConcept<std::conditional_t<E::codimension == 0,EntityExtended,EntityGeneral>,E>()
        );
      };
    } // nampespace Fallback
  } // nampespace Concept


  //! @expectConcept{ConceptEntitySeed,S}
  template <class S>
  constexpr void expectEntitySeed()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::EntitySeed<S>);
#else
    static_assert(models<Concept::Fallback::EntitySeed, S>());
#endif
  }

  //! @expectConcept{ConceptEntityGeneral,E}
  template <class E>
  constexpr void expectEntityGeneral()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::EntityGeneral<E>);
#else
    static_assert(models<Concept::Fallback::EntityGeneral, E>());
#endif
  }

  //! @expectConcept{ConceptEntityExtended,E}
  template <class E>
  constexpr void expectEntityExtended()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::EntityExtended<E>);
#else
    static_assert(models<Concept::Fallback::EntityExtended, E>());
#endif
  }

  //! @expectConcept{ConceptEntity,E}
  //! @ingroup{ConceptEntity}
  template <class E>
  constexpr void expectEntity()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::Entity<E>);
#else
    static_assert(models<Concept::Fallback::Entity,E>());
#endif
  }
} // end namespace Dune

#endif
