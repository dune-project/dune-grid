#ifndef DUNE_GRID_CONCEPTS_ENTITY_HH
#define DUNE_GRID_CONCEPTS_ENTITY_HH

#include <dune/grid/concepts/anytype.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/geometry/type.hh>

#include <dune/common/concept.hh>

namespace Dune {
  namespace Concept
  {

    struct EntitySeed
    {
      template<class S>
      auto require(S&& seed) -> decltype(
        requireConvertible<int>(S::codimension),
        requireConvertible<bool>(seed.isValid()),
        S{}
      );
    };

    struct CodimensionEnum
    {
      template<class E>
      auto require(E&& e) -> decltype(
        requireConvertible<int>(E::codimension)
      );
    };

    struct EntityGeneral  : public Refines<CodimensionEnum>
    {
      template<class E>
      auto require(E&& e) -> decltype(
        requireConcept<Dune::Concept::Geometry,typename E::Geometry>(),
        requireConcept<Dune::Concept::EntitySeed,typename E::EntitySeed>(),
        requireConvertible<int>(E::dimension),
        requireConvertible<int>(E::mydimension),
        requireTrue<(int)E::mydimension==((int)E::dimension-(int)E::codimension)>(),
        requireConvertible<int>(e.level()),
        requireConvertible<Dune::PartitionType>(e.partitionType()),
        requireConvertible<typename E::Geometry>(e.geometry()),
        requireConvertible<Dune::GeometryType>(e.type()),
        requireConvertible<unsigned int>(e.subEntities(/*codim*/ (unsigned int){})),
        requireConvertible<typename E::EntitySeed>(e.seed()),
        requireConvertible<bool>(e==e),
        requireConvertible<bool>(e!=e),
        E{},              // default constructible
        E{e},             // copy constructible
        E{std::move(e)},  // move constructible
        e = e,            // copy assignable
        e = std::move(e)  // move assignable
      );
    };

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

    struct EntityExtended : public Refines<Dune::Concept::EntityGeneral>
    {
      template<class E>
      auto require(E&& e) -> decltype(
        requireTrue<E::codimension == 0>(),
        requireConcept<Dune::Concept::Geometry,typename E::LocalGeometry>(),
        requireType<typename E::HierarchicIterator>(),
        requireConvertible<E>(e.father()),
        requireConvertible<bool>(e.hasFather()),
        requireConvertible<bool>(e.isLeaf()),
        requireConvertible<bool>(e.isRegular()),
        requireConvertible<typename E::LocalGeometry>(e.geometryInFather()),
        requireConvertible<typename E::HierarchicIterator>(e.hbegin(/*maxLevel*/ int{})),
        requireConvertible<typename E::HierarchicIterator>(e.hend(/*maxLevel*/ int{})),
        requireConvertible<bool>(e.isNew()),
        requireConvertible<bool>(e.mightVanish()),
        requireConvertible<bool>(e.hasBoundaryIntersections()),
        requireConcept<EntityCodimExtended<E::dimension>,E>(),
        requireTrue<std::is_same<E,typename E::template Codim<0>::Entity>::value>()
      );
    };

    template<int codim>
    struct Entity : public Dune::Concept::EntityGeneral
    {};

    template<>
    struct Entity<0> : public Dune::Concept::EntityExtended
    {};

  }

  template <class S>
  constexpr bool isEntitySeed()
  {
    return models<Concept::EntitySeed, S>();
  }

  template <class E>
  constexpr bool isEntity()
  {
    if constexpr (models<Concept::CodimensionEnum, E>())
      return models<Concept::Entity<E::codimension>, E>();
    else
      return false;
  }

}  // end namespace Dune

#endif
