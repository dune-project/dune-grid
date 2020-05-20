#ifndef DUNE_GRID_CONCEPTS_INTERSECTION_HH
#define DUNE_GRID_CONCEPTS_INTERSECTION_HH

#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/geometry.hh>

#include <dune/common/concept.hh>

namespace Dune {
  namespace Concept
  {

    struct Intersection
    {
      template<class I>
      auto require(I&& i) -> decltype(
        requireConcept<Dune::Concept::EntityGeneral,typename I::Entity>(),
        // requireTrue<I::Entity::codimension == 1>(),
        requireConcept<Dune::Concept::Geometry,typename I::Geometry>(),
        requireConcept<Dune::Concept::Geometry,typename I::LocalGeometry>(),
        requireType<typename I::LocalCoordinate>(),
        requireType<typename I::GlobalCoordinate>(),
        requireConvertible<int>(I::mydimension),
        requireConvertible<int>(I::dimensionworld),
        requireType<typename I::ctype>(),
        requireConvertible<bool>(i.boundary()),
        requireConvertible<size_t>(i.boundarySegmentIndex()),
        requireConvertible<bool>(i.neighbor()),
        requireConvertible<typename I::Entity>(i.inside()),
        requireConvertible<typename I::Entity>(i.outside()),
        requireConvertible<bool>(i.conforming()),
        requireConvertible<typename I::LocalGeometry>(i.geometryInInside()),
        requireConvertible<typename I::LocalGeometry>(i.geometryInOutside()),
        requireConvertible<typename I::Geometry>(i.geometry()),
        requireConvertible<Dune::GeometryType>(i.type()),
        requireConvertible<int>(i.indexInInside()),
        requireConvertible<int>(i.indexInOutside()),
        requireConvertible<typename I::GlobalCoordinate>(i.outerNormal(/*local*/ typename I::LocalCoordinate{} )),
        requireConvertible<typename I::GlobalCoordinate>(i.integrationOuterNormal(/*local*/ typename I::LocalCoordinate{} )),
        requireConvertible<typename I::GlobalCoordinate>(i.unitOuterNormal(/*local*/ typename I::LocalCoordinate{} )),
        requireConvertible<typename I::GlobalCoordinate>(i.centerUnitOuterNormal()),
        i == i,
        i != i,
        I{},
        I{i},
        I{std::move(i)},
        i = i,
        i = std::move(i)
      );
    };

  }

  template <class I>
  constexpr bool isIntersection()
  {
    return models<Concept::Intersection, I>();
  }

}  // end namespace Dune

#endif
