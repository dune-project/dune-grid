#ifndef DUNE_GRID_CONCEPTS_GEOMETRY_HH
#define DUNE_GRID_CONCEPTS_GEOMETRY_HH

#include <dune/geometry/type.hh>

#include <dune/common/concept.hh>

namespace Dune {
  namespace Concept
  {

    struct Geometry
    {
      template<class G>
      auto require(G&& g) -> decltype(
        requireConvertible<int>(G::coorddimension),
        requireConvertible<int>(G::mydimension),
        requireType<typename G::ctype>(),
        requireType<typename G::LocalCoordinate>(),
        requireType<typename G::GlobalCoordinate>(),
        requireType<typename G::Volume>(),
        requireType<typename G::JacobianInverseTransposed>(),
        requireType<typename G::JacobianTransposed>(),
        requireConvertible<Dune::GeometryType>(g.type()),
        requireConvertible<bool>(g.affine()),
        requireConvertible<int>(g.corners()),
        requireConvertible<typename G::GlobalCoordinate>(g.corner(/* i */ int{})),
        requireConvertible<typename G::GlobalCoordinate>(g.global(/* local */ typename G::LocalCoordinate{})),
        requireConvertible<typename G::LocalCoordinate>(g.local(/* global */ typename G::GlobalCoordinate{})),
        requireConvertible<typename G::ctype>(g.integrationElement(/* local */ typename G::LocalCoordinate{})),
        requireConvertible<typename G::Volume>(g.volume()),
        requireConvertible<typename G::GlobalCoordinate>(g.center()),
        requireConvertible<typename G::JacobianTransposed>(g.jacobianTransposed(/* local */ typename G::LocalCoordinate{})),
        requireConvertible<typename G::JacobianInverseTransposed>(g.jacobianInverseTransposed(/* local */ typename G::LocalCoordinate{})),
        requireTrue<not std::is_default_constructible<G>::value>()
      );
    };
  }

  template <class G>
  constexpr bool isGeometry()
  {
    return models<Concept::Geometry, G>();;
  }

}  // end namespace Dune

#endif
