#ifndef DUNE_GRID_CONCEPTS_GEOMETRY_HH
#define DUNE_GRID_CONCEPTS_GEOMETRY_HH

#include <dune/geometry/type.hh>

namespace Dune::Concept {

template<class R>
concept ReferenceElement = true;

/**
 * @brief Model of a geometry object
 * @ingroup GridConcepts
 * @details Dune::Geometry is a template for this model
 */
template<class G>
concept Geometry = requires(const G& g, typename G::GlobalCoordinate global, typename G::LocalCoordinate local)
{
  typename G::ctype;
  { G::mydimension                      } -> std::convertible_to<int>;
  { G::coorddimension                   } -> std::convertible_to<int>;
  { g.type()                            } -> std::convertible_to<Dune::GeometryType>;
  { g.affine()                          } -> std::convertible_to<bool>;
  { g.corners()                         } -> std::convertible_to<int>;
  { g.corner(/*i*/ int{})               } -> std::convertible_to<typename G::GlobalCoordinate>;
  { g.global(local)                     } -> std::convertible_to<typename G::GlobalCoordinate>;
  { g.local(global)                     } -> std::convertible_to<typename G::LocalCoordinate>;
  { g.integrationElement(local)         } -> std::convertible_to<typename G::Volume>;
  { g.volume()                          } -> std::convertible_to<typename G::Volume>;
  { g.center()                          } -> std::convertible_to<typename G::GlobalCoordinate>;
  { g.jacobian(local)                   } -> std::convertible_to<typename G::Jacobian>;
  { g.jacobianInverse(local)            } -> std::convertible_to<typename G::JacobianInverse>;
  { g.jacobianTransposed(local)         } -> std::convertible_to<typename G::JacobianTransposed>;
  { g.jacobianInverseTransposed(local)  } -> std::convertible_to<typename G::JacobianInverseTransposed>;
  { referenceElement(g)                 } -> ReferenceElement;
};

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_GEOMETRY_HH
