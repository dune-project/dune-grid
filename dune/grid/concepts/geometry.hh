#ifndef DUNE_GRID_CONCEPTS_GEOMETRY_HH
#define DUNE_GRID_CONCEPTS_GEOMETRY_HH

#include <dune/geometry/type.hh>

#include <dune/common/concept.hh>

#if DUNE_HAVE_CXX_CONCEPTS
#include <dune/common/std/concepts.hh>
#endif

namespace Dune {
  namespace Concept{

#if DUNE_HAVE_CXX_CONCEPTS
namespace Concept{

    template<class G>
    concept Geometry = requires(G g)
    {
      not Std::default_initializable<G>;
      { G::coorddimension                                                       } -> Std::convertible_to<int>;
      { G::coorddimension                                                       } -> Std::convertible_to<int>;
      { g.type()                                                                } -> Std::convertible_to<Dune::GeometryType>;
      { g.affine()                                                              } -> Std::convertible_to<bool>;
      { g.corner(/*i*/ int{})                                                   } -> Std::convertible_to<typename G::GlobalCoordinate>;
      { g.global(/*local*/ typename G::LocalCoordinate{})                       } -> Std::convertible_to<typename G::GlobalCoordinate>;
      { g.local(/*global*/ typename G::GlobalCoordinate{})                      } -> Std::convertible_to<typename G::LocalCoordinate>;
      { g.volume()                                                              } -> Std::convertible_to<typename G::Volume>;
      { g.center()                                                              } -> Std::convertible_to<typename G::GlobalCoordinate>;
      { g.jacobianTransposed(/* local */ typename G::LocalCoordinate{})         } -> Std::convertible_to<typename G::JacobianTransposed>;
      { g.jacobianInverseTransposed(/* local */ typename G::LocalCoordinate{})  } -> Std::convertible_to<typename G::JacobianInverseTransposed>;
    };
}
#endif

    // namespace Fallback {
      struct Geometry
      {
        template<class G>
        auto require(G&& g) -> decltype(
          requireConvertible<int                                        >( G::coorddimension                                                       ),
          requireConvertible<int                                        >( G::mydimension                                                          ),
          requireConvertible<Dune::GeometryType                         >( g.type()                                                                ),
          requireConvertible<bool                                       >( g.affine()                                                              ),
          requireConvertible<int                                        >( g.corners()                                                             ),
          requireConvertible<typename G::GlobalCoordinate               >( g.corner(/* i */ int{})                                                 ),
          requireConvertible<typename G::GlobalCoordinate               >( g.global(/* local */ typename G::LocalCoordinate{})                     ),
          requireConvertible<typename G::LocalCoordinate                >( g.local(/* global */ typename G::GlobalCoordinate{})                    ),
          requireConvertible<typename G::ctype                          >( g.integrationElement(/* local */ typename G::LocalCoordinate{})         ),
          requireConvertible<typename G::Volume                         >( g.volume()                                                              ),
          requireConvertible<typename G::GlobalCoordinate               >( g.center()                                                              ),
          requireConvertible<typename G::JacobianTransposed             >( g.jacobianTransposed(/* local */ typename G::LocalCoordinate{})         ),
          requireConvertible<typename G::JacobianInverseTransposed      >( g.jacobianInverseTransposed(/* local */ typename G::LocalCoordinate{})  ),
          requireTrue<not std::is_default_constructible<G>::value>()
        );
      };
    }
  // }

  template <class G>
  constexpr void expectGeometry()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::Concept::Geometry<G>);
#else
    static_assert(models<Concept::Geometry, G>());
#endif
  }

}  // end namespace Dune

#endif
