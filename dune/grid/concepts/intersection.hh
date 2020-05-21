#ifndef DUNE_GRID_CONCEPTS_INTERSECTION_HH
#define DUNE_GRID_CONCEPTS_INTERSECTION_HH

#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/geometry.hh>

#include <dune/common/concept.hh>

#if DUNE_HAVE_CXX_CONCEPTS
#include <dune/common/std/concepts.hh>
#endif

namespace Dune {
  namespace Concept {

/*!@defgroup ConceptIntersection Intersection
 * @{
 *  @ingroup Concepts
 *  @par Description
 *    This concept models how an intersection object should look like at compilation time.
 *    Dune::Intersection is a template for this model.
 *  @snippet this intersection-concept
 *  @par Uses
 *    - @ref ConceptEntity
 *    - @ref ConceptGeometry
 * @}
 */


#if DUNE_HAVE_CXX_CONCEPTS

    //! [intersection-concept]
    template<class I>
    concept Intersection = requires(I i, typename I::LocalCoordinate local)
    {
      requires EntityGeneral<typename I::Entity>;
      requires Geometry<typename I::Geometry>;
      requires Geometry<typename I::LocalGeometry>;
      typename I::ctype;
      { I::mydimension                  } -> Std::convertible_to<int                            >;
      { I::dimensionworld               } -> Std::convertible_to<int                            >;
      { i.boundary()                    } -> Std::convertible_to<bool                           >;
      { i.boundarySegmentIndex()        } -> Std::convertible_to<size_t                         >;
      { i.neighbor()                    } -> Std::convertible_to<bool                           >;
      { i.inside()                      } -> Std::convertible_to<typename I::Entity             >;
      { i.outside()                     } -> Std::convertible_to<typename I::Entity             >;
      { i.conforming()                  } -> Std::convertible_to<bool                           >;
      { i.geometryInInside()            } -> Std::convertible_to<typename I::LocalGeometry      >;
      { i.geometryInOutside()           } -> Std::convertible_to<typename I::LocalGeometry      >;
      { i.geometry()                    } -> Std::convertible_to<typename I::Geometry           >;
      { i.type()                        } -> Std::convertible_to<Dune::GeometryType             >;
      { i.indexInInside()               } -> Std::convertible_to<int                            >;
      { i.indexInOutside()              } -> Std::convertible_to<int                            >;
      { i.outerNormal(local)            } -> Std::convertible_to<typename I::GlobalCoordinate   >;
      { i.integrationOuterNormal(local) } -> Std::convertible_to<typename I::GlobalCoordinate   >;
      { i.unitOuterNormal(local)        } -> Std::convertible_to<typename I::GlobalCoordinate   >;
      { i.centerUnitOuterNormal()       } -> Std::convertible_to<typename I::GlobalCoordinate   >;
      { i==i                            } -> Std::convertible_to<bool                           >;
      { i!=i                            } -> Std::convertible_to<bool                           >;
      requires Std::default_initializable<I>;
      requires Std::copy_constructible<I>;
      requires Std::move_constructible<I>;
      i = i;
      i = std::move(i);
    };
    //! [intersection-concept]

#endif

    namespace Fallback {
      struct Intersection
      {
        template<class I>
        auto require(I&& i) -> decltype(
          requireConcept<EntityGeneral,typename I::Entity>(),
          requireConcept<Geometry,typename I::Geometry>(),
          requireConcept<Geometry,typename I::LocalGeometry>(),
          requireType<typename I::ctype>(),
          requireConvertible<int                            >( I::mydimension                                                     ),
          requireConvertible<int                            >( I::dimensionworld                                                  ),
          requireConvertible<bool                           >( i.boundary()                                                       ),
          requireConvertible<size_t                         >( i.boundarySegmentIndex()                                           ),
          requireConvertible<bool                           >( i.neighbor()                                                       ),
          requireConvertible<typename I::Entity             >( i.inside()                                                         ),
          requireConvertible<typename I::Entity             >( i.outside()                                                        ),
          requireConvertible<bool                           >( i.conforming()                                                     ),
          requireConvertible<typename I::LocalGeometry      >( i.geometryInInside()                                               ),
          requireConvertible<typename I::LocalGeometry      >( i.geometryInOutside()                                              ),
          requireConvertible<typename I::Geometry           >( i.geometry()                                                       ),
          requireConvertible<Dune::GeometryType             >( i.type()                                                           ),
          requireConvertible<int                            >( i.indexInInside()                                                  ),
          requireConvertible<int                            >( i.indexInOutside()                                                 ),
          requireConvertible<typename I::GlobalCoordinate   >( i.outerNormal(/*local*/ typename I::LocalCoordinate{} )            ),
          requireConvertible<typename I::GlobalCoordinate   >( i.integrationOuterNormal(/*local*/ typename I::LocalCoordinate{} ) ),
          requireConvertible<typename I::GlobalCoordinate   >( i.unitOuterNormal(/*local*/ typename I::LocalCoordinate{} )        ),
          requireConvertible<typename I::GlobalCoordinate   >( i.centerUnitOuterNormal()                                          ),
          requireConvertible<bool                           >( i == i                                                             ),
          requireConvertible<bool                           >( i != i                                                             ),
          I{},
          I{i},
          I{std::move(i)},
          i = i,
          i = std::move(i)
        );
      };
    } // nampespace Fallback
  } // nampespace Concept


  //! @expectConcept{ConceptIntersection,I}
  template <class I>
  constexpr void expectIntersection()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::Intersection<I>);
#else
    static_assert(models<Concept::Fallback::Intersection, I>());
#endif
  }

} // end namespace Dune

#endif
