#ifndef DUNE_GRID_CONCEPTS_GRIDVIEW_HH
#define DUNE_GRID_CONCEPTS_GRIDVIEW_HH

#include <dune/grid/concepts/anytype.hh>
#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/entityiterator.hh>
#include <dune/grid/concepts/intersection.hh>
#include <dune/grid/concepts/intersectioniterator.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/concepts/indexset.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/common/concept.hh>

#if DUNE_HAVE_CXX_CONCEPTS
#include <dune/common/std/concepts.hh>
#endif

namespace Dune {
  namespace Concept {


#if DUNE_HAVE_CXX_CONCEPTS

    template<class GV, int codim>
    concept GridViewCodim = requires(GV gv)
    {
      requires EntityIterator<typename GV::template Codim<codim>::Iterator>;
      requires EntityIterator<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::Iterator>;
      requires EntityIterator<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::Iterator>;
      requires EntityIterator<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::Iterator>;
      requires EntityIterator<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::Iterator>;
      requires EntityIterator<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::Iterator>;
      { gv.template begin<codim>()                                                       } -> Std::convertible_to< typename GV::template Codim<codim>::Iterator                                                                              >;
      { gv.template begin<codim,Dune::PartitionIteratorType::Interior_Partition>()       } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator         >;
      { gv.template begin<codim,Dune::PartitionIteratorType::InteriorBorder_Partition>() } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::Iterator   >;
      { gv.template begin<codim,Dune::PartitionIteratorType::Overlap_Partition>()        } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::Iterator          >;
      { gv.template begin<codim,Dune::PartitionIteratorType::OverlapFront_Partition>()   } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::Iterator     >;
      { gv.template begin<codim,Dune::PartitionIteratorType::All_Partition>()            } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::Iterator              >;
      { gv.template begin<codim,Dune::PartitionIteratorType::Ghost_Partition>()          } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::Iterator            >;
      { gv.template end<codim>()                                                         } -> Std::convertible_to< typename GV::template Codim<codim>::Iterator                                                                              >;
      { gv.template end<codim,Dune::PartitionIteratorType::Interior_Partition>()         } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator         >;
      { gv.template end<codim,Dune::PartitionIteratorType::InteriorBorder_Partition>()   } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::Iterator   >;
      { gv.template end<codim,Dune::PartitionIteratorType::Overlap_Partition>()          } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::Iterator          >;
      { gv.template end<codim,Dune::PartitionIteratorType::OverlapFront_Partition>()     } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::Iterator     >;
      { gv.template end<codim,Dune::PartitionIteratorType::All_Partition>()              } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::Iterator              >;
      { gv.template end<codim,Dune::PartitionIteratorType::Ghost_Partition>()            } -> Std::convertible_to< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::Iterator            >;
    };

    template<class GV, int codim = GV::dimension>
    struct is_grid_view_codim : std::conjunction<std::bool_constant<GridViewCodim<GV,codim>>,is_grid_view_codim<GV,codim-1>> {};

    // Stop recursion
    template<class GV>
    struct is_grid_view_codim<GV,0> : std::bool_constant<GridViewCodim<GV,0>> {};

#endif
    namespace Fallback {

      template<int codim>
      struct GridViewCodim : public Refines<GridViewCodim<codim-1>>
      {
        template<class GV>
        auto require(GV&& gv) -> decltype(
          requireConcept<EntityIterator, typename GV::template Codim<codim>::Iterator>(),
          requireConcept<Entity,typename GV::template Codim<codim>::Entity>(),
          requireConcept<Geometry,typename GV::template Codim<codim>::Geometry>(),
          requireConcept<Geometry,typename GV::template Codim<codim>::LocalGeometry>(),
          requireConcept<EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator>(),
          requireConcept<EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::Iterator>(),
          requireConcept<EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::Iterator>(),
          requireConcept<EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::Iterator>(),
          requireConcept<EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::Iterator>(),
          requireConcept<EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::Iterator>(),
          requireConvertible< typename GV::template Codim<codim>::Iterator                                                                              >( gv.template begin<codim>()                                                       ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator         >( gv.template begin<codim,Dune::PartitionIteratorType::Interior_Partition>()       ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::Iterator   >( gv.template begin<codim,Dune::PartitionIteratorType::InteriorBorder_Partition>() ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::Iterator          >( gv.template begin<codim,Dune::PartitionIteratorType::Overlap_Partition>()        ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::Iterator     >( gv.template begin<codim,Dune::PartitionIteratorType::OverlapFront_Partition>()   ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::Iterator              >( gv.template begin<codim,Dune::PartitionIteratorType::All_Partition>()            ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::Iterator            >( gv.template begin<codim,Dune::PartitionIteratorType::Ghost_Partition>()          ),
          requireConvertible< typename GV::template Codim<codim>::Iterator                                                                              >( gv.template end<codim>()                                                         ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator         >( gv.template end<codim,Dune::PartitionIteratorType::Interior_Partition>()         ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::Iterator   >( gv.template end<codim,Dune::PartitionIteratorType::InteriorBorder_Partition>()   ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::Iterator          >( gv.template end<codim,Dune::PartitionIteratorType::Overlap_Partition>()          ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::Iterator     >( gv.template end<codim,Dune::PartitionIteratorType::OverlapFront_Partition>()     ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::Iterator              >( gv.template end<codim,Dune::PartitionIteratorType::All_Partition>()              ),
          requireConvertible< typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::Iterator            >( gv.template end<codim,Dune::PartitionIteratorType::Ghost_Partition>()            )
        );
      };

      // stop recursion
      template<>
      struct GridViewCodim<-1> : public AnyType {};

    } // nampespace Fallback

#if DUNE_HAVE_CXX_CONCEPTS

    template<class GV>
    concept GridView = requires(GV gv, int codim, Dune::GeometryType type, const typename GV::template Codim<0>::Entity& entity)
    {
      typename GV::Traits;
      typename GV::ctype;
      requires IndexSet<typename GV::IndexSet>;
      requires Intersection<typename GV::Intersection>;
      requires IntersectionIterator<typename GV::IntersectionIterator>;
      { GV::conforming        } -> Std::convertible_to< bool                                  >;
      { GV::dimension         } -> Std::convertible_to< int                                   >;
      { GV::dimensionworld    } -> Std::convertible_to< int                                   >;
      { gv.grid()             } -> Std::convertible_to< const typename GV::Grid&              >;
      { gv.indexSet()         } -> Std::convertible_to< const typename GV::IndexSet&          >;
      { gv.size(codim)        } -> Std::convertible_to< int                                   >;
      { gv.size(type)         } -> Std::convertible_to< int                                   >;
      { gv.ibegin(entity)     } -> Std::convertible_to< typename GV::IntersectionIterator     >;
      { gv.iend(entity)       } -> Std::convertible_to< typename GV::IntersectionIterator     >;
      { gv.comm()             } -> Std::convertible_to< typename GV::CollectiveCommunication  >;
      { gv.overlapSize(codim) } -> Std::convertible_to< int                                   >;
      { gv.ghostSize(codim)   } -> Std::convertible_to< int                                   >;
      requires GridViewCodim<GV,0>; // Force compiler to issue errors on codim 0
      requires is_grid_view_codim<GV>::value; // Start recursion on other codims
      // gv.communicate(std::declval<Handler&>(),std::declval<Dune::InterfaceType>(), std::declval<CommunicationDirection>()), // FIXME use a default handler to instantiate this function
      Std::copy_constructible<GV>;
      gv = gv;
    };

#endif
    namespace Fallback {

      struct GridView
      {
        template<class GV>
        auto require(GV&& gv) -> decltype(
          requireType<typename GV::Traits>(),
          requireType<typename GV::ctype>(),
          requireConcept<IndexSet,typename GV::IndexSet>(),
          requireConcept<Intersection,typename GV::Intersection>(),
          requireConcept<IntersectionIterator,typename GV::IntersectionIterator>(),
          requireConcept<GridViewCodim<GV::dimension>,GV>(),
          requireConvertible< bool                                  >( GV::conforming                                                                        ),
          requireConvertible< int                                   >( GV::dimension                                                                         ),
          requireConvertible< int                                   >( GV::dimensionworld                                                                    ),
          requireConvertible< const typename GV::Grid&              >( gv.grid()                                                                             ),
          requireConvertible< const typename GV::IndexSet&          >( gv.indexSet()                                                                         ),
          requireConvertible< int                                   >( gv.size(/*codim*/ int{} )                                                             ),
          requireConvertible< int                                   >( gv.size(/*type*/ Dune::GeometryType{} )                                               ),
          requireConvertible< typename GV::IntersectionIterator     >( gv.ibegin(/*entity*/ std::declval<const typename GV::template Codim<0>::Entity&>() )  ),
          requireConvertible< typename GV::IntersectionIterator     >( gv.iend(/*entity*/ std::declval<const typename GV::template Codim<0>::Entity&>() )    ),
          requireConvertible< typename GV::CollectiveCommunication  >( gv.comm()                                                                             ),
          requireConvertible< int                                   >( gv.overlapSize(/*codim*/ int{} )                                                      ),
          requireConvertible< int                                   >( gv.ghostSize(/*codim*/ int{} )                                                        ),
          // gv.communicate(std::declval<Handler&>(),std::declval<Dune::InterfaceType>(), std::declval<CommunicationDirection>()), // FIXME use a default handler to instantiate this function
          GV{gv},
          gv = gv
        );
      };
    } // nampespace Fallback
  } // nampespace Concept

  template <class GV>
  constexpr void expectGridView()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::GridView<GV>);
#else
    static_assert(models<Concept::Fallback::GridView, GV>());
#endif
  }

}  // end namespace Dune

#endif
