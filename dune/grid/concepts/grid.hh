#ifndef DUNE_GRID_CONCEPTS_GRID_HH
#define DUNE_GRID_CONCEPTS_GRID_HH

#include <dune/grid/concepts/anytype.hh>
#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/entityiterator.hh>
#include <dune/grid/concepts/intersection.hh>
#include <dune/grid/concepts/intersectioniterator.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/concepts/indexset.hh>
#include <dune/grid/concepts/gridview.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/common/concept.hh>

#if DUNE_HAVE_CXX_CONCEPTS
#include <dune/common/std/concepts.hh>
#endif

namespace Dune {
  namespace Concept {

#if DUNE_HAVE_CXX_CONCEPTS

    template<class G, int codim>
    concept GridCodim = requires(G g, const typename G::template Codim<codim>::EntitySeed& seed)
    {
      requires Entity< typename G::template Codim<codim>::Entity>;
      requires EntitySeed< typename G::template Codim<codim>::EntitySeed>;
      requires Geometry< typename G::template Codim<codim>::Geometry>;
      requires Geometry< typename G::template Codim<codim>::LocalGeometry>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::LevelIterator>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::LeafIterator>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::LevelIterator>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::LeafIterator>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::LevelIterator>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::LeafIterator>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::LevelIterator>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::LeafIterator>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::LevelIterator>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::LeafIterator>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::LevelIterator>;
      requires EntityIterator< typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::LeafIterator>;
      requires EntityIterator< typename G::template Codim<codim>::LevelIterator>;
      requires EntityIterator< typename G::template Codim<codim>::LeafIterator>;
      { g.entity(seed) } -> Std::convertible_to<typename G::template Codim<codim>::Entity                                                                                >;
    };


    template<class G, int codim = G::dimension>
    struct is_grid_codim : std::conjunction<std::bool_constant<GridCodim<G,codim>>,is_grid_codim<G,codim-1>> {};

    // Stop recursion
    template<class G>
    struct is_grid_codim<G,0> : std::bool_constant<GridCodim<G,0>> {};

#endif

    namespace Fallback {

      template<int codim>
      struct GridCodim : public Refines<GridCodim<codim-1>>
      {
        template<class G>
        auto require(G&& g) -> decltype(
          requireConcept<Entity,typename G::template Codim<codim>::Entity>(),
          requireConcept<EntitySeed,typename G::template Codim<codim>::EntitySeed>(),
          requireConcept<Geometry,typename G::template Codim<codim>::Geometry>(),
          requireConcept<Geometry,typename G::template Codim<codim>::LocalGeometry>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::LevelIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::LeafIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::LevelIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::LeafIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::LevelIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::LeafIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::LevelIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::LeafIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::LevelIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::LeafIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::LevelIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::LeafIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::LevelIterator>(),
          requireConcept<EntityIterator, typename G::template Codim<codim>::LeafIterator>(),
          requireConvertible<typename G::template Codim<codim>::Entity>( g.entity(/*seed*/ std::declval<const typename G::template Codim<codim>::EntitySeed&>())  )
        );
      };

      // stop recursion
      template<>
      struct GridCodim<-1> : public AnyType {};

    } // nampespace Fallback

#if DUNE_HAVE_CXX_CONCEPTS

    template<class G>
    concept Grid = requires(G g, int level, int codim, int refCount, Dune::GeometryType type, const typename G::template Codim<0>::Entity& entity)
    {
      { G::dimension      } -> Std::convertible_to<int>;
      { G::dimensionworld } -> Std::convertible_to<int>;
      requires GridView< typename G::LeafGridView>;
      requires GridView< typename G::LevelGridView>;
      requires Intersection< typename G::LeafIntersection>;
      requires Intersection< typename G::LevelIntersection>;
      requires IntersectionIterator< typename G::LeafIntersectionIterator>;
      requires IntersectionIterator< typename G::LevelIntersectionIterator>;
      requires IndexSet< typename G::LevelIndexSet>;
      requires IndexSet< typename G::LeafIndexSet>;
      requires IdSet< typename G::GlobalIdSet>;
      requires IdSet< typename G::LocalIdSet>;
      requires std::is_same<G,typename G::LeafGridView::Grid>::value;
      requires std::is_same<G,typename G::LevelGridView::Grid>::value;
      requires GridCodim<G,0>; // Force compiler to show errors on codim 0 if interface assertion fails
      requires is_grid_codim<G>::value; // Start recursion on codims
      typename G::ctype;
      typename G::HierarchicIterator;
      { g.maxLevel()              } -> Std::convertible_to< int                                 >;
      { g.size(level, codim)      } -> Std::convertible_to< int                                 >;
      { g.size(codim)             } -> Std::convertible_to< int                                 >;
      { g.size(level, type)       } -> Std::convertible_to< int                                 >;
      { g.size(type)              } -> Std::convertible_to< int                                 >;
      { g.numBoundarySegments()   } -> Std::convertible_to< std::size_t                         >;
      { g.levelGridView(level)    } -> Std::convertible_to< typename G::LevelGridView           >;
      { g.leafGridView()          } -> Std::convertible_to< typename G::LeafGridView            >;
      { g.globalIdSet()           } -> Std::convertible_to< const typename G::GlobalIdSet&      >;
      { g.localIdSet()            } -> Std::convertible_to< const typename G::LocalIdSet&       >;
      { g.levelIndexSet(level)    } -> Std::convertible_to< const typename G::LevelIndexSet&    >;
      { g.leafIndexSet()          } -> Std::convertible_to< const typename G::LeafIndexSet&     >;
      { g.mark(refCount,entity)   } -> Std::convertible_to< bool                                >;
      { g.getMark(entity)         } -> Std::convertible_to< int                                 >;
      { g.preAdapt()              } -> Std::convertible_to< bool                                >;
      { g.adapt()                 } -> Std::convertible_to< bool                                >;
      { g.comm()                  } -> Std::convertible_to< typename G::CollectiveCommunication >;
      { g.loadBalance()           } -> Std::convertible_to< bool                                >;
      // requireConvertible<bool>(g.loadBalance(/*data*/ std::declval<DataHandle&>())) // FIXME use a default handler to instantiate this function
      g.globalRefine(refCount);
      g.postAdapt();
    };

#endif

    namespace Fallback {

      struct Grid
      {
        template<class G>
        auto require(G&& g) -> decltype(
          requireConvertible<int>(G::dimension),
          requireConvertible<int>(G::dimensionworld),
          requireConcept<GridView,typename G::LeafGridView>(),
          requireConcept<GridView,typename G::LevelGridView>(),
          requireConcept<Intersection,typename G::LeafIntersection>(),
          requireConcept<Intersection,typename G::LevelIntersection>(),
          requireConcept<IntersectionIterator,typename G::LeafIntersectionIterator>(),
          requireConcept<IntersectionIterator,typename G::LevelIntersectionIterator>(),
          requireConcept<IndexSet,typename G::LevelIndexSet>(),
          requireConcept<IndexSet,typename G::LeafIndexSet>(),
          requireConcept<IdSet,typename G::GlobalIdSet>(),
          requireConcept<IdSet,typename G::LocalIdSet>(),
          requireConcept<GridCodim<G::dimension>,G>(), // Start recursion on codims
          requireTrue<std::is_same<G,typename G::LeafGridView::Grid>::value>(),
          requireTrue<std::is_same<G,typename G::LevelGridView::Grid>::value>(),
          requireType<typename G::ctype>(),
          requireType<typename G::HierarchicIterator>(),
          requireConvertible< int                                 >( g.maxLevel()                                                                                       ),
          requireConvertible< int                                 >( g.size(/*level*/ int{}, /*codim*/ int{})                                                           ),
          requireConvertible< int                                 >( g.size(/*codim*/ int{})                                                                            ),
          requireConvertible< int                                 >( g.size(/*level*/ int{}, /*type*/ GeometryType{})                                                   ),
          requireConvertible< int                                 >( g.size(/*type*/ GeometryType{})                                                                    ),
          requireConvertible< std::size_t                         >( g.numBoundarySegments()                                                                            ),
          requireConvertible< typename G::LevelGridView           >( g.levelGridView(/*level*/ int{})                                                                   ),
          requireConvertible< typename G::LeafGridView            >( g.leafGridView()                                                                                   ),
          requireConvertible< const typename G::GlobalIdSet&      >( g.globalIdSet()                                                                                    ),
          requireConvertible< const typename G::LocalIdSet&       >( g.localIdSet()                                                                                     ),
          requireConvertible< const typename G::LevelIndexSet&    >( g.levelIndexSet(/*level*/ int{})                                                                   ),
          requireConvertible< const typename G::LeafIndexSet&     >( g.leafIndexSet()                                                                                   ),
          requireConvertible< bool                                >( g.mark(/*refCount*/ int{}, /*entity*/std::declval<const typename G::template Codim<0>::Entity&>()) ),
          requireConvertible< int                                 >( g.getMark(std::declval<const typename G::template Codim<0>::Entity&>())                            ),
          requireConvertible< bool                                >( g.preAdapt()                                                                                       ),
          requireConvertible< bool                                >( g.adapt()                                                                                          ),
          requireConvertible< typename G::CollectiveCommunication >( g.comm()                                                                                           ),
          requireConvertible< bool                                >( g.loadBalance()                                                                                    ),
          g.globalRefine(/*refCount*/ int{}),
          g.postAdapt()
          // requireConvertible<bool>(g.loadBalance(/*data*/ std::declval<DataHandle&>())) // FIXME use a default handler to instantiate this function
        );
      };
    } // nampespace Fallback
  } // nampespace Concept

  template <class G>
  constexpr void expectGrid()
  {
#if DUNE_HAVE_CXX_CONCEPTS
    static_assert(Concept::Grid<G>);
#else
    static_assert(models<Concept::Fallback::Grid, G>());
#endif
  }

} // end namespace Dune

#endif
