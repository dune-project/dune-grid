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


namespace Dune {
  namespace Concept
  {
    template<int codim>
    struct GridCodim : public Refines<GridCodim<codim-1>>
    {
      template<class G>
      auto require(G&& g) -> decltype(
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::Iterator>(),
        requireConcept<Dune::Concept::Entity,typename G::template Codim<codim>::Entity>(),
        requireConcept<Dune::Concept::EntitySeed,typename G::template Codim<codim>::EntitySeed>(),
        requireConcept<Dune::Concept::Geometry,typename G::template Codim<codim>::Geometry>(),
        requireConcept<Dune::Concept::Geometry,typename G::template Codim<codim>::LocalGeometry>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::LevelIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::LeafIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::LevelIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::LeafIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::LevelIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::LeafIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::LevelIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::LeafIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::LevelIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::LeafIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::LevelIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::LeafIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::LevelIterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename G::template Codim<codim>::LeafIterator>(),
        requireConvertible<typename G::template Codim<codim>::Iterator>(g.template begin<codim>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator>(g.template begin<codim,Dune::PartitionIteratorType::Interior_Partition>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::Iterator>(g.template begin<codim,Dune::PartitionIteratorType::InteriorBorder_Partition>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::Iterator>(g.template begin<codim,Dune::PartitionIteratorType::Overlap_Partition>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::Iterator>(g.template begin<codim,Dune::PartitionIteratorType::OverlapFront_Partition>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::Iterator>(g.template begin<codim,Dune::PartitionIteratorType::All_Partition>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::Iterator>(g.template begin<codim,Dune::PartitionIteratorType::Ghost_Partition>()),
        requireConvertible<typename G::template Codim<codim>::Iterator>(g.template end<codim>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator>(g.template end<codim,Dune::PartitionIteratorType::Interior_Partition>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::Iterator>(g.template end<codim,Dune::PartitionIteratorType::InteriorBorder_Partition>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::Iterator>(g.template end<codim,Dune::PartitionIteratorType::Overlap_Partition>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::Iterator>(g.template end<codim,Dune::PartitionIteratorType::OverlapFront_Partition>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::Iterator>(g.template end<codim,Dune::PartitionIteratorType::All_Partition>()),
        requireConvertible<typename G::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::Iterator>(g.template end<codim,Dune::PartitionIteratorType::Ghost_Partition>()),
        requireConvertible<typename G::template Codim<codim>::Entity>(g.entity(/*seed*/ std::declval<const typename G::template Codim<codim>::EntitySeed&>() ))
      );
    };

    // stop recursion
    template<>
    struct GridCodim<-1> : public AnyType {};

    struct Grid
    {
      template<class G>
      auto require(G&& g) -> decltype(
        requireConvertible<int>(G::dimension),
        requireConvertible<int>(G::dimensionworld),
        requireConcept<Dune::Concept::GridView,typename G::LeafGridView>(),
        requireConcept<Dune::Concept::GridView,typename G::LevelGridView>(),
        requireTrue<std::is_same<G,typename G::LeafGridView::Grid>::value>(),
        requireTrue<std::is_same<G,typename G::LevelGridView::Grid>::value>(),
        requireConcept<Dune::Concept::Intersection,typename G::LeafIntersection>(),
        requireConcept<Dune::Concept::Intersection,typename G::LevelIntersection>(),
        requireConcept<Dune::Concept::IntersectionIterator,typename G::LeafIntersectionIterator>(),
        requireConcept<Dune::Concept::IntersectionIterator,typename G::LevelIntersectionIterator>(),
        requireType<typename G::HierarchicIterator>(),
        requireConcept<Dune::Concept::IndexSet,typename G::LevelIndexSet>(),
        requireConcept<Dune::Concept::IndexSet,typename G::LeafIndexSet>(),
        requireConcept<Dune::Concept::IdSet,typename G::GlobalIdSet>(),
        requireConcept<Dune::Concept::IdSet,typename G::LocalIdSet>(),
        requireType<typename G::CollectiveCommunication>(),
        requireType<typename G::ctype>(),
        requireConvertible<int>(g.maxLevel()),
        requireConvertible<int>(g.size(/*level*/ int{}, /*codim*/ int{})),
        requireConvertible<int>(g.size(/*codim*/ int{})),
        requireConvertible<int>(g.size(/*level*/ int{}, /*type*/ GeometryType{})),
        requireConvertible<int>(g.size(/*type*/ GeometryType{})),
        requireConvertible<std::size_t>(g.numBoundarySegments()),
        requireConvertible<typename G::LevelGridView>(g.levelGridView(/*level*/ int{})),
        requireConvertible<typename G::LeafGridView>(g.leafGridView()),
        requireConvertible<const typename G::GlobalIdSet&>(g.globalIdSet()),
        requireConvertible<const typename G::LocalIdSet&>(g.localIdSet()),
        requireConvertible<const typename G::LevelIndexSet&>(g.levelIndexSet(/*level*/ int{})),
        requireConvertible<const typename G::LeafIndexSet&>(g.leafIndexSet()),
        g.globalRefine(/*refCount*/ int{}),
        requireConvertible<bool>(g.mark(/*refCount*/ int{}, /*entity*/std::declval<const typename G::template Codim<0>::Entity&>())),
        requireConvertible<int>(g.getMark(std::declval<const typename G::template Codim<0>::Entity&>())),
        requireConvertible<bool>(g.preAdapt()),
        requireConvertible<bool>(g.adapt()),
        g.postAdapt(),
        requireConvertible<typename G::CollectiveCommunication>(g.comm()),
        requireConvertible<bool>(g.loadBalance())
        // requireConvertible<bool>(g.loadBalance(/*data*/ std::declval<DataHandle&>())) // FIXME use a default handler to instantiate this function
      );
    };

  }

  template <class GV>
  constexpr void expectGrid()
  {
    static_assert(models<Concept::Grid, GV>());
  }

}  // end namespace Dune

#endif
