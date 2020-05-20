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


namespace Dune {
  namespace Concept
  {
    template<int codim>
    struct GridViewCodim : public Refines<GridViewCodim<codim-1>>
    {
      template<class GV>
      auto require(GV&& gv) -> decltype(
        requireConcept<Dune::Concept::EntityIterator, typename GV::template Codim<codim>::Iterator>(),
        requireConcept<Dune::Concept::Entity,typename GV::template Codim<codim>::Entity>(),
        requireConcept<Dune::Concept::Geometry,typename GV::template Codim<codim>::Geometry>(),
        requireConcept<Dune::Concept::Geometry,typename GV::template Codim<codim>::LocalGeometry>(),
        requireConcept<Dune::Concept::EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::Iterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::Iterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::Iterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::Iterator>(),
        requireConcept<Dune::Concept::EntityIterator, typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::Iterator>(),
        requireConvertible<typename GV::template Codim<codim>::Iterator>(gv.template begin<codim>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator>(gv.template begin<codim,Dune::PartitionIteratorType::Interior_Partition>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::Iterator>(gv.template begin<codim,Dune::PartitionIteratorType::InteriorBorder_Partition>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::Iterator>(gv.template begin<codim,Dune::PartitionIteratorType::Overlap_Partition>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::Iterator>(gv.template begin<codim,Dune::PartitionIteratorType::OverlapFront_Partition>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::Iterator>(gv.template begin<codim,Dune::PartitionIteratorType::All_Partition>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::Iterator>(gv.template begin<codim,Dune::PartitionIteratorType::Ghost_Partition>()),
        requireConvertible<typename GV::template Codim<codim>::Iterator>(gv.template end<codim>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Interior_Partition>::Iterator>(gv.template end<codim,Dune::PartitionIteratorType::Interior_Partition>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>::Iterator>(gv.template end<codim,Dune::PartitionIteratorType::InteriorBorder_Partition>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Overlap_Partition>::Iterator>(gv.template end<codim,Dune::PartitionIteratorType::Overlap_Partition>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>::Iterator>(gv.template end<codim,Dune::PartitionIteratorType::OverlapFront_Partition>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::All_Partition>::Iterator>(gv.template end<codim,Dune::PartitionIteratorType::All_Partition>()),
        requireConvertible<typename GV::template Codim<codim>::template Partition<Dune::PartitionIteratorType::Ghost_Partition>::Iterator>(gv.template end<codim,Dune::PartitionIteratorType::Ghost_Partition>())
      );
    };

    // stop recursion
    template<>
    struct GridViewCodim<-1> : public AnyType {};

    struct GridView
    {
      template<class GV>
      auto require(GV&& gv) -> decltype(
        requireType<typename GV::Traits>(),
        requireType<typename GV::Grid>(),
        requireType<typename GV::CollectiveCommunication>(),
        requireType<typename GV::ctype>(),
        requireConcept<Dune::Concept::IndexSet,typename GV::IndexSet>(),
        requireConcept<Dune::Concept::Intersection,typename GV::Intersection>(),
        requireConcept<Dune::Concept::IntersectionIterator,typename GV::IntersectionIterator>(),
        requireConcept<Dune::Concept::GridViewCodim<GV::dimension>,GV>(),
        requireConvertible<bool>(GV::conforming),
        requireConvertible<int>(GV::dimension),
        requireConvertible<int>(GV::dimensionworld),
        requireConvertible<const typename GV::Grid&>(gv.grid()),
        requireConvertible<const typename GV::IndexSet&>(gv.indexSet()),
        requireConvertible<int>(gv.size(/*codim*/ int{} )),
        requireConvertible<int>(gv.size(/*type*/ Dune::GeometryType{} )),
        requireConvertible<typename GV::IntersectionIterator>(gv.ibegin(/*entity*/ std::declval<const typename GV::template Codim<0>::Entity&>() )),
        requireConvertible<typename GV::IntersectionIterator>(gv.iend(/*entity*/ std::declval<const typename GV::template Codim<0>::Entity&>() )),
        requireConvertible<typename GV::CollectiveCommunication>(gv.comm()),
        requireConvertible<int>(gv.overlapSize(/*codim*/ int{} )),
        requireConvertible<int>(gv.ghostSize(/*codim*/ int{} )),
        // gv.communicate(std::declval<Handler&>(),std::declval<Dune::InterfaceType>(), std::declval<CommunicationDirection>()), // FIXME use a default handler to instantiate this function
        GV{gv},
        gv = gv
      );
    };

  }

  template <class GV>
  constexpr void expectGridView()
  {
    static_assert(models<Concept::GridView, GV>());
  }

}  // end namespace Dune

#endif
