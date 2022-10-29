// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_GRID_HH
#define DUNE_GRID_CONCEPTS_GRID_HH

#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

#include <dune/common/indices.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/entityiterator.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/concepts/gridview.hh>
#include <dune/grid/concepts/indexidset.hh>
#include <dune/grid/concepts/intersection.hh>
#include <dune/grid/concepts/intersectioniterator.hh>
#include <dune/grid/concepts/archetypes/datahandle.hh>


/*!@defgroup GridConcepts Grid Concepts
 * @{
 * @par Description
 *  This group gathers several concepts related to grids.
 * @}
 */

namespace Dune::Concept {
namespace Impl {

  template<class G, int codim,
    class Traits   = typename G::template Codim<codim>,
    class Interior = typename Traits::template Partition<Dune::PartitionIteratorType::Interior_Partition>,
    class IBorder  = typename Traits::template Partition<Dune::PartitionIteratorType::InteriorBorder_Partition>,
    class Overlap  = typename Traits::template Partition<Dune::PartitionIteratorType::Overlap_Partition>,
    class OFront   = typename Traits::template Partition<Dune::PartitionIteratorType::OverlapFront_Partition>,
    class All      = typename Traits::template Partition<Dune::PartitionIteratorType::All_Partition>,
    class Ghost    = typename Traits::template Partition<Dune::PartitionIteratorType::Ghost_Partition>>
  concept GridCodim = requires(const G cg, const typename Traits::EntitySeed& seed)
  {
    requires Geometry<typename Traits::Geometry>;
    requires Geometry<typename Traits::LocalGeometry>;
    requires EntityIterator<typename Traits::LevelIterator>;
    requires EntityIterator<typename Traits::LeafIterator>;
    requires EntityIterator<typename Interior::LevelIterator>;
    requires EntityIterator<typename Interior::LeafIterator>;
    requires EntityIterator<typename IBorder::LevelIterator>;
    requires EntityIterator<typename IBorder::LeafIterator>;
    requires EntityIterator<typename Overlap::LevelIterator>;
    requires EntityIterator<typename Overlap::LeafIterator>;
    requires EntityIterator<typename OFront::LevelIterator>;
    requires EntityIterator<typename OFront::LeafIterator>;
    requires EntityIterator<typename All::LevelIterator>;
    requires EntityIterator<typename All::LeafIterator>;
    requires EntityIterator<typename Ghost::LevelIterator>;
    requires EntityIterator<typename Ghost::LeafIterator>;

    requires Entity<typename Traits::Entity>;
    requires EntitySeed<typename Traits::EntitySeed>;
    { cg.entity(seed) } -> std::convertible_to<typename Traits::Entity>;

    requires (not Dune::Capabilities::canCommunicate<G,codim>::v) ||
      requires(G g, Archetypes::CommDataHandle<std::byte>& handle)
      {
        { g.loadBalance(handle) } -> std::convertible_to<bool>;
      };
  };


  template<class G, int dim>
  concept AllGridCodims = requires(std::make_integer_sequence<int,dim+1> codims)
  {
    []<int... cc>(std::integer_sequence<int,cc...>)
        requires (GridCodim<G,cc> &&...) {} (codims);
  };

} // end namespace Impl


/**
 * @brief Model of a grid
 * @ingroup GridConcepts
 * @details Dune::Grid is a template for this model
 */
template<class G>
concept Grid = requires(const G cg, int level, int codim, Dune::GeometryType type)
{
  // static constants
  { G::dimension      } -> std::convertible_to<int>;
  { G::dimensionworld } -> std::convertible_to<int>;

  // type and concepts requirements
  requires GridView<typename G::LeafGridView>;
  requires GridView<typename G::LevelGridView>;
  requires Intersection<typename G::LeafIntersection>;
  requires Intersection<typename G::LevelIntersection>;
  requires IntersectionIterator<typename G::LeafIntersectionIterator>;
  requires IntersectionIterator<typename G::LevelIntersectionIterator>;
  requires IndexSet<typename G::LevelIndexSet>;
  requires IndexSet<typename G::LeafIndexSet>;
  requires IdSet<typename G::GlobalIdSet>;
  requires IdSet<typename G::LocalIdSet>;
  requires std::same_as<G,typename G::LeafGridView::Grid>;
  requires std::same_as<G,typename G::LevelGridView::Grid>;
  requires Impl::AllGridCodims<G,G::dimension>;
  typename G::ctype;
  typename G::HierarchicIterator;

  // const methods
  { cg.maxLevel()            } -> std::convertible_to<int>;
  { cg.size(level, codim)    } -> std::convertible_to<int>;
  { cg.size(codim)           } -> std::convertible_to<int>;
  { cg.size(level, type)     } -> std::convertible_to<int>;
  { cg.size(type)            } -> std::convertible_to<int>;
  { cg.numBoundarySegments() } -> std::convertible_to<std::size_t>;
  { cg.levelGridView(level)  } -> std::convertible_to<typename G::LevelGridView>;
  { cg.leafGridView()        } -> std::convertible_to<typename G::LeafGridView>;
  { cg.globalIdSet()         } -> std::convertible_to<const typename G::GlobalIdSet&>;
  { cg.localIdSet()          } -> std::convertible_to<const typename G::LocalIdSet&>;
  { cg.levelIndexSet(level)  } -> std::convertible_to<const typename G::LevelIndexSet&>;
  { cg.leafIndexSet()        } -> std::convertible_to<const typename G::LeafIndexSet&>;
  { cg.comm()                } -> std::convertible_to<typename G::CollectiveCommunication>;

  // mutable methods
  requires requires(G g, int refCount, const typename G::template Codim<0>::Entity& entity)
  {
    { g.mark(refCount,entity)  } -> std::convertible_to<bool>;
    { g.getMark(entity)        } -> std::convertible_to<int>;
    { g.preAdapt()             } -> std::convertible_to<bool>;
    { g.adapt()                } -> std::convertible_to<bool>;
    { g.loadBalance()          } -> std::convertible_to<bool>;
    g.globalRefine(refCount);
    g.postAdapt();
  };
};

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_GRID_HH
