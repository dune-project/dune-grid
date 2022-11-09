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
 * @ingroup Grid
 * @{
 * @par Description
 *  This group gathers several concepts related to grids.
 * @}
 */

namespace Dune::Concept {
namespace Impl {

  template<class G, int codim, Dune::PartitionIteratorType partition>
  concept GridCodimPartition =
    EntityIterator<typename G::template Codim<codim>::template Partition<partition>::LevelIterator> &&
    EntityIterator<typename G::template Codim<codim>::template Partition<partition>::LeafIterator>;

  template<class G, int codim>
  concept GridCodimAllPartitions =
    GridCodimPartition<G,codim,Dune::PartitionIteratorType::Interior_Partition> &&
    GridCodimPartition<G,codim,Dune::PartitionIteratorType::InteriorBorder_Partition> &&
    GridCodimPartition<G,codim,Dune::PartitionIteratorType::Overlap_Partition> &&
    GridCodimPartition<G,codim,Dune::PartitionIteratorType::OverlapFront_Partition> &&
    GridCodimPartition<G,codim,Dune::PartitionIteratorType::All_Partition> &&
    GridCodimPartition<G,codim,Dune::PartitionIteratorType::Ghost_Partition>;

  template<class G, int codim>
  concept GridCodim =
    Geometry<typename G::template Codim<codim>::Geometry> &&
    Geometry<typename G::template Codim<codim>::LocalGeometry> &&
    Entity<typename G::template Codim<codim>::Entity> &&
    EntitySeed<typename G::template Codim<codim>::EntitySeed> &&
  requires(const G cg, const typename G::template Codim<codim>::EntitySeed& seed)
  {
    { cg.entity(seed) } -> std::convertible_to<typename G::template Codim<codim>::Entity>;

    requires (not Dune::Capabilities::canCommunicate<G,codim>::v) ||
      requires(G g, Archetypes::CommDataHandle<std::byte>& handle)
      {
        { g.loadBalance(handle) } -> std::convertible_to<bool>;
      };
  } && GridCodimAllPartitions<G,codim>;

  template<class G, int first, int last>
  concept GridAllCodims = requires(std::make_integer_sequence<int,last-first> codims)
  {
    []<int... c>(std::integer_sequence<int,c...>)
        requires (GridCodim<G,(first+c)> &&...) {} (codims);
  };

} // end namespace Impl


/**
 * \brief Requirements for implementations of the Dune::Grid interface.
 * \ingroup GridConcepts
 *
 * The `Grid` concept defines interface requirements of a parallel, in general
 * nonconforming, locally refined and hierarchical finite element mesh.
 * It consists of sub-concepts for `Dune::Concept::GridView`, `Dune::Concept::IndexSet`,
 * `Dune::Concept::IdSet`, and `Dune::Concept::Intersection`.
 *
 * See \ref Dune::Grid for an "abstract" interface definition of this concept.
 *
 * \par Models:
 * - `Dune::AlbertaGrid<dim,dow>`
 * - `Dune::GeometryGrid<G,F>` if `G` is a model of `Dune::Concept::Grid`.
 * - `Dune::IdentityGrid<G>` if `G` is a model of `Dune::Concept::Grid`.
 * - `Dune::OneDGrid`
 * - `Dune::UGGrid<dim>`
 * - `Dune::YaspGrid<dim, Coordinates>`
 *
 * \hideinitializer
 */
template<class G>
concept Grid =
  GridView<typename G::LeafGridView> &&
  GridView<typename G::LevelGridView> &&
  Intersection<typename G::LeafIntersection> &&
  Intersection<typename G::LevelIntersection> &&
  IntersectionIterator<typename G::LeafIntersectionIterator> &&
  IntersectionIterator<typename G::LevelIntersectionIterator> &&
  IndexSet<typename G::LevelIndexSet> &&
  IndexSet<typename G::LeafIndexSet> &&
  IdSet<typename G::GlobalIdSet> &&
  IdSet<typename G::LocalIdSet> &&
requires(const G cg, int level, int codim, Dune::GeometryType type)
{
  // static constants
  { G::dimension      } -> std::convertible_to<int>;
  { G::dimensionworld } -> std::convertible_to<int>;

  // type and concepts requirements
  requires std::same_as<G,typename G::LeafGridView::Grid>;
  requires std::same_as<G,typename G::LevelGridView::Grid>;
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
} &&
Impl::GridCodim<G,0> &&
Impl::GridAllCodims<G,1,G::dimension+1>;

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_GRID_HH
