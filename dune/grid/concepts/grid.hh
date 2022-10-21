// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_GRID_HH
#define DUNE_GRID_CONCEPTS_GRID_HH

#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/entityiterator.hh>
#include <dune/grid/concepts/intersection.hh>
#include <dune/grid/concepts/intersectioniterator.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/concepts/indexset.hh>
#include <dune/grid/concepts/gridview.hh>

#include <dune/grid/concepts/archetypes/datahandle.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/common/indices.hh>

#include <type_traits>
#include <cuchar>

/*!@defgroup GridConcepts Grid Concepts
 * @{
 * @par Description
 *  This group gathers several concepts related to grids.
 * @}
 */

namespace Dune::Concept {
  namespace Impl {

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
    { g.entity(seed) } -> std::convertible_to<typename G::template Codim<codim>::Entity                                                                                >;
  };


  template<class G, std::size_t dim>
  concept AllGridCodims = requires(std::make_index_sequence<dim+1> codims)
  {
    []<std::size_t... cc>(std::index_sequence<cc...>)
        requires (GridCodim<G,cc> &&...) {} (codims);
  };
}

/**
 * @brief Model of a grid
 * @ingroup GridConcepts
 * @details Dune::Grid is a template for this model
 */
template<class G>
concept Grid = requires(G g, int level, int codim, int refCount,
                        Dune::GeometryType type,
                        const typename G::template Codim<0>::Entity& entity)
{
  { G::dimension      } -> std::convertible_to<int>;
  { G::dimensionworld } -> std::convertible_to<int>;
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
  requires Impl::AllGridCodims<G,G::dimension>;
  typename G::ctype;
  typename G::HierarchicIterator;
  { g.maxLevel()              } -> std::convertible_to< int                                 >;
  { g.size(level, codim)      } -> std::convertible_to< int                                 >;
  { g.size(codim)             } -> std::convertible_to< int                                 >;
  { g.size(level, type)       } -> std::convertible_to< int                                 >;
  { g.size(type)              } -> std::convertible_to< int                                 >;
  { g.numBoundarySegments()   } -> std::convertible_to< std::size_t                         >;
  { g.levelGridView(level)    } -> std::convertible_to< typename G::LevelGridView           >;
  { g.leafGridView()          } -> std::convertible_to< typename G::LeafGridView            >;
  { g.globalIdSet()           } -> std::convertible_to< const typename G::GlobalIdSet&      >;
  { g.localIdSet()            } -> std::convertible_to< const typename G::LocalIdSet&       >;
  { g.levelIndexSet(level)    } -> std::convertible_to< const typename G::LevelIndexSet&    >;
  { g.leafIndexSet()          } -> std::convertible_to< const typename G::LeafIndexSet&     >;
  { g.mark(refCount,entity)   } -> std::convertible_to< bool                                >;
  { g.getMark(entity)         } -> std::convertible_to< int                                 >;
  { g.preAdapt()              } -> std::convertible_to< bool                                >;
  { g.adapt()                 } -> std::convertible_to< bool                                >;
  { g.comm()                  } -> std::convertible_to< typename G::CollectiveCommunication >;
  { g.loadBalance()           } -> std::convertible_to< bool                                >;
  requires requires(Archetypes::CommDataHandle<typename G::ctype>& handle,
                    InterfaceType iface, CommunicationDirection dir)
  {
    { g.loadBalance(handle) } -> std::convertible_to< bool >;
  };
  g.globalRefine(refCount);
  g.postAdapt();
};

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_GRID_HH
