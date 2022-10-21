// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_GRIDVIEW_HH
#define DUNE_GRID_CONCEPTS_GRIDVIEW_HH

#include <dune/grid/concepts/entityset.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/common/indices.hh>

namespace Dune::Concept {

namespace Impl {

  template<class ES, int codim,  Dune::PartitionIteratorType partition>
  concept EntityPartitionSpan = requires(ES es) {
    requires EntityIterator<typename ES::template Codim<codim>::template Partition<partition>::Iterator>;
    { es.template begin<codim,partition>()  } -> std::convertible_to< typename ES::template Codim<codim>::template Partition<partition>::Iterator>;
    { es.template end<codim,partition>()    } -> std::convertible_to< typename ES::template Codim<codim>::template Partition<partition>::Iterator>;
  };


  template<class ES, int codim>
  concept EntitySetAllPartitions = requires(ES es) {
    requires EntitySet<ES,codim>;
    requires EntityPartitionSpan<ES,codim,Dune::PartitionIteratorType::InteriorBorder_Partition>;
    requires EntityPartitionSpan<ES,codim,Dune::PartitionIteratorType::Overlap_Partition>;
    requires EntityPartitionSpan<ES,codim,Dune::PartitionIteratorType::OverlapFront_Partition>;
    requires EntityPartitionSpan<ES,codim,Dune::PartitionIteratorType::All_Partition>;
    requires EntityPartitionSpan<ES,codim,Dune::PartitionIteratorType::Ghost_Partition>;
  };

  template<class ES, class Grid, int codim>
    requires Dune::Capabilities::hasEntityIterator<Grid,codim>::v
  void requireEntitySetAllPartitions()
    requires EntitySetAllPartitions<ES,codim> {}

  template<class ES, class Grid, int codim>
    requires (not Dune::Capabilities::hasEntityIterator<Grid,codim>::v)
  void requireEntitySetAllPartitions() {}

}

/**
 * @brief Model of a grid view
 * @ingroup GridConcepts
 * @details Dune::GridView is a template for this model
 */
template<class GV>
concept GridView = requires(std::make_integer_sequence<int,GV::dimension+1> dims)
{
  requires EntitySet<GV,0>;

  []<int... d>(std::integer_sequence<int,d...>) requires
    requires { (Impl::requireEntitySetAllPartitions<GV,typename GV::Grid,(GV::dimension-d)>(),...); }{} (dims);
};

}  // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_GRIDVIEW_HH
