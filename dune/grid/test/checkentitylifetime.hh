// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKENTITYLIFETIME_HH
#define DUNE_GRID_TEST_CHECKENTITYLIFETIME_HH

/** \file
    \brief Tests that make sure range-based entity iteration works correctly
           and that copied entities have the correct lifetime

 */

#include <cassert>
#include <limits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/rangegenerators.hh>

#if not defined(DUNE_ENTITY_LIFETIME_CHECK_ELEMENT_COUNT)
#define DUNE_ENTITY_LIFETIME_CHECK_ELEMENT_COUNT 32
#endif

template<typename GV, int codim>
bool checkEntityLifetimeForCodim(GV gv, std::size_t check_element_count, Dune::Codim<codim>, std::true_type)
{

  using namespace Dune;

  std::cout << "Lifetime / consistency check for entities, codim " << codim << std::endl;

  if (check_element_count > static_cast<std::size_t>(gv.size(codim)))
    {
      std::cout << "WARNING! Requested check of first " << check_element_count
                << " entities, but grid view only contains " << gv.size(codim)
                << " entities" << std::endl;
      check_element_count = gv.size(codim);
    }

  auto& index_set = gv.indexSet();
  auto& id_set = gv.grid().localIdSet();

  std::vector<typename GV::IndexSet::IndexType> indices;
  std::vector<typename GV::Grid::LocalIdSet::IdType> ids;
  std::vector<typename GV::template Codim<codim>::Entity> entity_list;
  std::vector<typename GV::template Codim<codim>::Entity::Geometry::GlobalCoordinate> coords;

  // store indices + ids + entities + coordinates
  {
    std::size_t i = 0;
    for (const auto& e : entities(gv,Dune::Codim<codim>()))
      {
        if (++i > check_element_count)
          break;
        indices.push_back(index_set.index(e));
        ids.push_back(id_set.id(e));
        entity_list.push_back(e);
        coords.push_back(e.geometry().corner(0));
      }
  }

  // check for consistency
  for (std::size_t i = 0; i < check_element_count; ++i)
    {
      if (index_set.index(entity_list[i]) != indices[i])
        DUNE_THROW(
          InvalidStateException,
          "ERROR! inconsistent index for entity " << i <<
          " (" << index_set.index(entity_list[i]) << " != " << indices[i] << ")");
      if (id_set.id(entity_list[i]) != ids[i])
        DUNE_THROW(
          InvalidStateException,
          "ERROR! inconsistent id for entity " << i <<
          " (" << id_set.id(entity_list[i]) << " != " << ids[i] << ")");
      if ((entity_list[i].geometry().corner(0) - coords[i]).two_norm() > std::numeric_limits<typename GV::ctype>::epsilon())
        DUNE_THROW(
          InvalidStateException,
          "ERROR! inconsistent corner(0) coordinate for entity " << i <<
          " (" << entity_list[i].geometry().corner(0) << " != " << coords[i] << ")");
    }

  return true;
}

template<typename GV, int codim>
bool checkEntityLifetimeForCodim(GV, const std::size_t,
                                 Dune::Codim<codim>, std::false_type)
{
  std::cout << "SKIPPING lifetime / consistency check for missing entities, codim " << codim << std::endl;
  return false;
}

namespace {

  // helper stuff for checking all codims

  template<std::size_t... i>
  struct index_pack
  {};

  //! TMP to build an index_pack containing the sequence 0,...,n-1.
  template<std::size_t n, std::size_t... i>
  struct index_pack_builder
    : public index_pack_builder<n-1,n-1,i...>
  {};

  // end of recursion
  template<std::size_t... i>
  struct index_pack_builder<0,i...>
  {
    typedef index_pack<i...> type;
  };

  template<typename GV>
  typename index_pack_builder<GV::dimension + 1>::type
  create_codims(GV)
  {
    return {};
  }

  template<typename... T>
  void invoke(T&&...)
  {}

  template<typename GV, std::size_t... codim>
  void do_check_entity_lifetime(GV gv, const std::size_t check_element_count, index_pack<codim...>)
  {
    invoke(
      checkEntityLifetimeForCodim(
        gv,
        check_element_count,
        Dune::Codim<codim>(),
        std::integral_constant<
        bool,
        Dune::Capabilities::hasEntity<typename GV::Grid,codim>::v && Dune::Capabilities::hasEntityIterator<typename GV::Grid,codim>::v
        >()
        )...
      );
  }

}


template<typename GV>
void checkEntityLifetime(GV gv, const std::size_t check_element_count = 32)
{
  do_check_entity_lifetime(gv,check_element_count,create_codims(gv));
}

#endif // #ifndef DUNE_GRID_TEST_CHECKENTITYLIFETIME_HH
