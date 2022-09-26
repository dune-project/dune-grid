// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKINTERSECTIONLIFETIME_HH
#define DUNE_GRID_TEST_CHECKINTERSECTIONLIFETIME_HH

/** \file
    \brief Tests that make sure range-based intersection iteration works correctly
           and that copied intersections have the correct lifetime

 */

#include <cassert>
#include <iostream>
#include <vector>

#include <dune/common/exceptions.hh>

template<typename GV>
void checkIntersectionLifetime(GV gv, std::size_t check_element_count = 32)
{

  using namespace Dune;

  std::cout << "Intersection Lifetime / consistency check" << std::endl;

  if (check_element_count > static_cast<std::size_t>(gv.size(0)))
    {
      std::cout << "WARNING! Requested check of intersections for first " << check_element_count
                << " elements, but grid view only contains " << gv.size(0)
                << " elements" << std::endl;
      check_element_count = gv.size(0);
    }

  auto& index_set = gv.indexSet();

  std::vector<std::vector<int> > indices;
  std::vector<typename GV::template Codim<0>::Entity> entity_list;
  std::vector<std::vector<typename GV::Intersection> > intersection_list;
  std::vector<std::vector<typename GV::Intersection::Geometry::GlobalCoordinate> > coords;

  // store indices + entities + intersections + coordinates
  {
    std::size_t i = 0;
    for (const auto& e : elements(gv))
      {
        if (++i > check_element_count)
          break;
        entity_list.push_back(e);
        indices.push_back({});
        intersection_list.push_back({});
        coords.push_back({});
        for (const auto& is : intersections(gv,e))
          {
            indices.back().push_back(is.indexInInside());
            intersection_list.back().push_back(is);
            coords.back().push_back(is.geometry().corner(0));
          }
      }
  }

  // check for consistency
  {
    std::size_t i = 0;
    for (const auto& e : elements(gv))
      {
        if (i >= check_element_count)
          break;
        if (e != entity_list[i])
          DUNE_THROW(
            InvalidStateException,
            "ERROR! Got different entity on second iteration for entity " << i <<
            " (" << index_set.index(entity_list[i]) << " != " << index_set.index(e) << ")");
        std::size_t j = 0;
        for (const auto& is : intersections(gv,e))
          {
            if (indices[i][j] != is.indexInInside())
              DUNE_THROW(
                InvalidStateException,
                "ERROR! inconsistent indexInInside() for intersection " << j << " of entity " << i <<
                " (" << indices[i][j] << " != " << is.indexInInside() << ")");
            if (intersection_list[i][j] != is)
              DUNE_THROW(
                InvalidStateException,
                "ERROR! inconsistent intersection equality for intersection " << j << " of entity " << i);
            if (entity_list[i] != is.inside())
              DUNE_THROW(
                InvalidStateException,
                "ERROR! inconsistent inside() for intersection " << j << " of entity " << i);
            if (coords[i][j] != is.geometry().corner(0))
              DUNE_THROW(
                InvalidStateException,
                "ERROR! inconsistent corner(0) coordinate for intersection " << j << " of entity " << i <<
                " (" << coords[i][j] << " != " << is.geometry().corner(0) << ")");
            ++j;
          }
        ++i;
      }
  }

}

#endif // #ifndef DUNE_GRID_TEST_CHECKENTITYLIFETIME_HH
