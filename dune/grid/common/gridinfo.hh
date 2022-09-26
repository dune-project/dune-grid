// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_COMMON_GRIDINFO_HH
#define DUNE_GRID_COMMON_GRIDINFO_HH

#include <iostream>
#include <typeinfo>
#include <dune/common/exceptions.hh>
#include <dune/common/classname.hh>
#include <dune/geometry/referenceelements.hh>
#include "grid.hh"

/** @file
 * @author Peter Bastian
 * @brief Some functions to list information about a grid
 *
 */

namespace Dune
{
  /**
   * @addtogroup Grid
   *
   * @{
   */

  /** @brief A function to print some information about the grid as a whole.
   */
  template<class G>
  void gridinfo (const G& grid, std::string prefix="")
  {
    // first we extract the dimensions of the grid
    const int dim = G::dimension;
    const int dimworld = G::dimensionworld;

    // grid type and dimension
    std::cout << prefix << "=> " << className(grid)
              << " (dim=" << dim
              << ", dimworld=" << dimworld
              << ")" << std::endl;

    // level information
    for (int level=0; level<=grid.maxLevel(); level++)
    {
      std::cout << prefix << "level " << level;
      for (int cd=0; cd<=dim; cd++)
      {
        std::cout << " codim[" << cd << "]=" << grid.size(level,cd);
      }
      std::cout << std::endl;
    }

    // leaf information
    std::cout << prefix << "leaf   ";
    for (int cd=0; cd<=dim; cd++)
    {
      std::cout << " codim[" << cd << "]=" << grid.size(cd);
    }
    std::cout << std::endl;

    std::cout << prefix << "leaf"
              << " dim=" << dim
              << " types=(";
    bool first=true;
    for (int c=0; c<=dim; c++)
    {
      for (std::size_t i=0; i<grid.leafIndexSet().types(c).size(); i++)
      {
        if (!first) std::cout << ",";
        std::cout << grid.leafIndexSet().types(c)[i]
                  << "[" << c << "]"
                  << "=" << grid.leafIndexSet().size(grid.leafIndexSet().types(c)[i]);
        first=false;
      }
    }
    std::cout << ")" << std::endl;
  }


  /** @brief A function to print info about a grid level and its entities
   */
  template<class G>
  void gridlevellist (const G& grid, int level, std::string prefix)
  {
    // first we extract the dimensions of the grid
    const int dim = G::dimension;

    // type used for coordinates in the grid
    typedef typename G::ctype ct;

    // print info about this level
    std::cout << prefix << "level=" << level
              << " dim=" << dim
              << " types=(";
    bool first=true;
    for (unsigned i=0; i<grid.levelIndexSet(level).types(0).size(); i++)
    {
      if (!first) std::cout << ",";
      std::cout << grid.levelIndexSet(level).types(0)[i]
                << "=" << grid.levelIndexSet(level).size(grid.levelIndexSet(level).types(0)[i]);
      first=false;
    }
    std::cout << ")" << std::endl;

    // print info about each element on given level
    for (const auto& element : elements(levelGridView(grid, level)))
    {
      const auto& geometry = element.geometry();
      std::cout << prefix << "level=" << element.level()
                << " " << element.type() << "[" << dim << "]"
                << " index=" << grid.levelIndexSet(level).index(element)
                << " gid=" << grid.globalIdSet().template id<0>(element)
                << " leaf=" << element.isLeaf()
                << " partition=" << PartitionName(element.partitionType())
                << " center=("
                << geometry.global(Dune::ReferenceElements<ct,dim>::general(element.type()).position(0,0))
                << ")"
                << " first=(" << geometry.corner(0) << ")"
                << std::endl;

      std::cout << prefix << "codim " << dim << " subindex";
      for (unsigned int i=0; i < element.subEntities(dim); i++)
      {
        std::cout << " " << i << ":" << grid.levelIndexSet(level).subIndex(element,i,dim);
      }
      std::cout << std::endl;

      std::cout << prefix << "codim " << dim-1 << " subindex";
      for (unsigned int i=0; i < element.subEntities(dim-1); i++)
      {
        std::cout << " " << i << ":" << grid.levelIndexSet(level).subIndex(element,i,dim-1);
      }
      std::cout << std::endl;

    }
  }


  /** @brief A function to print info about a leaf grid and its entities
   */
  template<class G>
  void gridleaflist (const G& grid, std::string prefix)
  {
    // first we extract the dimensions of the grid
    const int dim = G::dimension;

    // type used for coordinates in the grid
    typedef typename G::ctype ct;

    // print info about the leaf grid
    std::cout << prefix << "leaf"
              << " dim=" << dim
              << " types=(";
    bool first=true;
    for (int c=0; c<=dim; c++)
    {
      for (unsigned i=0; i<grid.leafIndexSet().types(c).size(); i++)
      {
        if (!first) std::cout << ",";
        std::cout << grid.leafIndexSet().types(c)[i]
                  << "[" << c << "]"
                  << "=" << grid.leafIndexSet().size(grid.leafIndexSet().types(c)[i]);
        first=false;
      }
    }
    std::cout << ")" << std::endl;

    // print info about nodes in leaf grid
    for (const auto& vertex : vertices(leafGridView(grid)))
    {
      std::cout << prefix << "level=" << vertex.level()
                << " " << vertex.type() << "[" << dim << "]"
                << " index=" << grid.leafIndexSet().index(vertex)
                << " gid=" << grid.globalIdSet().template id<dim>(vertex)
                << " partition=" << PartitionName(vertex.partitionType())
                << " pos=(" << vertex.geometry().corner(0) << ")"
                << std::endl;
    }

    // print info about each element in leaf grid
    for (const auto& element : elements(leafGridView(grid)))
    {
      const auto& geometry = element.geometry();
      std::cout << prefix << "level=" << element.level()
                << " " << element.type() << "[" << dim << "]"
                << " index=" << grid.leafIndexSet().index(element)
                << " gid=" << grid.globalIdSet().template id<0>(element)
                << " leaf=" << element.isLeaf()
                << " partition=" << PartitionName(element.partitionType())
                << " center=("
                << geometry.global(Dune::ReferenceElements<ct,dim>::general(element.type()).position(0,0))
                << ")"
                << " first=(" << geometry.corner(0) << ")"
                << std::endl;

      std::cout << prefix << "codim " << dim << " subindex";
      for (unsigned int i=0; i < element.subEntities(dim); i++)
      {
        std::cout << " " << i << ":" << grid.leafIndexSet().subIndex(element,i,dim);
      }
      std::cout << std::endl;

      std::cout << prefix << "codim " << dim-1 << " subindex";
      for (unsigned int i=0; i < element.subEntities(dim-1); i++)
      {
        std::cout << " " << i << ":" << grid.leafIndexSet().subIndex(element,i,dim-1);
      }
      std::cout << std::endl;

    }
  }


  /** @} */

}
#endif
