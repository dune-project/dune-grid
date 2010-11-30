// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#ifndef DUNE_GRIDINFO_HH
#define DUNE_GRIDINFO_HH

#include <iostream>
#include <typeinfo>
#include <dune/common/exceptions.hh>
#include <dune/common/classname.hh>
#include "grid.hh"
#include "genericreferenceelements.hh"

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
              << " geomTypes=(";
    bool first=true;
    for (int c=0; c<=dim; c++)
    {
      for (std::size_t i=0; i<grid.leafIndexSet().geomTypes(c).size(); i++)
      {
        if (!first) std::cout << ",";
        std::cout << grid.leafIndexSet().geomTypes(c)[i]
                  << "[" << c << "]"
                  << "=" << grid.leafIndexSet().size(grid.leafIndexSet().geomTypes(c)[i]);
        first=false;
      }
    }
    std::cout << ")" << std::endl;


    return;
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

    // the grid has an iterator providing the access to
    // all elements (better codim 0 entities) on a grid level
    // Note the use of the typename and template keywords.
    typedef typename G::Traits::template Codim<0>::LevelIterator LevelIterator;

    // print info about this level
    std::cout << prefix << "level=" << level
              << " dim=" << dim
              << " geomTypes=(";
    bool first=true;
    for (unsigned i=0; i<grid.levelIndexSet(level).geomTypes(0).size(); i++)
    {
      if (!first) std::cout << ",";
      std::cout << grid.levelIndexSet(level).geomTypes(0)[i]
                << "=" << grid.levelIndexSet(level).size(grid.levelIndexSet(level).geomTypes(0)[i]);
      first=false;
    }
    std::cout << ")" << std::endl;

    // print info about each element on given level
    LevelIterator eendit = grid.template lend<0>(level);
    for (LevelIterator it = grid.template lbegin<0>(level); it!=eendit; ++it)
    {
      std::cout << prefix << "level=" << it->level()
                << " " << it->type() << "[" << dim << "]"
                << " index=" << grid.levelIndexSet(level).index(*it)
                << " gid=" << grid.globalIdSet().template id<0>(*it)
                << " leaf=" << it->isLeaf()
                << " partition=" << PartitionName(it->partitionType())
                << " center=("
                << it->geometry().global(Dune::GenericReferenceElements<ct,dim>::general(it->type()).position(0,0))
                << ")"
                << " first=(" << it->geometry().corner(0) << ")"
                << std::endl;

      std::cout << prefix << "codim " << dim << " subindex";
      for (int i=0; i<it->template count<dim>(); i++)
      {
        std::cout << " " << i << ":" << grid.levelIndexSet(level).subIndex(*it,i,dim);
      }
      std::cout << std::endl;

      std::cout << prefix << "codim " << dim-1 << " subindex";
      for (int i=0; i<it->template count<dim-1>(); i++)
      {
        std::cout << " " << i << ":" << grid.levelIndexSet(level).subIndex(*it,i,dim-1);
      }
      std::cout << std::endl;

    }

    return;
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

    // the grid has an iterator providing the access to
    // all elements (better codim 0 entities) on a grid level
    // Note the use of the typename and template keywords.
    typedef typename G::Traits::template Codim<0>::LeafIterator LeafIterator;
    typedef typename G::Traits::template Codim<dim>::LeafIterator VLeafIterator;

    // print info about the leaf grid
    std::cout << prefix << "leaf"
              << " dim=" << dim
              << " geomTypes=(";
    bool first=true;
    for (int c=0; c<=dim; c++)
    {
      for (unsigned i=0; i<grid.leafIndexSet().geomTypes(c).size(); i++)
      {
        if (!first) std::cout << ",";
        std::cout << grid.leafIndexSet().geomTypes(c)[i]
                  << "[" << c << "]"
                  << "=" << grid.leafIndexSet().size(grid.leafIndexSet().geomTypes(c)[i]);
        first=false;
      }
    }
    std::cout << ")" << std::endl;

    // print info about nodes in leaf grid
    VLeafIterator veendit = grid.template leafend<dim>();
    for (VLeafIterator it = grid.template leafbegin<dim>(); it!=veendit; ++it)
    {
      std::cout << prefix << "level=" << it->level()
                << " " << it->type() << "[" << dim << "]"
                << " index=" << grid.leafIndexSet().index(*it)
                << " gid=" << grid.globalIdSet().template id<dim>(*it)
                << " partition=" << PartitionName(it->partitionType())
                << " pos=(" << it->geometry().corner(0) << ")"
                << std::endl;
    }

    // print info about each element in leaf grid
    LeafIterator eendit = grid.template leafend<0>();
    for (LeafIterator it = grid.template leafbegin<0>(); it!=eendit; ++it)
    {
      std::cout << prefix << "level=" << it->level()
                << " " << it->type() << "[" << dim << "]"
                << " index=" << grid.leafIndexSet().index(*it)
                << " gid=" << grid.globalIdSet().template id<0>(*it)
                << " leaf=" << it->isLeaf()
                << " partition=" << PartitionName(it->partitionType())
                << " center=("
                << it->geometry().global(Dune::GenericReferenceElements<ct,dim>::general(it->type()).position(0,0))
                << ")"
                << " first=(" << it->geometry().corner(0) << ")"
                << std::endl;

      std::cout << prefix << "codim " << dim << " subindex";
      for (int i=0; i<it->template count<dim>(); i++)
      {
        std::cout << " " << i << ":" << grid.leafIndexSet().subIndex(*it,i,dim);
      }
      std::cout << std::endl;

      std::cout << prefix << "codim " << dim-1 << " subindex";
      for (int i=0; i<it->template count<dim-1>(); i++)
      {
        std::cout << " " << i << ":" << grid.leafIndexSet().subIndex(*it,i,dim-1);
      }
      std::cout << std::endl;

    }

    return;
  }


  /** @} */

}
#endif
