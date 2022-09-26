// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDHIERARCHICITERATOR_HH
#define DUNE_GRID_YASPGRIDHIERARCHICITERATOR_HH

/** \file
 * The YaspHierarchicIterator class
 *
 * Enables iteration over son entities of codim 0
 */

namespace Dune {

  /** \brief YaspHierarchicIterator enables iteration over son entities of codim 0
   */
  template<class GridImp>
  class YaspHierarchicIterator
  {
    constexpr static int dim = GridImp::dimension;

    typedef YaspEntity<0,GridImp::dimension,GridImp> YaspEntityImp;

  public:
    // types used from grids
    typedef typename GridImp::YGridLevelIterator YGLI;
    typedef typename GridImp::YGrid::Iterator I;
    typedef typename GridImp::template Codim<0>::Entity Entity;

    //! default constructor creating empty iterator
    YaspHierarchicIterator () : _entity(), _maxlevel(-1), stack() {}

    //! constructor
    YaspHierarchicIterator (const YGLI& g, const I& it, int maxlevel) :
      _entity(YaspEntity<0, dim, GridImp>(g,it))
    {
      // store reference to entity implementation for better readability
      YaspEntityImp& entity = _entity.impl();
      // now iterator points to current cell
      StackElem se(entity._g);
      std::copy(entity._it.coord().begin(), entity._it.coord().end(), se.coord.begin());
      stack.push(se);

      // determine maximum level
      _maxlevel = std::min(maxlevel,entity._g->mg->maxLevel());

      // if maxlevel not reached then push yourself and sons
      if (entity._g->level()<_maxlevel)
      {
        push_sons();
      }

      // and make iterator point to first son if stack is not empty
      if (!stack.empty())
        pop_tos();
    }

    //! increment
    void increment ()
    {
      // sanity check: do nothing when stack is empty
      if (stack.empty()) return;

      // if maxlevel not reached then push sons
      if (_entity.impl()._g->level()<_maxlevel)
        push_sons();

      // in any case pop one element
      pop_tos();
    }

    //! equality
    bool equals (const YaspHierarchicIterator& rhs) const
    {
      return (_entity == rhs._entity);
    }

    //! dereferencing
    const Entity& dereference() const
    {
      return _entity;
    }

    void print (std::ostream& s) const
    {
      // store reference to entity implementation for better readability
      YaspEntityImp& entity = _entity.impl();
      s << "HIER: " << "level=" << entity._g.level()
        << " position=" << entity._it.coord()
        << " superindex=" << entity._it.superindex()
        << " maxlevel=" << entity._maxlevel
        << " stacksize=" << stack.size()
        << std::endl;
    }

  private:
    Entity _entity; //!< entity

    int _maxlevel;       //!< maximum level of elements to be processed

    struct StackElem {
      YGLI g;         // grid level of the element
      std::array<int,dim> coord;   // and the coordinates
      StackElem(YGLI gg) : g(gg) {}
    };
    std::stack<StackElem> stack;    //!< stack holding elements to be processed

    // push sons of current element on the stack
    void push_sons ()
    {
      // store reference to entity implementation for better readability
      YaspEntityImp& entity = _entity.impl();

      // yes, process all 1<<dim sons
      YGLI finer = entity._g;
      ++finer;
      StackElem se(finer);
      for (int i=0; i<(1<<dim); i++)
      {
        for (int k=0; k<dim; k++)
          if (i&(1<<k))
            se.coord[k] = entity._it.coord(k)*2+1;
          else
            se.coord[k] = entity._it.coord(k)*2;
        // not all entities have 2^d subentities due to refineOptions with keep_ovlp==false
        bool exists = true;
        for (int k=0; k<dim; k++)
          if ((se.coord[k] < finer->overlap[0].dataBegin()->origin(k)) || (se.coord[k] >= finer->overlap[0].dataBegin()->origin(k)+finer->overlap[0].dataBegin()->size(k)))
            exists = false;
        if (exists)
          stack.push(se);
      }
    }

    // make TOS the current element
    void pop_tos ()
    {
      StackElem se = stack.top();
      stack.pop();
      YaspEntityImp& entity = _entity.impl();
      entity._g = se.g;
      entity._it.reinit(entity._g->overlap[0],se.coord);
    }
  };

} // namespace Dune

#endif  //  DUNE_GRID_YASPGRIDHIERARCHICITERATOR_HH
