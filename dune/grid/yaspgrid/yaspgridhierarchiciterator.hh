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
  class YaspHierarchicIterator :
    public YaspEntityPointer<0,GridImp>
  {
    enum { dim=GridImp::dimension };

    typedef YaspEntity<0,GridImp::dimension,GridImp> YaspEntityImp;

  public:
    // types used from grids
    typedef typename GridImp::YGridLevelIterator YGLI;
    typedef typename GridImp::YGrid::Iterator I;
    typedef typename GridImp::template Codim<0>::Entity Entity;

    //! constructor
    YaspHierarchicIterator (const YGLI& g, const I& it, int maxlevel) :
      YaspEntityPointer<0,GridImp>(g,it)
    {
      // store reference to entity implementation for better readability
      YaspEntityImp& entity = entityImplementation();
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

    //! constructor
    YaspHierarchicIterator (const YaspHierarchicIterator& it) :
      YaspEntityPointer<0,GridImp>(it),
      _maxlevel(it._maxlevel), stack(it.stack)
    {}

    //! increment
    void increment ()
    {
      // sanity check: do nothing when stack is empty
      if (stack.empty()) return;

      // if maxlevel not reached then push sons
      if (entityImplementation()._g->level()<_maxlevel)
        push_sons();

      // in any case pop one element
      pop_tos();
    }

    void print (std::ostream& s) const
    {
      // store reference to entity implementation for better readability
      YaspEntityImp& entity = entityImplementation();
      s << "HIER: " << "level=" << entity._g.level()
        << " position=" << entity._it.coord()
        << " superindex=" << entity._it.superindex()
        << " maxlevel=" << entity._maxlevel
        << " stacksize=" << stack.size()
        << std::endl;
    }

  private:
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
      YaspEntityImp& entity = entityImplementation();

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

    YaspEntityImp& entityImplementation()
    {
      return GridImp::getRealImplementation(this->_entity);
    }

    // make TOS the current element
    void pop_tos ()
    {
      StackElem se = stack.top();
      stack.pop();
      entityImplementation()._g = se.g;
      entityImplementation()._it.reinit(entityImplementation()._g->overlap[0],se.coord);
    }
  };

} // namespace Dune

#endif  //  DUNE_GRID_YASPGRIDHIERARCHICITERATOR_HH
