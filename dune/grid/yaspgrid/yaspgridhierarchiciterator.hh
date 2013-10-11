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
  public:
    // types used from grids
    typedef typename GridImp::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,typename GridImp::ctype>::TransformingSubIterator TSI;
    typedef typename GridImp::template Codim<0>::Entity Entity;

    //! constructor
    YaspHierarchicIterator (const GridImp* yg, const YGLI& g, const TSI& it, int maxlevel) :
      YaspEntityPointer<0,GridImp>(yg,g,it)
    {
      // now iterator points to current cell
      StackElem se(this->_g);
      std::copy(this->_it.coord().begin(), this->_it.coord().end(), se.coord.begin());
      stack.push(se);

      // determine maximum level
      _maxlevel = std::min(maxlevel,this->_g->mg->maxLevel());

      // if maxlevel not reached then push yourself and sons
      if (this->_g->level()<_maxlevel)
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
      if (this->_g->level()<_maxlevel)
        push_sons();

      // in any case pop one element
      pop_tos();
    }

    void print (std::ostream& s) const
    {
      s << "HIER: " << "level=" << this->_g.level()
        << " position=" << this->_it.coord()
        << " superindex=" << this->_it.superindex()
        << " maxlevel=" << this->_maxlevel
        << " stacksize=" << stack.size()
        << std::endl;
    }

  private:
    int _maxlevel;       //!< maximum level of elements to be processed

    struct StackElem {
      YGLI g;         // grid level of the element
      array<int,dim> coord;   // and the coordinates
      StackElem(YGLI gg) : g(gg) {}
    };
    std::stack<StackElem> stack;    //!< stack holding elements to be processed

    // push sons of current element on the stack
    void push_sons ()
    {
      // yes, process all 1<<dim sons
      YGLI finer = this->_g;
      ++finer;
      StackElem se(finer);
      for (int i=0; i<(1<<dim); i++)
      {
        for (int k=0; k<dim; k++)
          if (i&(1<<k))
            se.coord[k] = this->_it.coord(k)*2+1;
          else
            se.coord[k] = this->_it.coord(k)*2;
        stack.push(se);
      }
    }

    // make TOS the current element
    void pop_tos ()
    {
      StackElem se = stack.top();
      stack.pop();
      this->_g = se.g;
      this->_it.reinit(this->_g->cell_overlap,se.coord);
    }
  };

} // namespace Dune

#endif  //  DUNE_GRID_YASPGRIDHIERARCHICITERATOR_HH
