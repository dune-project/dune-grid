// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_INTERSECTIONITERATOR_HH
#define DUNE_GRID_INTERSECTIONITERATOR_HH

#include <dune/common/iteratorfacades.hh>

#include <dune/grid/common/intersection.hh>

namespace Dune
{

  /** \brief Mesh entities of codimension 0 ("elements") allow to visit all
     intersections with "neighboring" elements and with the domain
     boundary.

     Template parameters are:

     - <tt>GridImp</tt> Type that is a model of Dune::Grid
     - <tt>IntersectionIteratorImp</tt> Class template that is a model of
     Dune::IntersectionIterator

     @warning the IntersectionIterator used to be both, Intersection and IntersectionIterator,
     at the same time. The two concepts are now properly separated. The IntersectionIterator
     still offers the old methods, but these are forwarded to the Intersection. All these methods
     are now marked deprecated.

     \deprecated All Intersection methods on the IntersectionIterator are deprecated,
       dereference IntersectionIterator to get the Intersection and call methods there.

     @warning The number of neigbors may be different from the number of
     faces/edges of an element!

     <h2>Overview</h2>

     Intersections are codimension 1 objects. These
     intersections are accessed via an IntersectionIterator. This allows
     the implementation of non-matching grids, as one face can now
     consist of several intersections.
     In a conforming mesh such an intersection corresponds to an entity of
     codimension 1 but in the general non-conforming case there will be no entity
     in the mesh that directly corresponds to the intersection. Thus, the
     IntersectionIterator describes these intersections implicitly.

     <H2>Engine Concept</H2>

     The IntersectionIterator class template wraps an object of type IntersectionIteratorImp
     and forwards all member
     function calls to corresponding members of this class. In that sense IntersectionIterator
     defines the interface and IntersectionIteratorImp supplies the implementation.

     <h2>Intersections, leaf grid and level grid</h2>

     On an entity \b e of codimension zero there are two ways to create
     IntersectionIterators by either using ilevelbegin() / ilevelend() or
     ileafbegin()/ileafend(). In the first case intersections with
     neighboring entities having the same level as \b e are traversed; in
     the second case  ileafbegin()==ileafend() if \b e is not a leaf otherwise
     all intersections with neighboring leaf entities are traversed.

     Consider a situation where two elements \b a and \b b have a common intersection.
     %Element \b b has been refined into an element \b c and \b d, while \b a has not
     been refined.
     In one space dimension this situation is depicted in the figure below.

     \image html  islocalref.png "IntersectionIterator in a locally refined mesh."
     \image latex islocalref.eps "IntersectionIterator in a locally refined mesh." width=\textwidth

     Here the rule is the following: The %LevelIntersectionIterator
     delivers all intersections
     with elements on the same level, the %LeafIntersectionIterator delivers
     the intersections with all leaf elements
     if it has been started on a leaf element.  Both iterators also stop at intersections
     with the grid boundary.
     According to this rule the level intersection iterator started at element \b a
     in the example above delivers an intersection with \b b and the left grid boundary,
     whereas the leaf intersection iterator returns \b c instead of \b b.
     Starting on entity \b c the level intersection iterator returns \b d and the
     intersection with the left boundary of the level 1 grid,
     but the leaf intersection iterator returns both \b d and \b a.
     Finally, starting on \b b the level intersection
     iterator returns \b a and the right boundary, but the leaf intersection iterator is empty since
     \b b is not a leaf entity of the grid. Starting on \b d both the
     level and the leaf intersection iterators will return the element \b c
     together with the right grid boundary.

     @ingroup GIIntersectionIterator
   */
  template<class GridImp, template<class> class IntersectionIteratorImp, template<class> class IntersectionImp>
  class IntersectionIterator
  {
    IntersectionIteratorImp<const GridImp> realIterator;

  public:

    // type of real implementation
    typedef IntersectionIteratorImp<const GridImp> ImplementationType;

    /** \brief Type of Intersection this IntersectionIterator points to */
    typedef Dune::Intersection< const GridImp, IntersectionImp > Intersection;

    //===========================================================
    /** @name Dereferencing
     */
    //@{
    //===========================================================

    /** \brief Dereferencing operator. */
    const Intersection & operator*() const
    {
      return this->realIterator.dereference();
    }

    /** \brief Pointer operator. */
    const Intersection * operator->() const
    {
      return & this->realIterator.dereference();
    }
    //@}


    //===========================================================
    /** @name Compare methods
     */
    //@{
    //===========================================================

    /** @brief Checks for equality.
        Only Iterators pointing to the same intersection from the same Entity
        are equal. Pointing to the same intersection from neighbor is
        unequal as inside and outside are permuted.
     */
    bool operator==(const IntersectionIterator& rhs) const
    {
      return rhs.equals(*this);
    }

    /** @brief Checks for inequality.
        Only Iterators pointing to the same intersection from the same Entity
        are equal. Pointing to the same intersection from neighbor is
        unequal as inside and outside are permuted.
     */
    bool operator!=(const IntersectionIterator& rhs) const
    {
      return ! rhs.equals(*this);
    }
    //@}

    /** @brief Preincrement operator. Proceed to next intersection.*/
    IntersectionIterator& operator++()
    {
      this->realIterator.increment();
      return *this;
    }

    //===========================================================
    /** @name Implementor interface
     */
    //@{
    //===========================================================

    /** @brief forward equality check to realIterator */
    bool equals(const IntersectionIterator& rhs) const
    {
      return this->realIterator.equals(rhs.realIterator);
    }

    /** Copy Constructor from IntersectionIteratorImp */
    IntersectionIterator(const IntersectionIteratorImp<const GridImp> & i) :
      realIterator(i) {};

    /** Copy constructor */
    IntersectionIterator(const IntersectionIterator& i) :
      realIterator(i.realIterator) {}
    //@}

    typedef typename remove_const<GridImp>::type mutableGridImp;
  protected:
    // give the GridDefaultImplementation class access to the realImp
    friend class GridDefaultImplementation<
        GridImp::dimension, GridImp::dimensionworld,
        typename GridImp::ctype,
        typename GridImp::GridFamily> ;

    //! return reference to the real implementation
    ImplementationType & getRealImp() { return realIterator; }
    //! return reference to the real implementation
    const ImplementationType & getRealImp() const { return realIterator; }

  };

  /** \brief specialization of intersection iterator if your grid still follows the old IntersectionIterator semantics

      This class implements a pseudo IntersectionIterator which stores an Intersection<IntersectionIteratorImp> and forwards all calls there.

      For interface documentation please look at Dune::IntersectionIterator

      \deprecated

      \todo clean up this hack (and remove the friend decl in intersection.hh)
   */
  template<class GridImp, template<class> class IntersectionIteratorImp>
  class IntersectionIterator<GridImp, IntersectionIteratorImp, IntersectionIteratorImp>
  {

    Dune::Intersection<const GridImp, IntersectionIteratorImp> intersection;

  public:

    typedef IntersectionIteratorImp<const GridImp> ImplementationType;

    typedef Dune::Intersection<const GridImp, IntersectionIteratorImp> Intersection;

    const Intersection & operator*() const
    {
      return intersection;
    }

    const Intersection * operator->() const
    {
      return & intersection;
    }

    bool operator==(const IntersectionIterator& rhs) const
    {
      return rhs.equals(*this);
    }

    bool operator!=(const IntersectionIterator& rhs) const
    {
      return ! rhs.equals(*this);
    }

    IntersectionIterator& operator++()
    {
      getRealImp().increment();
      return *this;
    }

    bool equals(const IntersectionIterator& rhs) const
    {
      return getRealImp().equals(rhs.getRealImp());
    }

    DUNE_DEPRECATED IntersectionIterator(const IntersectionIteratorImp<const GridImp> & i) :
      intersection(i) {};

    DUNE_DEPRECATED IntersectionIterator(const IntersectionIterator& i) :
      intersection(i.intersection.getRealImp()) {}

    typedef typename remove_const<GridImp>::type mutableGridImp;
  protected:
    // give the GridDefaultImplementation class access to the realImp
    friend class GridDefaultImplementation<
        GridImp::dimension, GridImp::dimensionworld,
        typename GridImp::ctype,
        typename GridImp::GridFamily> ;

    //! return reference to the real implementation
    ImplementationType & getRealImp() { return intersection.getRealImp(); }
    //! return reference to the real implementation
    const ImplementationType & getRealImp() const { return intersection.getRealImp(); }

  };

} // namespace Dune

#include "intersection.hh"

#endif // DUNE_GRID_INTERSECTIONITERATOR_HH
