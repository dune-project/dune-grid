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

    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };

  public:

    // type of real implementation
    typedef IntersectionIteratorImp<const GridImp> ImplementationType;

    /** \brief Type of entity that this IntersectionIterator belongs to */
    typedef typename GridImp::template Codim<0>::Entity Entity;

    /** \brief Pointer to the type of entities that this IntersectionIterator belongs to */
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    /** \brief Type of Intersection this IntersectionIterator points to */
    typedef Dune::Intersection< const GridImp, IntersectionImp > Intersection;

    /** \brief Codim 1 geometry returned by intersectionGlobal() */
    typedef typename GridImp::template Codim<1>::Geometry Geometry;

    /** \brief Codim 1 geometry returned by intersectionSelfLocal()
        and intersectionNeighborLocal() */
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    //! @brief Export grid dimension
    enum { dimension=dim /*!< grid dimension */ };

    //! @brief Export dimension of world
    enum { dimensionworld=dimworld /*!< dimension of world */ };

    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

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
    /** @name Query methods
     */
    //@{
    //===========================================================

    //! return true if intersection is with interior or exterior boundary (see the cases above)
    bool boundary () const DUNE_DEPRECATED
    {
      return (*this)->boundary();
    }

    /**
       \brief Identifier for boundary segment from macro grid.

       One can attach a boundary Id to a boundary segment on the macro
       grid. This Id will also be used for all fragments of these
       boundary segments.

       The numbering is defined as:
       - Id==0 for all intersections without boundary()==false
       - Id>=0 for all intersections without boundary()==true

       The way the Identifiers are attached to the grid may differ
       between the different grid implementations.
     */
    int boundaryId () const DUNE_DEPRECATED
    {
      return (*this)->boundaryId();
    }

    //! @brief return true if intersection is shared with another element.
    bool neighbor () const DUNE_DEPRECATED
    {
      return (*this)->neighbor();
    }

    /*! @brief return EntityPointer to the Entity on the inside of this
       intersection. That is the Entity where we started this Iterator.
     */
    EntityPointer inside() const DUNE_DEPRECATED
    {
      return (*this)->inside();
    }

    /*! @brief return EntityPointer to the Entity on the outside of this
       intersection. That is the neighboring Entity.

       @warning Don't call this method if there is no neighboring Entity
       (neighbor() returns false). In this case the result is undefined.
     */
    EntityPointer outside() const DUNE_DEPRECATED
    {
      return (*this)->outside();
    }

    /*! @brief geometrical information about this intersection in local
       coordinates of the inside() entity.

       This method returns a Geometry object that provides a mapping from
       local coordinates of the intersection to local coordinates of the
       inside() entity.
     */
    const LocalGeometry& intersectionSelfLocal () const DUNE_DEPRECATED
    {
      return (*this)->intersectionSelfLocal();
    }
    /*! @brief geometrical information about this intersection in local
       coordinates of the outside() entity.

       This method returns a Geometry object that provides a mapping from
       local coordinates of the intersection to local coordinates of the
       outside() entity.
     */
    const LocalGeometry& intersectionNeighborLocal () const DUNE_DEPRECATED
    {
      return (*this)->intersectionNeighborLocal();
    }

    /*! @brief geometrical information about this intersection in global coordinates.

       This method returns a Geometry object that provides a mapping from
       local coordinates of the intersection to global (world) coordinates.
     */
    const Geometry& intersectionGlobal () const DUNE_DEPRECATED
    {
      return (*this)->intersectionGlobal();
    }

    //! Local number of codim 1 entity in the inside() Entity where intersection is contained in
    int numberInSelf () const DUNE_DEPRECATED
    {
      return (*this)->numberInSelf();
    }

    //! Local number of codim 1 entity in outside() Entity where intersection is contained in
    int numberInNeighbor () const DUNE_DEPRECATED
    {
      return (*this)->numberInNeighbor();
    }

    /*! @brief Return an outer normal (length not necessarily 1)

       The returned vector may depend on local position within the intersection.
     */
    FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dim-1>& local) const DUNE_DEPRECATED
    {
      return (*this)->outerNormal( local );
    }

    /*! @brief return outer normal scaled with the integration element
          @copydoc outerNormal
       The normal is scaled with the integration element of the intersection. This
          method is redundant but it may be more efficent to use this function
          rather than computing the integration element via intersectionGlobal().
     */
    FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const DUNE_DEPRECATED
    {
      return (*this)->integrationOuterNormal( local );
    }

    /*! @brief Return unit outer normal (length == 1)

       The returned vector may depend on the local position within the intersection.
       It is scaled to have unit length.
     */
    FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dim-1>& local) const DUNE_DEPRECATED
    {
      return (*this)->unitOuterNormal( local );
    }
    //@}

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

} // namespace Dune

#include "intersection.hh"

#endif // DUNE_GRID_INTERSECTIONITERATOR_HH
