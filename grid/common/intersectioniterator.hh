// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_INTERSECTIONITERATOR_HH
#define DUNE_GRID_INTERSECTIONITERATOR_HH

#include <dune/common/iteratorfacades.hh>

namespace Dune
{

  /** \brief Mesh entities of codimension 0 ("elements") allow to visit all
     intersections with "neighboring" elements and with the domain
     boundary.

     Template parameters are:

     - <tt>GridImp</tt> Type that is a model of Dune::Grid
     - <tt>IntersectionIteratorImp</tt> Class template that is a model of
     Dune::IntersectionIterator

     @warning The name IntersectionIterator may be somewhat misleading. This
     class has neither an operator* nor an operator->. It iterates over codimension
     1 intersections with other entities and the (sub-)domain boundary.

     @warning The number of neigbors may be different from the number of
     faces/edges of an element!

     <h2>Overview</h2>

     Intersections are codimension 1 objects. These
     intersections are accessed via a IntersectionIterator. This allows
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
     IntersectionIterators by either using ilevelbegin()/ilevelend() or
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
     with elements on the same level, the %LeafIntersectionIterator
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

     <h2>Methods neighbor and boundary </h2>

     The following holds for both the level and the leaf intersection
     iterator:
     The %intersection iterator is started on a codimension 0 entity of the grid.
     If this entity belongs to the interior or the overlap region
     (see. ???) then the union of all intersections is identical to the
     boundary of the entity. On ghost elements the iterator only stops
     on the border of the domain, i.e., only on the intersections with
     entities in the interior region. Depending on the boolean values returned by
     the methods %boundary() and %neighbor()
     one can detect the position of an intersection
     relative to the boundary. In general
     the method boundary() returns true if and only if the intersection is
     part of the physical boundary of the domain. The method neighbor() returns
     true only if the method outside() has a well defined return value.

     The following cases are possible if the intersection iterator is
     started on an entity in the interior or overlap region. More
     details are given below:

     <table>
     <tr>
     <td></td><td>intersection</td><td>neighbor()</td><td>boundary()</td><td>outside()</td></tr>
     <tr>
     <td>1</td><td>with inner, overlap <br>
                  or ghost entity</td>
     <td>true</td><td>false</td>
     <td>the neighbor entity</td></tr>
     <tr>
     <td>2</td><td>on domain boundary</td>
     <td>false</td><td>true</td><td><em>undefined</em></td></tr>
     <tr>
     <td>3</td><td>on periodic boundary</td>
     <td>true</td><td>true</td><td>Ghost-/Overlap cell@br (with transformed geometry)</td></tr>
     <tr>
     <td>4</td><td>on processor boundary</td>
     <td>false <em>if grid has no ghosts</em><br>true <em>otherwise</em></td><td>false </td>
     <td>ghost entity <em>(if it exists)</em></td></tr>
     </table>

     -# <b> Inner Intersections: </b> \n
       The type of the neighboring entity can be determined through
       methods defined on the outside entity.
     -# <b>  Handling physical boundaries: </b>
       Different types of physical boundaries can be modeled using either
       the global coordinates of the intersection or by using the
       boundaryID method. On some grids (AluGrid, AlbertaGrid) this
       method returns an integer value which can be individually assigned
       to each boundary intersection of the macro grid and which is
       prolonged to higher levels during grid refinement. \br
       A more general concept will be included in latter releases along the
       following guidelines:
       - We require differently constructed geometries outside the domain
       - The kind of construction depends on the discrete problem
       - Therefor these constructions can't be part of the Grid interface
       - Utility classes are required to do this construction
       - The utility classes must be parameterized with the intersection (in our
         case the IntersectionIterator)
       - The utility classes return a suitable transformation of the inner()
         entitys geometry (with respect to the intersection), e.g.,
         reflection at the intersection
         point reflection
         reflection combined with translation...
       .
     -# <b> Handling periodic boundaries: </b>
       - The IntersectionIterator stops at periodic boundaries
       - periodic grids are handled in correspondence to parallel grids
       - %At the periodic boundary one can adjust an overlap- or ghost-layer.
       - outer() returns a ghost or overlap cell (for ghost and overlap look into
         the documentation of the parallel grid interface)
       - outer() cell has a periodically transformed geometry (so that one does
         not see a jump or something like that)
       - outer() cell has its own index
       - outer() cell has the same id as the corresponding "original" cell
     -# <b> Processor boundaries: </b> \n
       At processor boundaries, i.e. when an element has an intersection with
       another element
       in the sequential grid but this element is only stored in other processors
       the intersection iterator stops but neither
       neighbor() nor boundary()
       are true.
     .


     <h2>Geometry of an intersection</h2>

     The method intersectionGlobal returns a geometry mapping the intersection
     as a codim one structure to global coordinates. The methods
     intersectionSelfLocal and intersectionNeighborLocal return geometries
     mapping the intersection into the reference elements of the
     originating entity and the neighboring entity, respectively.
     The numberInSelf and numberInNeighbor methods return the codim one
     subentities which contain the intersection.


     @ingroup GIIntersectionIterator
   */
  template<class GridImp, template<class> class IntersectionIteratorImp>
  class IntersectionIterator
  {
    IntersectionIteratorImp<const GridImp> realIterator;

    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };
    typedef typename GridImp::ctype ct;
  public:

    // type of real implementation
    typedef IntersectionIteratorImp<const GridImp> ImplementationType;

    /** \brief Type of entity that this IntersectionIterator belongs to */
    typedef typename GridImp::template Codim<0>::Entity Entity;

    /** \brief Pointer to the type of entities that this IntersectionIterator belongs to */
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

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
    typedef ct ctype;

    /** @brief Preincrement operator. Proceed to next intersection.*/
    IntersectionIterator& operator++()
    {
      this->realIterator.increment();
      return *this;
    }

    //! return true if intersection is with interior or exterior boundary (see the cases above)
    bool boundary () const
    {
      return this->realIterator.boundary();
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
    int boundaryId () const
    {
      return this->realIterator.boundaryId();
    }

    //! @brief return true if intersection is shared with another element.
    bool neighbor () const
    {
      return this->realIterator.neighbor();
    }

    /*! @brief return EntityPointer to the Entity on the inside of this
       intersection. That is the Entity where we started this Iterator.
     */
    EntityPointer inside() const
    {
      return this->realIterator.inside();
    }

    /*! @brief return EntityPointer to the Entity on the outside of this
       intersection. That is the neighboring Entity.

       @warning Don't call this method if there is no neighboring Entity
       (neighbor() returns false). In this case the result is undefined.
     */
    EntityPointer outside() const
    {
      return this->realIterator.outside();
    }

    /*! @brief geometrical information about this intersection in local
       coordinates of the inside() entity.
       This method returns a Geometry object that provides a mapping from
       local coordinates of the intersection to local coordinates of the
       inside() entity.
     */
    const LocalGeometry& intersectionSelfLocal () const
    {
      return this->realIterator.intersectionSelfLocal();
    }
    /*! @brief geometrical information about this intersection in local
       coordinates of the outside() entity.
       This method returns a Geometry object that provides a mapping from
       local coordinates of the intersection to local coordinates of the
       outside() entity.
     */
    const LocalGeometry& intersectionNeighborLocal () const
    {
      return this->realIterator.intersectionNeighborLocal();
    }

    /*! @brief geometrical information about this intersection in global coordinates.
       This method returns a Geometry object that provides a mapping from
       local coordinates of the intersection to global (world) coordinates.
     */
    const Geometry& intersectionGlobal () const
    {
      return this->realIterator.intersectionGlobal();
    }

    //! Local number of codim 1 entity in the inside() Entity where intersection is contained in
    int numberInSelf () const
    {
      return this->realIterator.numberInSelf ();
    }

    //! Local number of codim 1 entity in outside() Entity where intersection is contained in
    int numberInNeighbor () const
    {
      return this->realIterator.numberInNeighbor ();
    }

    /*! @brief Return an outer normal (length not necessarily 1)
       The returned Vector is copied to take advantage from the return
       type optimization.  Usually one will use the normal in several
       calculations, so that one stores it before using it. The
       returned vector may depend on local position within the intersection.
     */
    FieldVector<ct, dimworld> outerNormal (const FieldVector<ct, dim-1>& local) const
    {
      return this->realIterator.outerNormal(local);
    }

    /*! @brief return outer normal scaled with the integration element
          @copydoc outerNormal
       The normal is scaled with the integration element of the intersection. This
          method is redundant but it may be more efficent to use this function
          rather than computing the integration element via intersectionGlobal.
     */
    FieldVector<ct, dimworld> integrationOuterNormal (const FieldVector<ct, dim-1>& local) const
    {
      return this->realIterator.integrationOuterNormal(local);
    }

    /*! @brief Return unit outer normal (length == 1)

       The returned vector may depend on the local position within the intersection.
       It is scaled to have unit length.
     */
    FieldVector<ct, dimworld> unitOuterNormal (const FieldVector<ct, dim-1>& local) const
    {
      return this->realIterator.unitOuterNormal(local);
    }

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
            @copydoc operator==
     */
    bool operator!=(const IntersectionIterator& rhs) const
    {
      return ! rhs.equals(*this);
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

    typedef typename RemoveConst<GridImp>::Type mutableGridImp;
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

    // private, so that nobody will use this method
    /** @brief Postincrement operator. Deprecated, do not use it anymore.*/
    IntersectionIterator operator++(int) DUNE_DEPRECATED
    {
      const IntersectionIterator tmp(*this);
      this->realIterator.increment();
      return tmp;
    }

  };

  //**********************************************************************
  /**
     @brief Default Implementations for IntersectionIteratorImp

     @ingroup GridDevel
   */
  template<class GridImp, template<class> class IntersectionIteratorImp>
  class IntersectionIteratorDefaultImplementation
  {
    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };
    typedef typename GridImp::ctype ct;
  public:
    //! return unit outer normal, this should be dependent on
    //! local coordinates for higher order boundary
    //! the normal is scaled with the integration element of the intersection.
    FieldVector<ct, dimworld> integrationOuterNormal (const FieldVector<ct, dim-1>& local) const
    {
      FieldVector<ct, dimworld> n = unitOuterNormal(local);
      n *= asImp().intersectionGlobal().integrationElement(local);
      return n;
    }
    //! return unit outer normal
    FieldVector<ct, dimworld> unitOuterNormal (const FieldVector<ct, dim-1>& local) const
    {
      FieldVector<ct, dimworld> n = asImp().outerNormal(local);
      n /= n.two_norm();
      return n;
    }

  private:
    //!  Barton-Nackman trick
    IntersectionIteratorImp<GridImp>& asImp ()
    {return static_cast<IntersectionIteratorImp<GridImp>&>(*this);}
    const IntersectionIteratorImp<GridImp>& asImp () const
    {return static_cast<const IntersectionIteratorImp<GridImp>&>(*this);}
  };

}

#endif // DUNE_GRID_INTERSECTIONITERATOR_HH
