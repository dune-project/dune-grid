// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_COMMON_INTERSECTION_HH
#define DUNE_GRID_COMMON_INTERSECTION_HH

#include <dune/grid/common/grid.hh>

namespace Dune
{

  /** \brief %Intersection of a mesh entity of codimension 0 ("element")
      with a "neighboring" element or with the domain
      boundary.

     \tparam GridImp Type that is a model of Dune::Grid
     \tparam IntersectionImp Class template that is a model of Dune::Intersection

     <h2>Overview</h2>

     Intersections are codimension 1 objects. These
     intersections are accessed via an Intersection. This allows
     the implementation of non-matching grids, as one facet can now
     consist of several intersections.
     In a conforming mesh such an intersection corresponds to an entity of
     codimension 1 but in the general non-conforming case there will be no entity
     in the mesh that directly corresponds to the intersection. Thus, the
     Intersection describes these intersections implicitly.

     <H2>Engine Concept</H2>

     The Intersection class template wraps an object of type IntersectionImp
     and forwards all member
     function calls to corresponding members of this class. In that sense Intersection
     defines the interface and IntersectionImp supplies the implementation.

     <h2>Methods neighbor and boundary </h2>

     The following holds for both the level and the leaf intersection
     :
     The %intersection  is started on a codimension 0 entity of the grid.
     If this entity belongs to the interior or the overlap region
     (see. ???) then the union of all intersections is identical to the
     boundary of the entity. On ghost elements the  only stops
     on the border of the domain, i.e., only on the intersections with
     entities in the interior region. Depending on the boolean values returned by
     the methods %boundary() and %neighbor()
     one can detect the position of an intersection
     relative to the boundary. In general
     the method boundary() returns true if and only if the intersection is
     part of the physical boundary of the domain. The method neighbor() returns
     true only if the method outside() has a well defined return value.

     The following cases are possible if the intersection  is
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
     <td>true</td><td>true</td><td>Ghost-/Overlap cell <br>
                                   (with transformed geometry)</td></tr>
     <tr>
     <td>4</td><td>on processor boundary</td>
     <td>false <em>if grid has no ghosts</em><br>true <em>otherwise</em></td><td>false </td>
     <td>ghost entity <em>(if it exists)</em></td></tr>
     </table>

     <dl>
     <dt>Inner Intersections: </dt>
     <dd>
       The type of the neighboring entity can be determined through
       methods defined on the outside entity.
     </dd>
     <dt>Handling physical boundaries: </dt>
     <dd>
       Data (like the type of boundary) can be attached to physical boundaries
       either using global coordinates or the intersection's boundary segment
       index.<br>
       The boundary segment indices form a local, zero-based, contiguous set of
       integer values.
       Each boundary segment on the macro level is assigned a unique index from
       this set, which is inherited by child boundary segments.
       The size of the boundary segment index set (i.e., the number of boundary
       indices on the macro level) can be determined through the method
       Grid::numBoundarySegments.<br>
       Note that the boundary segment index is not persistent over dynamic load
       balancing.
     \if old_documentation
       Different types of physical boundaries can be modeled using either
       the global coordinates of the intersection or by using the
       boundaryID method. On some grids (dune-ALUGrid, AlbertaGrid) this
       method returns an integer value which can be individually assigned
       to each boundary intersection of the macro grid and which is
       prolonged to higher levels during grid refinement.<br>
       A more general concept will be included in latter releases along the
       following guidelines:
       - We require differently constructed geometries outside the domain
       - The kind of construction depends on the discrete problem
       - Therefor these constructions can't be part of the Grid interface
       - Utility classes are required to do this construction
       - The utility classes must be parameterized with the intersection (in our
         case the Intersection)
       - The utility classes return a suitable transformation of the inner()
         entitys geometry (with respect to the intersection), e.g.,
         reflection at the intersection
         point reflection
         reflection combined with translation...
       .
     \endif
     </dd>
     <dt>Handling periodic boundaries: </dt>
     <dd>
       - The Intersection stops at periodic boundaries
       - periodic grids are handled in correspondence to parallel grids
       - %At the periodic boundary one can adjust an overlap- or ghost-layer.
       - outer() returns a ghost or overlap cell (for ghost and overlap look into
         the documentation of the parallel grid interface)
       - outer() cell has a periodically transformed geometry (so that one does
         not see a jump or something like that)
       - outer() cell has its own index
       - outer() cell has the same id as the corresponding "original" cell
       .
     </dd>
     <dt>Processor boundaries: </dt>
     <dd>
       At processor boundaries, i.e. when an element has an intersection with
       another element
       in the sequential grid but this element is only stored in other processors
       the intersection  stops but neither
       neighbor() nor boundary()
       are true.
     </dd>
     </dl>


     <h2>Geometry of an intersection</h2>

     The method geometry returns a geometry mapping the intersection
     as a codim one structure to global coordinates. The methods
     geometryInInside and geometryInOutside return geometries
     mapping the intersection into the reference elements of the
     originating entity and the neighboring entity, respectively.
     The indexInInside and indexInOutside methods return the codim one
     subentities which contain the intersection.


     @ingroup GIIntersectionIterator
   */
  template< class GridImp, class IntersectionImp >
  class Intersection
  {
  public:
    /**
     * \brief type of underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    typedef IntersectionImp Implementation;

    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    Implementation &impl () { return real; }

    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    const Implementation &impl () const { return real; }

  protected:
    Implementation real;

  public:
    /** \brief Type of entity that this Intersection belongs to */
    typedef typename GridImp::template Codim<0>::Entity Entity;

    /** \brief Codim 1 geometry returned by geometry() */
    typedef typename GridImp::template Codim<1>::Geometry Geometry;

    /** \brief Type for vectors of coordinates on the intersection */
    typedef typename Geometry::LocalCoordinate LocalCoordinate;

    /** \brief Type for normal vectors */
    typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

    /** \brief Codim 1 geometry returned by geometryInInside() and geometryInOutside() */
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    //! @brief dimension of the intersection
    constexpr static int mydimension = GridImp::dimension - 1;

    //! @brief dimension of world
    constexpr static int dimensionworld = GridImp::dimensionworld;

    //! Type of individual coefficients of coordinate vectors
    typedef typename GridImp::ctype ctype;

    //! Return true if intersection is with interior or exterior boundary (see the cases above)
    bool boundary () const
    {
      return this->real.boundary();
    }

    /** \brief index of the boundary segment within the macro grid
     *
     *  In many applications, special data needs to be attached to the boundary
     *  segments of the macro grid (e.g., a function selecting the boundary
     *  condition).
     *  Usually, this data is inherited by the children of the boundary segment.
     *
     *  In the DUNE framework, data is stored in arrays, addressed by an index,
     *  in this case the boundarySegmentIndex. The size of these arrays can be
     *  obtained by the Grid::numBoundarySegments method.
     *
     *  The indices returned by this method are consecutive, zero based, and local to the
     *  processor. Notice that these indices do not necessarily coincide with the insertion
     *  indices of the corresponding boundary segments as provided by the GridFactory.
     */
    size_t boundarySegmentIndex () const
    {
      return this->real.boundarySegmentIndex();
    }

    //! @brief return true if intersection is shared with another element.
    bool neighbor () const
    {
      return this->real.neighbor();
    }

    /*! @brief return Entity on the inside of this
       intersection. That is the Entity where we started this.
     */
    Entity inside() const
    {
      return this->real.inside();
    }

    /*! @brief return Entity on the outside of this
       intersection. That is the neighboring Entity.

       @warning Don't call this method if there is no neighboring Entity
       (neighbor() returns false). In this case the result is undefined.
     */
    Entity outside() const
    {
      return this->real.outside();
    }

    /*! @brief Return true if intersection is conforming.
     */
    bool conforming () const
    {
      return this->real.conforming();
    }

    /** \brief geometrical information about this intersection in local
     *         coordinates of the inside() entity.
     *
     *  This method returns a Geometry object that provides a mapping from
     *  local coordinates of the intersection to local coordinates of the
     *  inside() entity.
     *
     *  \note Previously, the geometry was encapsulated in the intersection object
     *        and a const reference was returned.
     *
     *  \note The returned geometry object is guaranteed to remain valid until the
     *        grid is modified (or deleted).
     */
    LocalGeometry geometryInInside () const
    {
      return this->real.geometryInInside();
    }

    /** \brief geometrical information about this intersection in local
     *         coordinates of the outside() entity.
     *
     *  This method returns a Geometry object that provides a mapping from
     *  local coordinates of the intersection to local coordinates of the
     *  outside() entity.
     *
     *  \note Previously, the geometry was encapsulated in the intersection object
     *        and a const reference was returned.
     *
     *  \note The returned geometry object is guaranteed to remain valid until the
     *        grid is modified (or deleted).
     */
    LocalGeometry geometryInOutside () const
    {
      return this->real.geometryInOutside();
    }

    /** \brief geometrical information about the intersection in global coordinates.
     *
     *  This method returns a Geometry object that provides a mapping from
     *  local coordinates of the intersection to global (world) coordinates.
     *
     *  \note If the returned geometry has type <b>none</b> then only a limited set of features
     *        is available for the geometry, i.e. center and volume.
     *
     *  \note Previously, the geometry was encapsulated in the intersection object
     *        and a const reference was returned.
     *
     *  \note The returned geometry object is guaranteed to remain valid until the
     *        grid is modified (or deleted).
     */
    Geometry geometry () const
    {
      return this->real.geometry();
    }

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return this->real.type();
    }

    /** \brief Local index of codim 1 entity in the inside() entity where
     *         intersection is contained in
     *
     *  This method returns the facet number with respect to the generic reference
     *  element.
     *
     *  \note This index can be used with the inside entity's
     *        `subEntity<1>()` method to obtain the facet.
     *
     *  \returns the index of the inside entity's facet containing this
     *           intersection (with respect to the generic reference element)
     */
    int indexInInside () const
    {
      return this->real.indexInInside();
    }

    /** \brief Local index of codim 1 entity in outside() entity where
     *         intersection is contained in
     *
     *  This method returns the facet number with respect to the generic reference
     *  element.
     *
     *  \note This index can be used with the outside entity's
     *        `subEntity<1>()` method to obtain the facet.
     *
     *  \returns the index of the outside entity's facet containing this
     *           intersection (with respect to the generic reference element)
     */
    int indexInOutside () const
    {
      return this->real.indexInOutside();
    }

    /*! @brief Return an outer normal (length not necessarily 1)

       The returned vector may depend on local position within the intersection.
     */
    GlobalCoordinate outerNormal (const LocalCoordinate& local) const
    {
      return this->real.outerNormal(local);
    }

    /*! @brief return unit outer normal scaled with the integration element

       The normal is scaled with the integration element of the intersection. This
          method is redundant but it may be more efficient to use this function
          rather than computing the integration element via geometry().

       The returned vector may depend on local position within the intersection.
     */
    GlobalCoordinate integrationOuterNormal (const LocalCoordinate& local) const
    {
      return this->real.integrationOuterNormal(local);
    }

    /*! @brief Return unit outer normal (length == 1)

       The returned vector may depend on the local position within the intersection.
       It is scaled to have unit length.
     */
    GlobalCoordinate unitOuterNormal (const LocalCoordinate& local) const
    {
      return this->real.unitOuterNormal(local);
    }

    /*! @brief Return unit outer normal (length == 1)

       The returned vector is the normal at the center() of the
       intersection's geometry.
       It is scaled to have unit length.
     */
    GlobalCoordinate centerUnitOuterNormal () const
    {
      return this->real.centerUnitOuterNormal();
    }

    //! Compares two intersections for equality.
    bool operator==(const Intersection& other) const
    {
      return real.equals(other.real);
    }

    //! Compares two intersections for inequality.
    bool operator!=(const Intersection& other) const
    {
      return !real.equals(other.real);
    }

    //! Default constructor.
    Intersection()
    {}

    //! Copy constructor from an existing intersection.
    Intersection(const Intersection& other)
      : real(other.real)
    {}

    //! Move constructor from an existing intersection.
    Intersection(Intersection&& other)
      : real(std::move(other.real))
    {}

    //! Copy assignment operator from an existing intersection.
    Intersection& operator=(const Intersection& other)
    {
      real = other.real;
      return *this;
    }

    //! Move assignment operator from an existing intersection.
    Intersection& operator=(Intersection&& other)
    {
      real = std::move(other.real);
      return *this;
    }

    //===========================================================
    /** @name Implementor interface
     */
    //@{
    //===========================================================

    /** Copy Constructor from IntersectionImp */
    Intersection ( const Implementation &impl )
      : real( impl )
    {}

    /** Move Constructor from IntersectionImp */
    Intersection ( Implementation&& impl )
      : real( std::move(impl) )
    {}

    //@}

  protected:
    //! give the pseudo IntersectionIterator class access to the realImp
    //! \todo cleanup this hack
    friend class IntersectionIterator<GridImp, IntersectionImp, IntersectionImp>;

  };

  //**********************************************************************
  /**
     @brief Default Implementations of integrationOuterNormal and unitOuterNormal for IntersectionImp

     @ingroup GridDevel
   */
  template< class GridImp, class IntersectionImp >
  class IntersectionDefaultNormalVectors
  {
    constexpr static int dim = GridImp::dimension;
    constexpr static int dimworld = GridImp::dimensionworld;
    typedef typename GridImp::ctype ct;
  public:

    //! return unit outer normal, this should be dependent on
    //! local coordinates for higher order boundary
    //! the normal is scaled with the integration element of the intersection.
    FieldVector<ct, dimworld> integrationOuterNormal (const FieldVector<ct, dim-1>& local) const
    {
      FieldVector<ct, dimworld> n = asImp().unitOuterNormal(local);
      n *= asImp().geometry().integrationElement(local);
      return n;
    }

    //! return unit outer normal
    FieldVector<ct, dimworld> unitOuterNormal (const FieldVector<ct, dim-1>& local) const
    {
      FieldVector<ct, dimworld> n = asImp().outerNormal(local);
      n /= n.two_norm();
      return n;
    }

    //! return unit outer normal at center of intersection geometry
    FieldVector<ct, dimworld> centerUnitOuterNormal () const
    {
      // For now, we do this...
      GeometryType type = asImp().geometry().type();
      auto refElement = referenceElement<ct, dim-1>(type);
      return asImp().unitOuterNormal(refElement.position(0,0));
      // But later, if we change the meaning of center(),
      // we may have to change to this...
      // return asImp().unitOuterNormal(asImp().local(asImp().center()));
    }

  private:
    // CRTP (curiously recurring template pattern)
    IntersectionImp &asImp () { return static_cast< IntersectionImp & >( *this ); }
    const IntersectionImp &asImp () const { return static_cast< const IntersectionImp & >( *this ); }
  };

} // namespace Dune

#endif // DUNE_GRID_COMMON_INTERSECTION_HH
