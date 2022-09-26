// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRID_INTERSECTIONS_HH
#define DUNE_IDENTITYGRID_INTERSECTIONS_HH

#include "identitygridleafiterator.hh"
#include <dune/grid/identitygrid/identitygridentity.hh>

/** \file
 * \brief The IdentityGridLeafIntersection and IdentityGridLevelIntersection classes
 */

namespace Dune {


  // External forward declarations
  template< class Grid >
  struct HostGridAccess;


  /** \brief An intersection with a leaf neighbor element
   * \ingroup IdentityGrid
   * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
   * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
   * These neighbors are accessed via a IntersectionIterator. This allows the implement
   * non-matching meshes. The number of neighbors may be different from the number
   * of an element!
   */
  template<class GridImp>
  class IdentityGridLeafIntersection
  {

    friend class IdentityGridLeafIntersectionIterator<GridImp>;

    friend struct HostGridAccess< typename std::remove_const< GridImp >::type >;

    constexpr static int dim = GridImp::dimension;

    constexpr static int dimworld = GridImp::dimensionworld;

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::HostGridType::LeafGridView::Intersection HostLeafIntersection;

  public:

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef FieldVector<ctype, dimworld> NormalVector;

    IdentityGridLeafIntersection()
    {}

    IdentityGridLeafIntersection(const GridImp* identityGrid,
                                 const HostLeafIntersection& hostIntersection)
      : identityGrid_(identityGrid)
      , hostIntersection_(hostIntersection)
    {}

    IdentityGridLeafIntersection(const GridImp* identityGrid,
                                 HostLeafIntersection&& hostIntersection)
      : identityGrid_(identityGrid)
      , hostIntersection_(std::move(hostIntersection))
    {}

    bool equals(const IdentityGridLeafIntersection& other) const
    {
      return hostIntersection_ == other.hostIntersection_;
    }

    //! return Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    Entity inside() const {
      return IdentityGridEntity<0,dim,GridImp>(identityGrid_,hostIntersection_.inside());
    }


    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    Entity outside() const {
      return IdentityGridEntity<0,dim,GridImp>(identityGrid_,hostIntersection_.outside());
    }


    //! return true if intersection is with boundary.
    bool boundary () const {
      return hostIntersection_.boundary();
    }

    /** \brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal () const {
      return hostIntersection_.centerUnitOuterNormal();
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return hostIntersection_.neighbor();
    }

    //! return the boundary segment index
    size_t boundarySegmentIndex() const {
      return hostIntersection_.boundarySegmentIndex();
    }

    //! Return true if this is a conforming intersection
    bool conforming () const {
      return hostIntersection_.conforming();
    }

    //! Geometry type of an intersection
    GeometryType type () const {
      return hostIntersection_.type();
    }


    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( hostIntersection_.geometryInInside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const
    {
      return LocalGeometry( hostIntersection_.geometryInOutside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const
    {
      return Geometry( hostIntersection_.geometry() );
    }


    //! local number of codim 1 entity in self where intersection is contained in
    int indexInInside () const {
      return hostIntersection_.indexInInside();
    }


    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const {
      return hostIntersection_.indexInOutside();
    }


    //! return outer normal
    FieldVector<ctype, GridImp::dimensionworld> outerNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      return hostIntersection_.outerNormal(local);
    }

    //! return outer normal multiplied by the integration element
    FieldVector<ctype, GridImp::dimensionworld> integrationOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      return hostIntersection_.integrationOuterNormal(local);
    }

    //! return unit outer normal
    FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      return hostIntersection_.unitOuterNormal(local);
    }


  private:
    //**********************************************************
    //  private methods
    //**********************************************************

    const GridImp* identityGrid_;

    HostLeafIntersection hostIntersection_;
  };




  //! \todo Please doc me !
  template<class GridImp>
  class IdentityGridLevelIntersection
  {

    friend class IdentityGridLevelIntersectionIterator<GridImp>;

    friend struct HostGridAccess< typename std::remove_const< GridImp >::type >;

    constexpr static int dim = GridImp::dimension;

    constexpr static int dimworld = GridImp::dimensionworld;

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::HostGridType::LevelGridView::Intersection HostLevelIntersection;

  public:

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef FieldVector<ctype, dimworld> NormalVector;

    IdentityGridLevelIntersection()
    {}

    IdentityGridLevelIntersection(const GridImp* identityGrid,
                                  const HostLevelIntersection& hostIntersection)
      : identityGrid_(identityGrid)
      , hostIntersection_(hostIntersection)
    {}

    IdentityGridLevelIntersection(const GridImp* identityGrid,
                                  HostLevelIntersection&& hostIntersection)
      : identityGrid_(identityGrid)
      , hostIntersection_(std::move(hostIntersection))
    {}

    bool equals(const IdentityGridLevelIntersection& other) const
    {
      return hostIntersection_ == other.hostIntersection_;
    }

    //! return Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    Entity inside() const {
      return IdentityGridEntity<0,dim,GridImp>(identityGrid_,hostIntersection_.inside());
    }


    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    Entity outside() const {
      return IdentityGridEntity<0,dim,GridImp>(identityGrid_,hostIntersection_.outside());
    }


    /** \brief return true if intersection is with boundary.
     */
    bool boundary () const {
      return hostIntersection_.boundary();
    }

    /** \brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal () const {
      return hostIntersection_.centerUnitOuterNormal();
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return hostIntersection_.neighbor();
    }

    //! return the boundary segment index
    size_t boundarySegmentIndex() const {
      return hostIntersection_.boundarySegmentIndex();
    }

    //! Return true if this is a conforming intersection
    bool conforming () const {
      return hostIntersection_.conforming();
    }

    //! Geometry type of an intersection
    GeometryType type () const {
      return hostIntersection_.type();
    }


    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( hostIntersection_.geometryInInside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const
    {
      return LocalGeometry( hostIntersection_.geometryInOutside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const
    {
      return Geometry( hostIntersection_.geometry() );
    }


    //! local number of codim 1 entity in self where intersection is contained in
    int indexInInside () const {
      return hostIntersection_.indexInInside();
    }


    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const {
      return hostIntersection_.indexInOutside();
    }


    //! return outer normal
    FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dim-1>& local) const {
      return hostIntersection_.outerNormal(local);
    }

    //! return outer normal multiplied by the integration element
    FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const {
      return hostIntersection_.integrationOuterNormal(local);
    }

    //! return unit outer normal
    FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dim-1>& local) const {
      return hostIntersection_.unitOuterNormal(local);
    }

  private:

    const GridImp* identityGrid_;

    HostLevelIntersection hostIntersection_;

  };


}  // namespace Dune

#endif
