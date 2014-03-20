// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRID_INTERSECTIONS_HH
#define DUNE_IDENTITYGRID_INTERSECTIONS_HH

#include "identitygridleafiterator.hh"

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

    friend struct HostGridAccess< typename remove_const< GridImp >::type >;

    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::HostGridType::template Codim<0>::Entity::LeafIntersectionIterator HostLeafIntersectionIterator;

  public:

    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef FieldVector<ctype, dimworld> NormalVector;

    IdentityGridLeafIntersection(const GridImp* identityGrid,
                                 const HostLeafIntersectionIterator& hostIterator)
      : identityGrid_(identityGrid),
        hostIterator_(hostIterator)
    {}

    //! return EntityPointer to the Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    EntityPointer inside() const {
      return IdentityGridEntityPointer<0,GridImp> (identityGrid_, hostIterator_->inside());
    }


    //! return EntityPointer to the Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    EntityPointer outside() const {
      return IdentityGridEntityPointer<0,GridImp> (identityGrid_, hostIterator_->outside());
    }


    //! return true if intersection is with boundary.
    bool boundary () const {
      return hostIterator_->boundary();
    }

    /** \brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal () const {
      return hostIterator_->centerUnitOuterNormal();
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return hostIterator_->neighbor();
    }


    //! return information about the Boundary
    int boundaryId () const {
      return hostIterator_->boundaryId();
    }

    //! return the boundary segment index
    size_t boundarySegmentIndex() const {
      return hostIterator_->boundarySegmentIndex();
    }

    //! Return true if this is a conforming intersection
    bool conforming () const {
      return hostIterator_->conforming();
    }

    //! Geometry type of an intersection
    GeometryType type () const {
      return hostIterator_->type();
    }


    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( hostIterator_->geometryInInside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const
    {
      return LocalGeometry( hostIterator_->geometryInOutside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const
    {
      return Geometry( hostIterator_->geometry() );
    }


    //! local number of codim 1 entity in self where intersection is contained in
    int indexInInside () const {
      return hostIterator_->indexInInside();
    }


    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const {
      return hostIterator_->indexInOutside();
    }


    //! return outer normal
    FieldVector<ctype, GridImp::dimensionworld> outerNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      return hostIterator_->outerNormal(local);
    }

    //! return outer normal multiplied by the integration element
    FieldVector<ctype, GridImp::dimensionworld> integrationOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      return hostIterator_->integrationOuterNormal(local);
    }

    //! return unit outer normal
    FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      return hostIterator_->unitOuterNormal(local);
    }


  private:
    //**********************************************************
    //  private methods
    //**********************************************************

    const GridImp* identityGrid_;

    HostLeafIntersectionIterator hostIterator_;
  };




  //! \todo Please doc me !
  template<class GridImp>
  class IdentityGridLevelIntersection
  {

    friend class IdentityGridLevelIntersectionIterator<GridImp>;

    friend struct HostGridAccess< typename remove_const< GridImp >::type >;

    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::HostGridType::template Codim<0>::Entity::LevelIntersectionIterator HostLevelIntersectionIterator;

  public:

    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef FieldVector<ctype, dimworld> NormalVector;

    IdentityGridLevelIntersection(const GridImp* identityGrid,
                                  const HostLevelIntersectionIterator& hostIterator)
      : identityGrid_(identityGrid), hostIterator_(hostIterator)
    {}

    //! return EntityPointer to the Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    EntityPointer inside() const {
      return IdentityGridEntityPointer<0,GridImp> (identityGrid_, hostIterator_->inside());
    }


    //! return EntityPointer to the Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    EntityPointer outside() const {
      return IdentityGridEntityPointer<0,GridImp> (identityGrid_, hostIterator_->outside());
    }


    /** \brief return true if intersection is with boundary.
     */
    bool boundary () const {
      return hostIterator_->boundary();
    }

    /** \brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal () const {
      return hostIterator_->centerUnitOuterNormal();
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return hostIterator_->neighbor();
    }


    //! return information about the Boundary
    int boundaryId () const {
      return hostIterator_->boundaryId();
    }

    //! return the boundary segment index
    size_t boundarySegmentIndex() const {
      return hostIterator_->boundarySegmentIndex();
    }

    //! Return true if this is a conforming intersection
    bool conforming () const {
      return hostIterator_->conforming();
    }

    //! Geometry type of an intersection
    GeometryType type () const {
      return hostIterator_->type();
    }


    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( hostIterator_->geometryInInside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const
    {
      return LocalGeometry( hostIterator_->geometryInOutside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const
    {
      return Geometry( hostIterator_->geometry() );
    }


    //! local number of codim 1 entity in self where intersection is contained in
    int indexInInside () const {
      return hostIterator_->indexInInside();
    }


    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const {
      return hostIterator_->indexInOutside();
    }


    //! return outer normal
    FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dim-1>& local) const {
      return hostIterator_->outerNormal(local);
    }

    //! return outer normal multiplied by the integration element
    FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const {
      return hostIterator_->integrationOuterNormal(local);
    }

    //! return unit outer normal
    FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dim-1>& local) const {
      return hostIterator_->unitOuterNormal(local);
    }

  private:

    const GridImp* identityGrid_;

    HostLevelIntersectionIterator hostIterator_;

  };


}  // namespace Dune

#endif
