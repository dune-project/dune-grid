// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_INTERSECTIONS_HH
#define DUNE_UGGRID_INTERSECTIONS_HH

#include <dune/common/sllist.hh>
#include <dune/common/shared_ptr.hh>

/** \file
 * \brief The UGGridLeafIntersection and UGGridLevelIntersection classes
 */

namespace Dune {

  /** \brief Implementation class for an intersection with an element on the same level */
  template<class GridImp>
  class UGGridLevelIntersection
  {
  public:
    enum {dim=GridImp::dimension};
    enum {dimworld=GridImp::dimensionworld};

  private:
    friend class UGGridEntity<0,dim,GridImp>;

    // The type used to store coordinates
    typedef typename GridImp::ctype UGCtype;

    // The corresponding iterator needs to access all members
    friend class UGGridLevelIntersectionIterator<GridImp>;

    typedef FieldVector<UGCtype, dimworld> WorldVector;
    typedef FieldVector<UGCtype, dim-1> FaceVector;

  public:

    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;

    /** The default Constructor makes empty Iterator
        \todo Should be private
     */
    UGGridLevelIntersection(typename UG_NS<dim>::Element* center, int nb)
      : selfLocal_(UGGridGeometry<dim-1,dimworld,GridImp>()),
        neighLocal_(UGGridGeometry<dim-1,dimworld,GridImp>()),
        neighGlob_(UGGridGeometry<dim-1,dimworld,GridImp>()),
        center_(center), neighborCount_(nb)
    {}

    //! equality
    bool equals(const UGGridLevelIntersection<GridImp>& i) const {
      return center_==i.center_ && neighborCount_ == i.neighborCount_;
    }

    //! return EntityPointer to the Entity on the inside of this intersection
    //! (that is the entity where we started this iterator)
    EntityPointer inside() const {
      return UGGridEntityPointer<0,GridImp>(center_);
    }

    //! return EntityPointer to the Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    EntityPointer outside() const {
      typename UG_NS<dim>::Element* otherelem = UG_NS<dim>::NbElem(center_, neighborCount_);

      if (otherelem==0)
        DUNE_THROW(GridError,"no neighbor found in outside()");

      return UGGridEntityPointer<0,GridImp>(otherelem);
    }

    //! return true if intersection is with boundary.
    bool boundary () const {
      return UG_NS<dim>::Side_On_Bnd(center_, neighborCount_);
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return UG_NS<dim>::NbElem(center_, neighborCount_) != NULL;
    }

    //! return information about the Boundary
    int boundaryId () const DUNE_DEPRECATED {
      return boundarySegmentIndex();
    }

    /** \brief return index of the corresponding coarse grid boundary segment */
    int boundarySegmentIndex () const {
#ifndef NDEBUG
      if (!boundary())
        DUNE_THROW(GridError, "Calling boundarySegmentIndex() for a non-boundary intersection!");
#endif
      return UG_NS<dim>::boundarySegmentIndex(center_, neighborCount_);
    }

    /** \brief Returns true, because UG level intersections are always conforming */
    bool conforming() const {
      return true;
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    const LocalGeometry &geometryInInside () const;

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    const Geometry &geometry () const;

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return geometryInInside().type();
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    const LocalGeometry &geometryInOutside () const;

    //! local index of codim 1 entity in self where intersection is contained in
    int indexInInside () const
    {
      const int number = UGGridRenumberer<dim>::facesUGtoDUNE(neighborCount_, UG_NS<dimworld>::Sides_Of_Elem(center_));
      typedef GenericGeometry::MapNumberingProvider< dim > Numbering;
      return Numbering::template dune2generic< 1 >( inside()->type().id(), number );
    }

    //! local index of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const;

    //! return outer normal
    const WorldVector&
    outerNormal (const FaceVector& local) const;

    //! return outer normal
    const FieldVector<UGCtype, dimworld>&
    integrationOuterNormal (const FieldVector<UGCtype, dim-1>& local) const
    {
      integrationOuterNormal_ = outerNormal(local);

      const UGCtype scale = geometry().integrationElement( local ) / integrationOuterNormal_.two_norm();
      integrationOuterNormal_ *= scale;

      return integrationOuterNormal_;
    }

    //! return outer normal
    const FieldVector<UGCtype, GridImp::dimensionworld>&
    unitOuterNormal (const FieldVector<UGCtype, dim-1>& local) const
    {
      unitOuterNormal_ = outerNormal(local);
      unitOuterNormal_ /= unitOuterNormal_.two_norm();
      return unitOuterNormal_;
    }

    //! return outer normal
    const FieldVector<UGCtype, GridImp::dimensionworld>&
    centerUnitOuterNormal () const
    {
      GeometryType type = geometry().type();
      const GenericReferenceElement<UGCtype, dim-1> & refElement =
        GenericReferenceElements<UGCtype, dim-1>::general(type);
      return unitOuterNormal(refElement.position(0,0));
    }

  private:

    //! vector storing the outer normal
    mutable FieldVector<UGCtype, dimworld> outerNormal_;
    mutable FieldVector<UGCtype, dimworld> integrationOuterNormal_;
    mutable FieldVector<UGCtype, dimworld> unitOuterNormal_;

    //! pointer to element holding the self_local and self_global information.
    //! This element is created on demand.
    mutable MakeableInterfaceObject<LocalGeometry> selfLocal_;
    mutable MakeableInterfaceObject<LocalGeometry> neighLocal_;

    //! pointer to element holding the neighbor_global and neighbor_local
    //! information.
    mutable MakeableInterfaceObject<Geometry> neighGlob_;

    //! The UG element the iterator was created from
    typename UG_NS<dim>::Element *center_;

    //! count on which neighbor we are looking at. Note that this is interpreted in UG's ordering!
    int neighborCount_;

  };


  /** \brief Implementation class for a leaf intersection in a UGGrid */
  template<class GridImp>
  class UGGridLeafIntersection
  {

    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    friend class UGGridEntity<0,dim,GridImp>;

    // The type used to store coordinates
    typedef typename GridImp::ctype UGCtype;

    // An element face identfied by the element and a face number
    typedef std::pair<const typename UG_NS<dim>::Element*, int> Face;

    // The corresponding iterator needs to access all members
    friend class UGGridLeafIntersectionIterator<GridImp>;

    typedef FieldVector<UGCtype, dimworld> WorldVector;
    typedef FieldVector<UGCtype, dim-1> FaceVector;

  public:

    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;

    UGGridLeafIntersection(typename UG_NS<dim>::Element* center, int nb)
      : selfLocal_(UGGridGeometry<dim-1,dimworld,GridImp>()),
        neighLocal_(UGGridGeometry<dim-1,dimworld,GridImp>()),
        neighGlob_(UGGridGeometry<dim-1,dimworld,GridImp>()),
        center_(center), neighborCount_(nb), subNeighborCount_(0)
    {
      if (neighborCount_ < UG_NS<dim>::Sides_Of_Elem(center_))
        constructLeafSubfaces();
    }

    //! equality
    bool equals(const UGGridLeafIntersection<GridImp>& other) const {
      return center_           == other.center_
             && neighborCount_    == other.neighborCount_
             && subNeighborCount_ == other.subNeighborCount_;
    }

    //! return EntityPointer to the Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    EntityPointer inside() const {
      return UGGridEntityPointer<0,GridImp>(center_);
    }

    //! return EntityPointer to the Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    EntityPointer outside() const {

      const typename UG_NS<dim>::Element* otherelem = leafSubFaces_[subNeighborCount_].first;

      if (otherelem==0)
        DUNE_THROW(GridError,"no neighbor found in outside()");

      /** \todo Remove the const_cast */
      return UGGridEntityPointer<0,GridImp>(const_cast<typename UG_NS<dim>::Element*>(otherelem));
    }

    //! return true if intersection is with boundary.
    bool boundary () const {
      return UG_NS<dim>::Side_On_Bnd(center_, neighborCount_);
    }

    //! return true if a neighbor element exists across this intersection
    bool neighbor () const {
      return leafSubFaces_[subNeighborCount_].first != NULL;
    }

    //! return information about the Boundary
    int boundaryId () const DUNE_DEPRECATED {
      return boundarySegmentIndex();
    }

    /** \brief Return index of corresponding coarse grid boundary segment */
    int boundarySegmentIndex () const {
#ifndef NDEBUG
      if (!boundary())
        DUNE_THROW(GridError, "Calling boundarySegmentIndex() for a non-boundary intersection!");
#endif
      return UG_NS<dim>::boundarySegmentIndex(center_, neighborCount_);
    }

    /** \brief Is this intersection conforming? */
    bool conforming() const {

      const typename UG_NS<dim>::Element* outside = leafSubFaces_[subNeighborCount_].first;

      if (outside == NULL         // boundary intersection
          // inside and outside are on the same level
          || UG_NS<dim>::myLevel(outside) == UG_NS<dim>::myLevel(center_)
          // outside is on a higher level, but there is only one intersection
          || (UG_NS<dim>::myLevel(outside) > UG_NS<dim>::myLevel(center_)
              && leafSubFaces_.size()==1))
        return true;

      // outside is on a lower level.  we have to check whether vertices match
      int numInsideIntersectionVertices  = UG_NS<dim>::Corners_Of_Side(center_, neighborCount_);
      int numOutsideIntersectionVertices = UG_NS<dim>::Corners_Of_Side(outside, leafSubFaces_[subNeighborCount_].second);
      if (numInsideIntersectionVertices != numOutsideIntersectionVertices)
        return false;

      // Loop over all vertices of the face of this element that corresponds to this intersection
      for (int i=0; i<numInsideIntersectionVertices; i++) {

        const typename UG_NS<dim>::Vertex* insideVertex = UG_NS<dim>::Corner(center_, UG_NS<dim>::Corner_Of_Side(center_, neighborCount_, i))->myvertex;

        // Loop over all vertices of the corresponding element side of the outside element
        bool vertexFound = false;
        for (int j=0; j<numOutsideIntersectionVertices; j++) {

          // get vertex
          const typename UG_NS<dim>::Vertex* outsideVertex = UG_NS<dim>::Corner(outside, UG_NS<dim>::Corner_Of_Side(outside, leafSubFaces_[subNeighborCount_].second, j))->myvertex;

          // Stop if we have found corresponding vertices
          if (insideVertex==outsideVertex) {
            vertexFound = true;
            break;
          }

        }

        // One of this face's vertices has not been found in the face of the outside element
        if (vertexFound == false)
          return false;

      }

      return true;
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    const LocalGeometry &geometryInInside () const;

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    const Geometry &geometry () const;

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    const LocalGeometry &geometryInOutside () const;

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return geometryInInside().type();
    }

    //! local index of codim 1 entity in self where intersection is contained in
    int indexInInside () const
    {
      const int number = UGGridRenumberer<dim>::facesUGtoDUNE(neighborCount_, UG_NS<dimworld>::Sides_Of_Elem(center_));

      typedef GenericGeometry::MapNumberingProvider< dim > Numbering;
      return Numbering::template dune2generic< 1 >( inside()->type().id(), number );
    }

    //! local index of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const;

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    const WorldVector&
    outerNormal (const FaceVector& local) const;

    //! return outer normal
    const FieldVector<UGCtype, dimworld>&
    integrationOuterNormal (const FieldVector<UGCtype, dim-1>& local) const
    {
      integrationOuterNormal_ = outerNormal(local);

      //integrationOuterNormal_ /= integrationOuterNormal_.two_norm();
      //integrationOuterNormal_ *= geometry().integrationElement(local);

      const UGCtype scale = geometry().integrationElement( local ) / integrationOuterNormal_.two_norm();
      integrationOuterNormal_ *= scale;

      return integrationOuterNormal_;
    }

    //! return outer normal
    const FieldVector<UGCtype, dimworld>&
    unitOuterNormal (const FieldVector<UGCtype, dim-1>& local) const {
      unitOuterNormal_ = outerNormal(local);
      unitOuterNormal_ /= unitOuterNormal_.two_norm();
      return unitOuterNormal_;
    }

    //! return outer normal
    const FieldVector<UGCtype, dimworld>&
    centerUnitOuterNormal () const
    {
      GeometryType type = geometry().type();
      const GenericReferenceElement<UGCtype, dim-1> & refElement =
        GenericReferenceElements<UGCtype, dim-1>::general(type);
      return unitOuterNormal(refElement.position(0,0));
    }

  private:
    //**********************************************************
    //  private methods
    //**********************************************************

    int numberInNeighbor(const typename UG_NS<dim>::Element* me, const typename UG_NS<dim>::Element* other) const {
      const int nSides = UG_NS<dim>::Sides_Of_Elem(other);

      for (int i=0; i<nSides; i++)
        if (UG_NS<dim>::NbElem(other,i) == me)
          return i;

      // this point should not be reached, otherwise throw exception
      DUNE_THROW(InvalidStateException,"no consistency in numberInNeighbor");
      return -1;
    }

    /** \brief Find the topological father face of a given fact*/
    int getFatherSide(const Face& currentFace) const;

    /** \brief Precompute list of all leaf intersections of the current element face */
    void constructLeafSubfaces();

    //! vector storing the outer normal
    mutable FieldVector<UGCtype, dimworld> outerNormal_;
    mutable FieldVector<UGCtype, dimworld> integrationOuterNormal_;
    mutable FieldVector<UGCtype, dimworld> unitOuterNormal_;

    //! pointer to element holding the self_local and self_global information.
    //! This element is created on demand.
    mutable MakeableInterfaceObject<LocalGeometry> selfLocal_;
    mutable MakeableInterfaceObject<LocalGeometry> neighLocal_;

    //! pointer to element holding the neighbor_global and neighbor_local
    //! information.
    mutable MakeableInterfaceObject<Geometry> neighGlob_;

    //! The UG element the iterator was created from
    typename UG_NS<dim>::Element *center_;

    //! count on which neighbor we are lookin' at. Note that this is interpreted in UG's ordering!
    int neighborCount_;

    /** \brief List of precomputed intersections */
    std::vector<Face> leafSubFaces_;

    /** \brief Current position in the leafSubFaces_ array */
    unsigned int subNeighborCount_;

  };

}  // namespace Dune

#endif
