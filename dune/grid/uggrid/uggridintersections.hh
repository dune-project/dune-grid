// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_INTERSECTIONS_HH
#define DUNE_UGGRID_INTERSECTIONS_HH

#include <memory>

#include <dune/grid/uggrid/uggridrenumberer.hh>

/** \file
 * \brief The UGGridLeafIntersection and UGGridLevelIntersection classes
 */

namespace Dune {

  /** \brief Implementation class for an intersection with an element on the same level */
  template<class GridImp>
  class UGGridLevelIntersection
  {
  public:
    constexpr static int dim = GridImp::dimension;
    constexpr static int dimworld = GridImp::dimensionworld;

  private:
    friend class UGGridEntity<0,dim,GridImp>;

    // The type used to store coordinates
    typedef typename GridImp::ctype UGCtype;

    // The corresponding iterator needs to access all members
    friend class UGGridLevelIntersectionIterator<GridImp>;

    typedef FieldVector<UGCtype, dimworld> WorldVector;
    typedef FieldVector<UGCtype, dim-1> FaceVector;

    typedef typename GridImp::Traits::template Codim<1>::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim<1>::LocalGeometryImpl LocalGeometryImpl;

  public:
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;

    UGGridLevelIntersection()
      : center_(nullptr)
      , neighborCount_(-1) // fixed marker value for invalid intersections to make equals() work
      , gridImp_(nullptr)
    {}

    /** The default Constructor makes empty Iterator
        \todo Should be private
     */
    UGGridLevelIntersection(typename UG_NS<dim>::Element* center, int nb, const GridImp* gridImp)
      : center_(center), neighborCount_(nb),
        gridImp_(gridImp)
    {}

    //! equality
    bool equals(const UGGridLevelIntersection<GridImp>& i) const {
      return center_==i.center_ && neighborCount_ == i.neighborCount_;
    }

    //! return Entity on the inside of this intersection
    //! (that is the entity where we started this iterator)
    Entity inside() const {
      return Entity(UGGridEntity<0,dim,GridImp>(center_,gridImp_));
    }

    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    Entity outside() const {
      typename UG_NS<dim>::Element* otherelem = UG_NS<dim>::NbElem(center_, neighborCount_);

      if (otherelem==0)
        DUNE_THROW(GridError,"no neighbor found in outside()");

      return Entity(UGGridEntity<0,dim,GridImp>(otherelem,gridImp_));
    }

    //! return true if intersection is with boundary.
    bool boundary () const {
      return UG_NS<dim>::Side_On_Bnd(center_, neighborCount_);
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return UG_NS<dim>::NbElem(center_, neighborCount_) != nullptr;
    }

    /** \brief return index of the corresponding coarse grid boundary segment */
    size_t boundarySegmentIndex () const {
#ifndef NDEBUG
      if (!boundary())
        DUNE_THROW(GridError, "Calling boundarySegmentIndex() for a non-boundary intersection!");
#endif
      UG_NS<dim>::Set_Current_BVP(gridImp_->multigrid_->theBVP);
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
    LocalGeometry geometryInInside () const;

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const;

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return geometryInInside().type();
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const;

    //! local index of codim 1 entity in self where intersection is contained in
    int indexInInside () const
    {
      return UGGridRenumberer<dim>::facesUGtoDUNE(neighborCount_, UG_NS<dim>::Tag(center_));
    }

    //! local index of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const;

    //! return outer normal
    const WorldVector&
    outerNormal (const FaceVector& local) const;

    //! return outer normal
    const WorldVector&
    integrationOuterNormal (const FaceVector& local) const
    {
      integrationOuterNormal_ = outerNormal(local);

      const UGCtype scale = geometry().integrationElement( local ) / integrationOuterNormal_.two_norm();
      integrationOuterNormal_ *= scale;

      return integrationOuterNormal_;
    }

    //! return outer normal
    const WorldVector&
    unitOuterNormal (const FaceVector& local) const
    {
      unitOuterNormal_ = outerNormal(local);
      unitOuterNormal_ /= unitOuterNormal_.two_norm();
      return unitOuterNormal_;
    }

    //! return outer normal
    const WorldVector&
    centerUnitOuterNormal () const
    {
      auto refElement = referenceElement(geometry());
      return unitOuterNormal(refElement.position(0,0));
    }

  private:

    //! vector storing the outer normal
    mutable WorldVector outerNormal_;
    mutable WorldVector integrationOuterNormal_;
    mutable WorldVector unitOuterNormal_;

    //! pointers holding the global and local geometries
    mutable std::shared_ptr<GeometryImpl>      geometry_;
    mutable std::shared_ptr<LocalGeometryImpl> geometryInInside_;
    mutable std::shared_ptr<LocalGeometryImpl> geometryInOutside_;

    //! The UG element the iterator was created from
    typename UG_NS<dim>::Element *center_;

    //! count on which neighbor we are looking at. Note that this is interpreted in UG's ordering!
    int neighborCount_;

    /** \brief The grid we belong to.  We need it to call set_Current_BVP */
    const GridImp* gridImp_;

  };


  /** \brief Implementation class for a leaf intersection in a UGGrid */
  template<class GridImp>
  class UGGridLeafIntersection
  {

    constexpr static int dim = GridImp::dimension;

    constexpr static int dimworld = GridImp::dimensionworld;

    friend class UGGridEntity<0,dim,GridImp>;

    // The type used to store coordinates
    typedef typename GridImp::ctype UGCtype;

    // An element face identfied by the element and a face number
    typedef std::pair<const typename UG_NS<dim>::Element*, int> Face;

    // The corresponding iterator needs to access all members
    friend class UGGridLeafIntersectionIterator<GridImp>;

    typedef FieldVector<UGCtype, dimworld> WorldVector;
    typedef FieldVector<UGCtype, dim-1> FaceVector;

    typedef typename GridImp::Traits::template Codim<1>::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim<1>::LocalGeometryImpl LocalGeometryImpl;

  public:
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;

    UGGridLeafIntersection()
      : center_(nullptr)
      , neighborCount_(-1)                  // fixed marker value for invalid intersections to make equals() work
      , subNeighborCount_(~unsigned(0))     // fixed marker value for invalid intersections to make equals() work
      , gridImp_(nullptr)
    {}

    UGGridLeafIntersection(typename UG_NS<dim>::Element* center, int nb, const GridImp* gridImp)
      : center_(center), neighborCount_(nb), subNeighborCount_(0),
        gridImp_(gridImp)
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

    //! return Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    Entity inside() const {
      return Entity(UGGridEntity<0,dim,GridImp>(center_,gridImp_));
    }

    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    Entity outside() const {

      const typename UG_NS<dim>::Element* otherelem = leafSubFaces_[subNeighborCount_].first;

      if (otherelem==0)
        DUNE_THROW(GridError,"no neighbor found in outside()");

      /** \todo Remove the const_cast */
      return Entity(UGGridEntity<0,dim,GridImp>(const_cast<typename UG_NS<dim>::Element*>(otherelem),gridImp_));
    }

    //! return true if intersection is with boundary.
    bool boundary () const {
      return UG_NS<dim>::Side_On_Bnd(center_, neighborCount_);
    }

    //! return true if a neighbor element exists across this intersection
    bool neighbor () const {
      return leafSubFaces_[subNeighborCount_].first != nullptr;
    }

    /** \brief Return index of corresponding coarse grid boundary segment */
    size_t boundarySegmentIndex () const {
#ifndef NDEBUG
      if (!boundary())
        DUNE_THROW(GridError, "Calling boundarySegmentIndex() for a non-boundary intersection!");
#endif
      UG_NS<dim>::Set_Current_BVP(gridImp_->multigrid_->theBVP);
      return UG_NS<dim>::boundarySegmentIndex(center_, neighborCount_);
    }

    /** \brief Is this intersection conforming? */
    bool conforming() const {

      const typename UG_NS<dim>::Element* outside = leafSubFaces_[subNeighborCount_].first;

      if (outside == nullptr         // boundary intersection
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
    LocalGeometry geometryInInside () const;

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const;

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const;

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return geometryInInside().type();
    }

    //! local index of codim 1 entity in self where intersection is contained in
    int indexInInside () const
    {
      return UGGridRenumberer<dim>::facesUGtoDUNE(neighborCount_, UG_NS<dimworld>::Tag(center_));
    }

    //! local index of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const;

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    const WorldVector&
    outerNormal (const FaceVector& local) const;

    //! return outer normal
    const WorldVector&
    integrationOuterNormal (const FaceVector& local) const
    {
      integrationOuterNormal_ = outerNormal(local);

      //integrationOuterNormal_ /= integrationOuterNormal_.two_norm();
      //integrationOuterNormal_ *= geometry().integrationElement(local);

      const UGCtype scale = geometry().integrationElement( local ) / integrationOuterNormal_.two_norm();
      integrationOuterNormal_ *= scale;

      return integrationOuterNormal_;
    }

    //! return outer normal
    const WorldVector&
    unitOuterNormal (const FaceVector& local) const {
      unitOuterNormal_ = outerNormal(local);
      unitOuterNormal_ /= unitOuterNormal_.two_norm();
      return unitOuterNormal_;
    }

    //! return outer normal
    const WorldVector&
    centerUnitOuterNormal () const
    {
      auto refElement = referenceElement(geometry());
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
    mutable WorldVector outerNormal_;
    mutable WorldVector integrationOuterNormal_;
    mutable WorldVector unitOuterNormal_;

    //! pointer to global and local intersection geometries
    mutable std::shared_ptr<GeometryImpl>      geometry_;
    mutable std::shared_ptr<LocalGeometryImpl> geometryInInside_;
    mutable std::shared_ptr<LocalGeometryImpl> geometryInOutside_;

    //! The UG element the iterator was created from
    typename UG_NS<dim>::Element *center_;

    //! count on which neighbor we are lookin' at. Note that this is interpreted in UG's ordering!
    int neighborCount_;

    /** \brief List of precomputed intersections */
    std::vector<Face> leafSubFaces_;

    /** \brief Current position in the leafSubFaces_ array */
    unsigned int subNeighborCount_;

    /** \brief The grid we belong to.  We need it to call set_Current_BVP */
    const GridImp* gridImp_;

  };

}  // namespace Dune

#endif
