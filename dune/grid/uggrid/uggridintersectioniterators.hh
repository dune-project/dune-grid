// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGINTERSECTIONIT_HH
#define DUNE_UGINTERSECTIONIT_HH

/** \file
 * \brief The UGGridIntersectionIterator class
 */

namespace Dune {

  /** \brief Implementation of the UGGrid LevelIntersectionIterator
   */
  template<class GridImp>
  class UGGridLevelIntersectionIterator
  {

    enum {dim=GridImp::dimension};

    friend class UGGridEntity<0,dim,GridImp>;

  public:

    typedef Dune::Intersection<const GridImp, Dune::UGGridLevelIntersection> Intersection;

    /** The default Constructor makes empty Iterator
     */
    UGGridLevelIntersectionIterator(typename UG_NS<dim>::Element* center, int nb)
      : intersection_(UGGridLevelIntersection<GridImp>(center,nb))
    {}

    //! equality
    bool equals(const UGGridLevelIntersectionIterator<GridImp>& other) const {
      return GridImp::getRealImplementation(intersection_).equals(GridImp::getRealImplementation(other.intersection_));
    }

    //! prefix increment
    void increment() {
      GridImp::getRealImplementation(intersection_).neighborCount_++;

      GridImp::getRealImplementation(intersection_).geometryIsUpToDate_          = false;
      GridImp::getRealImplementation(intersection_).geometryInInsideIsUpToDate_  = false;
      GridImp::getRealImplementation(intersection_).geometryInOutsideIsUpToDate_ = false;
    }

    //! \brief dereferencing
    const Intersection & dereference() const {
      return intersection_;
    }

  private:

    mutable MakeableInterfaceObject<Intersection> intersection_;

  };


  /** \brief Implementation of the UGGrid LeafIntersectionIterator
   */
  template<class GridImp>
  class UGGridLeafIntersectionIterator
  {

    enum {dim=GridImp::dimension};

    friend class UGGridEntity<0,dim,GridImp>;

  public:

    typedef Dune::Intersection<const GridImp, Dune::UGGridLeafIntersection> Intersection;

    UGGridLeafIntersectionIterator(typename UG_NS<dim>::Element* center, int nb)
      : intersection_(UGGridLeafIntersection<GridImp>(center,nb))
    {}

    //! equality
    bool equals(const UGGridLeafIntersectionIterator<GridImp>& other) const {
      return GridImp::getRealImplementation(intersection_).equals(GridImp::getRealImplementation(other.intersection_));
    }

    /** \brief Prefix increment.

       The UG data structure does not directly contain information about leaf neighbors/intersections.
       Therefore getting that information is fairly expensive.  In particular, it is too expensive to
       start looking for the next intersection whenever 'increment()' is called.  Therefore, all
       intersections of one face of the 'inside' element are precomputed, and incrementing then traverses
       this precomputed list.  If the list is exhausted the iterator advances to the next faces
       and precomputes all intersections there.
     */
    void increment() {

      UGGridLeafIntersection<GridImp>& intersectionImp = GridImp::getRealImplementation(intersection_);

      intersectionImp.subNeighborCount_++;

      // are there no more intersections for the current element face?
      if (intersectionImp.subNeighborCount_ >= intersectionImp.leafSubFaces_.size() ) {

        // move to the next face
        intersectionImp.neighborCount_++;
        intersectionImp.subNeighborCount_ = 0;

        // if the next face is not the end iterator construct all intersections for it
        if (intersectionImp.neighborCount_ < UG_NS<dim>::Sides_Of_Elem(intersectionImp.center_))
          intersectionImp.constructLeafSubfaces();

      }

      // make sure geometries are not taken from the cache the next time
      // the geometry{|InInside|InOutside} methods are called.
      intersectionImp.geometryIsUpToDate_          = false;
      intersectionImp.geometryInInsideIsUpToDate_  = false;
      intersectionImp.geometryInOutsideIsUpToDate_ = false;
    }

    //! \brief dereferencing
    const Intersection & dereference() const {
      return intersection_;
    }

  private:

    mutable MakeableInterfaceObject<Intersection> intersection_;

  };

}  // namespace Dune

#endif
