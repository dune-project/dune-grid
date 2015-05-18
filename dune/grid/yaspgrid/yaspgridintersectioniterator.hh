// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDINTERSECTIONITERATOR_HH
#define DUNE_GRID_YASPGRIDINTERSECTIONITERATOR_HH

/** \file
 * \brief The YaspIntersectionIterator class

   YaspIntersectionIterator enables iteration over intersections with
   neighboring codim 0 entities.
 */

namespace Dune {

  /** \brief YaspIntersectionIterator enables iteration over intersections with
             neighboring codim 0 entities.
   */
  template<class GridImp>
  class YaspIntersectionIterator
  {
    enum { dim=GridImp::dimension };
  public:
    // types used from grids
    typedef Dune::YaspIntersection< GridImp > IntersectionImp;
    typedef Dune::Intersection< GridImp, IntersectionImp > Intersection;

    //! increment
    void increment()
    {
      GridImp::getRealImplementation(intersection_)._count += (GridImp::getRealImplementation(intersection_)._count < 2*dim);
    }

    //! equality
    bool equals (const YaspIntersectionIterator& other) const
    {
      return intersection_ == other.intersection_;
    }

    //! \brief dereferencing
    const Intersection & dereference() const
    {
      GridImp::getRealImplementation(intersection_).update();
      return intersection_;
    }

    YaspIntersectionIterator()
    {}

    //! make intersection iterator from entity
    YaspIntersectionIterator (const YaspEntity<0,dim,GridImp>& myself, bool toend)
      : intersection_(IntersectionImp(myself, toend))
    {}

    //! copy constructor
    YaspIntersectionIterator (const YaspIntersectionIterator& other)
      : intersection_(other.intersection_)
    {}

    //! assignment
    YaspIntersectionIterator & operator = (const YaspIntersectionIterator& other)
    {
      intersection_ = other.intersection_;
      return *this;
    }

  private:
    // The intersection this iterator points to
    mutable Intersection intersection_;
  };

}        // namespace Dune

#endif   // DUNE_GRID_YASPGRIDINTERSECTIONITERATOR_HH
