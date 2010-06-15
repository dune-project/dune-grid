// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_BOUNDARYITERATORS_HH
#define DUNE_GRID_IO_FILE_VTK_BOUNDARYITERATORS_HH

#include <iterator>

#include <dune/common/iteratorfacades.hh>
#include <dune/common/shared_ptr.hh>

namespace Dune {

  //! \addtogroup VTK
  //! \{

  /** @file
      @author Jö Fahlke
      @brief Functions for VTK output on the skeleton
   */

  namespace VTK {

    //! iterate over the GridViews boundary intersections
    /**
     * This will visit all intersections for which boundary() is true and
     * neighbor() is false.
     */
    template<typename GV>
    class BoundaryIterator
      : public ForwardIteratorFacade
        < BoundaryIterator<GV>,
            const typename GV::Intersection,
            const typename GV::Intersection&,
            typename std::iterator_traits<typename GV::template Codim<0>::
                Iterator>::difference_type>
    {
    public:
      // reiterator the facades typedefs here
      typedef BoundaryIterator<GV> DerivedType;
      typedef const typename GV::Intersection Value;
      typedef Value& Reference;
      typedef typename GV::template Codim<0>::Iterator ElementIterator;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename std::iterator_traits<ElementIterator>::difference_type
      DifferenceType;

    private:
      typedef ForwardIteratorFacade<DerivedType, Value, Reference,
          DifferenceType> Facade;

      const GV* gv;
      ElementIterator eit;
      shared_ptr<IntersectionIterator> iit;

      bool valid() const {
        // we're valid if we're passed-the-end
        if(eit == gv->template end<0>()) return true;
        // or if we're on a boundary
        if((*iit)->boundary() && !(*iit)->neighbor()) return true;
        // otherwise we're invalid
        return false;
      }

      void basic_increment() {
        ++*iit;
        if(*iit == gv->iend(*eit)) {
          iit.reset();
          ++eit;
          if(eit != gv->template end<0>())
            iit.reset(new IntersectionIterator(gv->ibegin(*eit)));
        }
      }

    public:
      Reference dereference() const {
        return **iit;
      }
      bool equals(const DerivedType& other) const {
        if(eit != other.eit) return false;

        // this is a bit tricky, since we may not compare iit if we are
        // passed-the-end
        bool mePassedTheEnd = eit == gv->template end<0>();
        bool otherPassedTheEnd = other.eit == other.gv->template end<0>();

        // both passed-the-end => consider them equal
        if(mePassedTheEnd && otherPassedTheEnd) return true;

        // one passed the end => not equal
        if(mePassedTheEnd || otherPassedTheEnd) return false;

        // none passed-the-end => do their iit iterators match?
        return *iit == *other.iit;
      }

      void increment() {
        basic_increment();
        while(!valid()) basic_increment();
      }

      //! construct a BoundaryIterator
      /**
       * The iterator will initially point to the intersection iit_.  If that
       * intersection is not valid, it will advance to the first valid one.
       */
      BoundaryIterator(const GV& gv_, const ElementIterator& eit_,
                       const IntersectionIterator& iit_)
        : gv(&gv_), eit(eit_), iit(new IntersectionIterator(iit_))
      {
        while(!valid()) basic_increment();
      }
      //! construct a BoundaryIterator
      /**
       * The iterator will initially point to eit_'s first intersection.  If
       * that intersection is not valid, it will advance to the first valid
       * one.
       */
      BoundaryIterator(const GV& gv_, const ElementIterator& eit_)
        : gv(&gv_), eit(eit_)
      {
        if(eit != gv->template end<0>())
          iit.reset(new IntersectionIterator(gv->ibegin(*eit)));

        while(!valid()) basic_increment();
      }
      //! construct a BoundaryIterator
      /**
       * If end == true, construct an end iterator for the given gridview.
       * Otherwise, construct a begin iterator.
       */
      BoundaryIterator(const GV& gv_, bool end = false)
        : gv(&gv_), eit(end ? gv->template end<0>() : gv->template begin<0>())
      {
        if(eit != gv->template end<0>())
          iit.reset(new IntersectionIterator(gv->ibegin(*eit)));

        while(!valid()) basic_increment();
      }
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_BOUNDARYITERATORS_HH
