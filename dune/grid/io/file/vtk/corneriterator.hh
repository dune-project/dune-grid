// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_CORNERITERATOR_HH
#define DUNE_GRID_IO_FILE_VTK_CORNERITERATOR_HH

#include <iterator>

#include <dune/common/iteratorfacades.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/io/file/vtk/corner.hh>

namespace Dune
{
  //! \addtogroup VTK
  //! \{

  namespace VTK {

    //! iterate over the corners of some cell range
    /**
     * This will visit all the corners of all elements visited by
     * CellIterator.
     */
    template<typename CellIterator>
    class CornerIterator
      : public ForwardIteratorFacade
        < CornerIterator<CellIterator>,
            const Corner<typename std::remove_const<typename std::iterator_traits<
                        CellIterator>::value_type>::type>,
            const Corner<typename std::remove_const<typename std::iterator_traits<
                        CellIterator>::value_type>::type>&,
            typename std::iterator_traits<CellIterator>::difference_type>
    {
    public:
      // reiterate the facades typedefs here
      typedef CornerIterator<CellIterator> DerivedType;
      typedef VTK::Corner<typename std::remove_const<typename std::iterator_traits<
                  CellIterator>::value_type>::type> Corner;
      typedef const Corner Value;
      typedef Value& Reference;
      typedef typename std::iterator_traits<CellIterator>::difference_type
      DifferenceType;

      typedef typename std::iterator_traits<CellIterator>::value_type::Geometry::ctype
      ctype;
      static const unsigned dim = std::iterator_traits<CellIterator>::
                                  value_type::mydimension;
      typedef ReferenceElements<ctype, dim> Refelems;

    private:
      typedef ForwardIteratorFacade<DerivedType, Value, Reference,
          DifferenceType> Facade;

      CellIterator cellit;
      CellIterator cellend;
      Corner corner;

    public:
      Reference dereference() const {
        return corner;
      }

      bool isDereferencable() const {
        return cellit != cellend;
      }

      bool equals(const DerivedType& other) const {
        bool mePassedTheEnd = !isDereferencable();
        bool otherPassedTheEnd = !other.isDereferencable();
        // both are passed the end => return true
        if(mePassedTheEnd && otherPassedTheEnd) return true;
        // one is passed the end => return false
        if(mePassedTheEnd || otherPassedTheEnd) return false;
        // none is passed the end, do their iterators and indices match?
        return cellit == other.cellit &&
               corner.duneIndex() == other.corner.duneIndex();
      }

      void increment() {
        int index = corner.vtkIndex();
        ++index;
        if(index == Refelems::general(cellit->type()).size(dim)) {
          ++cellit;
          if(cellit != cellend) {
            corner.cell(*cellit);
            corner.vtkIndex(0);
          }
        }
        else
          corner.vtkIndex(index);
      }

      //! construct a CornerIterator
      /**
       * \param cellit_  The begin iterator of the undelying range.
       * \param cellend_ The end iterator of the underlying range.
       * \param vtkIndex VTKIndex of the currently pointed to corner.
       */
      CornerIterator(const CellIterator& cellit_, const CellIterator& cellend_,
                     unsigned vtkIndex = 0)
        : cellit(cellit_), cellend(cellend_)
      {
        if(cellit != cellend) {
          corner.cell(*cellit);
          corner.vtkIndex(vtkIndex);
        }
      }
      //! construct a CornerIterator
      /**
       * This constructs a passed-the-end iterator value.
       */
      CornerIterator(const CellIterator& cellend_)
        : cellit(cellend_), cellend(cellend_)
      { }
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_CORNERITERATOR_HH
