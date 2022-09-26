// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_POINTITERATOR_HH
#define DUNE_GRID_IO_FILE_VTK_POINTITERATOR_HH

#include <iterator>
#include <vector>

#include <dune/common/iteratorfacades.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/io/file/vtk/corner.hh>
#include <dune/grid/io/file/vtk/corneriterator.hh>

namespace Dune
{
  //! \addtogroup VTK
  //! \{

  namespace VTK {

    //! iterate over the points of some corner range
    /**
     * This iterator essentially filters the output from CornerIterator to
     * visit each point only once.
     *
     * \tparam CellIterator Type of the iterators over cells.  The usual
     *                      codim0 EntityIterators should work here.
     * \tparam IS           Type of the indexset.
     *
     * The IndexSet must provide methods and member types:
     * <ul>
     * <li>size(unsigned codimension) should return the number of vertices in
     *     the indexset.  The dimension of the cells will be passed as the
     *     codimension parameter value.
     * <li>subIndex(const Cell& cell, unsigned subIndex, unsigned codimension)
     *     should return the index of the point with local Dune index subIndex
     *     in cell cell.  The codimension parameter will again be the
     *     dimension of the cell.
     * <li>IndexType should be the type of the indices returned by subIndex().
     * </ul>
     * The requirements are crafted in such a way that the indexsets provided
     * by the grid should fulfill them as long as the cells are codim=0
     * Entities.
     */
    template<typename CellIterator, typename IS>
    class PointIterator
      : public ForwardIteratorFacade
        < PointIterator<CellIterator, IS>,
            const Corner<typename std::remove_const<typename std::iterator_traits<
                        CellIterator>::value_type>::type>,
            const Corner<typename std::remove_const<typename std::iterator_traits<
                        CellIterator>::value_type>::type>&,
            typename std::iterator_traits<CellIterator>::difference_type>
    {
    public:
      typedef VTK::Corner<typename std::remove_const<typename std::iterator_traits<
                  CellIterator>::value_type>::type> Corner;

      // reiterate the facades typedefs here
      typedef PointIterator<CellIterator, IS> DerivedType;
      typedef const Corner Value;
      typedef Value& Reference;
      typedef typename std::iterator_traits<CellIterator>::difference_type
      DifferenceType;

      static const unsigned mydim = std::iterator_traits<CellIterator>::
                                    value_type::mydimension;

    private:
      typedef ForwardIteratorFacade<DerivedType, Value, Reference,
          DifferenceType> Facade;

      CornerIterator<CellIterator> cornerit;
      const IS* is;
      std::vector<bool> seen;

    public:
      Reference dereference() const {
        return *cornerit;
      }

      bool isDereferencable() const {
        return cornerit.isDereferencable();
      }

      bool equals(const DerivedType& other) const {
        return cornerit == other.cornerit;
      }

      void increment() {
        for(++cornerit; isDereferencable(); ++cornerit) {
          typename IS::IndexType index =
            is->subIndex(cornerit->cell(), cornerit->duneIndex(), mydim);

          if(!seen[index]) {
            seen[index] = true;
            break;
          }
        }
      }

      //! construct a CornerIterator
      /**
       * \param cellit   The begin iterator of the undelying range.
       * \param cellend  The end iterator of the underlying range.
       * \param is_      A reference to the indexset to use.
       */
      PointIterator(const CellIterator& cellit, const CellIterator& cellend,
                    const IS& is_)
        : cornerit(cellit, cellend), is(&is_), seen(is->size(mydim), false)
      { }
      //! construct a CornerIterator
      /**
       * This constructs a passed-the-end iterator value.
       */
      PointIterator(const CellIterator& cellend_)
        : cornerit(cellend_), is(0)
      { }
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_POINTITERATOR_HH
