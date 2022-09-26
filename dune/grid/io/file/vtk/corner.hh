// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_CORNER_HH
#define DUNE_GRID_IO_FILE_VTK_CORNER_HH

#include <dune/grid/io/file/vtk/common.hh>

namespace Dune
{
  //! \addtogroup VTK
  //! \{

  namespace VTK {

    //! simple class representing a corner of a cell
    /**
     * \tparam Cell Type of element this is a corner for.  This can
     *              essentially be anything with a type() method, in
     *              particular an Entity or an Intersection.
     */
    template<typename Cell>
    class Corner {
      // store a pointer to the element
      const Cell* cell_;
      // store index of the corner within element (Dune numbering)
      unsigned index;

    public:
      //! construct a Corner
      /**
       * \param cell      Reference to the cell
       * \param duneIndex Index of the corner within the element in
       *                  Dune-numbering
       */
      Corner(const Cell& cell, unsigned duneIndex)
        : cell_(&cell), index(duneIndex)
      { }

      //! construct an invalid Corner
      Corner() { }

      //! get reference to the cell
      const Cell& cell() const { return *cell_; }
      //! set a new cell
      /**
       * This also resets the index of the element to 0 (Dune-numbering), so
       * if you want to set both element and index, set the cell first.
       */
      void cell(const Cell& cell__) { cell_ = &cell__; index = 0; }

      //! get the index of the corner within the cell in Dune-numbering
      unsigned duneIndex() const { return index; }
      //! set the index of the corner within the cell in Dune-numbering
      void duneIndex(unsigned i) { index = i; }

      //! get the index of the corner within the cell in VTK-numbering
      /**
       * This requires that the cell is valid
       */
      unsigned vtkIndex() const { return renumber(cell_->type(), index); }
      //! set the index of the corner within the cell in VTK-numbering
      /**
       * This requires that the cell is valid
       */
      void vtkIndex(unsigned i) { index = renumber(cell_->type(), i); }
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_CORNER_HH
