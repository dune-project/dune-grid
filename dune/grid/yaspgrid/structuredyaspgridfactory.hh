// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_STRUCTUREDYASPGRIDFACTORY_HH
#define DUNE_GRID_YASPGRID_STRUCTUREDYASPGRIDFACTORY_HH

#include <memory>

#include <dune/grid/utility/structuredgridfactory.hh>

/** \file
 * \brief Specialization of the StructuredGridFactory class for YaspGrid
 */

namespace Dune
{
  /** \brief Specialization of the StructuredGridFactory for YaspGrid

      This allows a YaspGrid to be constructed using the
      StructuredGridFactory just like the unstructured Grids.  There are two
      limitations:
      \li YaspGrid does not support simplices
      \li If the lower left corner should not be at the origin, the second template parameter
          of Yaspgrid has to be chosen as Dune::EquidistantOffsetCoordinates<ctype,dim>.
   */
  template<class ctype, int dim>
  class StructuredGridFactory<YaspGrid<dim, EquidistantCoordinates<ctype,dim> > >
  {
    typedef YaspGrid<dim, EquidistantCoordinates<ctype,dim> > GridType;
    static const int dimworld = GridType::dimensionworld;

  public:
    /** \brief Create a structured cube grid

        \param lowerLeft  Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements   Number of elements in each coordinate direction

        \note The default variant of YaspGrid only supports lowerLeft at the origin.
              Use YaspGrid<dim, EquidistantOffsetCoordinates<ctype,dim> > instead
              for non-trivial origin.
     */
    static std::unique_ptr<GridType>
    createCubeGrid(const FieldVector<ctype,dimworld>& lowerLeft,
                   const FieldVector<ctype,dimworld>& upperRight,
                   const std::array<unsigned int,dim>& elements)
    {
      using std::abs;
      for(int d = 0; d < dimworld; ++d)
        if(abs(lowerLeft[d]) > abs(upperRight[d])*1e-10)
          DUNE_THROW(GridError, className<StructuredGridFactory>()
                     << "::createCubeGrid(): You have to use Yaspgrid<dim"
                     ", EquidistantOffsetCoordinates<ctype,dim> > as your"
                     "grid type for non-trivial origin." );

      // construct array of ints instead of unsigned ints
      std::array<int, dim> elem;
      std::copy(elements.begin(), elements.end(), elem.begin());

      return std::make_unique<GridType>(upperRight, elem,
                             std::bitset<dim>(), 1);  // default constructor of bitset sets to zero
    }

    /** \brief Create a structured simplex grid

        \note Simplices are not supported in YaspGrid, so this functions
              unconditionally throws a GridError.
     */
    static std::unique_ptr<GridType>
    createSimplexGrid(const FieldVector<ctype,dimworld>& lowerLeft,
                      const FieldVector<ctype,dimworld>& upperRight,
                      const std::array<unsigned int,dim>& elements)
    {
      DUNE_THROW(GridError, className<StructuredGridFactory>()
                 << "::createSimplexGrid(): Simplices are not supported "
                 "by YaspGrid.");
    }

  };

  /** \brief Specialization of the StructuredGridFactory for YaspGrid<EquidistantOffsetCoordinates>

      This allows a YaspGrid to be constructed using the
      StructuredGridFactory just like the unstructured Grids.  Only limitation:
      limitations:
      \li YaspGrid does not support simplices
   */
  template<class ctype, int dim>
  class StructuredGridFactory<YaspGrid<dim, EquidistantOffsetCoordinates<ctype,dim> > > {
    typedef YaspGrid<dim, EquidistantOffsetCoordinates<ctype,dim> > GridType;
    static const int dimworld = GridType::dimensionworld;

  public:
    /** \brief Create a structured cube grid

        \param lowerLeft  Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements   Number of elements in each coordinate direction
     */
    static std::unique_ptr<GridType>
    createCubeGrid(const FieldVector<ctype,dimworld>& lowerLeft,
                   const FieldVector<ctype,dimworld>& upperRight,
                   const std::array<unsigned int,dim>& elements)
    {
      // construct array of ints instead of unsigned ints
      std::array<int, dim> elem;
      std::copy(elements.begin(), elements.end(), elem.begin());

      return std::make_unique<GridType>(lowerLeft, upperRight, elem,
                             std::bitset<dim>(), 1);  // default constructor of bitset sets to zero
    }

    /** \brief Create a structured simplex grid

        \note Simplices are not supported in YaspGrid, so this functions
              unconditionally throws a GridError.
     */
    static std::unique_ptr<GridType>
    createSimplexGrid(const FieldVector<ctype,dimworld>& lowerLeft,
                      const FieldVector<ctype,dimworld>& upperRight,
                      const std::array<unsigned int,dim>& elements)
    {
      DUNE_THROW(GridError, className<StructuredGridFactory>()
                 << "::createSimplexGrid(): Simplices are not supported "
                 "by YaspGrid.");
    }

  };

}  // namespace Dune
#endif
