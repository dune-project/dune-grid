// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_STRUCTURED_GRID_FACTORY_HH
#define DUNE_STRUCTURED_GRID_FACTORY_HH

/** \file
    \brief A class to construct structured cube and simplex grids using the grid factory
 */

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <memory>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/multiindex.hh>

namespace Dune {

  /** \brief Construct structured cube and simplex grids in unstructured grid managers
   */
  template <class GridType>
  class StructuredGridFactory
  {
    typedef typename GridType::ctype ctype;

    static const int dim = GridType::dimension;

    static const int dimworld = GridType::dimensionworld;

    /** \brief Insert a structured set of vertices into the factory */
    static void insertVertices(GridFactory<GridType>& factory,
                               const FieldVector<ctype,dimworld>& lowerLeft,
                               const FieldVector<ctype,dimworld>& upperRight,
                               const std::array<unsigned int,dim>& vertices)
    {
      FactoryUtilities::MultiIndex<dim> index(vertices);

      // Compute the total number of vertices to be created
      int numVertices = index.cycle();

      // Create vertices
      for (int i=0; i<numVertices; i++, ++index) {

        // scale the multiindex to obtain a world position
        FieldVector<double,dimworld> pos(0);
        for (int j=0; j<dim; j++)
          pos[j] = lowerLeft[j] + index[j] * (upperRight[j]-lowerLeft[j])/(vertices[j]-1);
        for (int j=dim; j<dimworld; j++)
          pos[j] = lowerLeft[j];

        factory.insertVertex(pos);

      }

    }

    // Compute the index offsets needed to move to the adjacent vertices
    // in the different coordinate directions
    static std::array<unsigned int, dim> computeUnitOffsets(const std::array<unsigned int,dim>& vertices)
    {
      std::array<unsigned int, dim> unitOffsets;
      if (dim>0)        // paranoia
        unitOffsets[0] = 1;

      for (int i=1; i<dim; i++)
        unitOffsets[i] = unitOffsets[i-1] * vertices[i-1];

      return unitOffsets;
    }

  public:

    /** \brief insert structured cube grid into grid factory

        If the grid dimension is less than the world dimension, the coefficients (dim+1,...,dimworld) in
        the vertex coordinates are set to the corresponding values of the lowerLeft input argument.

        \param factory grid factory used for creating the grid
        \param lowerLeft Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements Number of elements in each coordinate direction
     */
    static void createCubeGrid(
      GridFactory<GridType>& factory,
      const FieldVector<ctype,dimworld>& lowerLeft,
      const FieldVector<ctype,dimworld>& upperRight,
      const std::array<unsigned int,dim>& elements)
    {
      if (factory.comm().rank() == 0)
      {
        // Insert uniformly spaced vertices
        std::array<unsigned int,dim> vertices = elements;
        for( size_t i = 0; i < vertices.size(); ++i )
          vertices[i]++;

        // Insert vertices for structured grid into the factory
        insertVertices(factory, lowerLeft, upperRight, vertices);

        // Compute the index offsets needed to move to the adjacent
        // vertices in the different coordinate directions
        std::array<unsigned int, dim> unitOffsets =
          computeUnitOffsets(vertices);

        // Compute an element template (the cube at (0,...,0).  All
        // other cubes are constructed by moving this template around
        unsigned int nCorners = 1<<dim;

        std::vector<unsigned int> cornersTemplate(nCorners,0);

        for (size_t i=0; i<nCorners; i++)
          for (int j=0; j<dim; j++)
            if ( i & (1<<j) )
              cornersTemplate[i] += unitOffsets[j];

        // Insert elements
        FactoryUtilities::MultiIndex<dim> index(elements);

        // Compute the total number of elementss to be created
        int numElements = index.cycle();

        for (int i=0; i<numElements; i++, ++index) {

          // 'base' is the index of the lower left element corner
          unsigned int base = 0;
          for (int j=0; j<dim; j++)
            base += index[j] * unitOffsets[j];

          // insert new element
          std::vector<unsigned int> corners = cornersTemplate;
          for (size_t j=0; j<corners.size(); j++)
            corners[j] += base;

          factory.insertElement(GeometryTypes::cube(dim), corners);

        }

      }       // if(rank == 0)
    }

    /** \brief Create a structured cube grid

        If the grid dimension is less than the world dimension, the coefficients (dim+1,...,dimworld) in
        the vertex coordinates are set to the corresponding values of the lowerLeft input argument.

        \param lowerLeft Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements Number of elements in each coordinate direction
     */
    static std::unique_ptr<GridType> createCubeGrid(
      const FieldVector<ctype,dimworld>& lowerLeft,
      const FieldVector<ctype,dimworld>& upperRight,
      const std::array<unsigned int,dim>& elements)
    {
      GridFactory<GridType> factory;
      createCubeGrid(factory, lowerLeft, upperRight, elements);
      return std::unique_ptr<GridType>(factory.createGrid());
    }

    /** \brief insert structured simplex grid into grid factory

        This works in all dimensions.  The Coxeter-Freudenthal-Kuhn triangulation is
        used, which splits each cube into dim! (i.e., dim faculty) simplices.  See Allgower and Georg,
        'Numerical Path Following' for a description.

        If the grid dimension is less than the world dimension, the coefficients (dim+1,...,dimworld) in
        the vertex coordinates are set to the corresponding values of the lowerLeft input argument.

        \param factory grid factory used for creating the grid
        \param lowerLeft Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements Number of elements in each coordinate direction
     */
    static void createSimplexGrid(
      GridFactory<GridType>& factory,
      const FieldVector<ctype,dimworld>& lowerLeft,
      const FieldVector<ctype,dimworld>& upperRight,
      const std::array<unsigned int,dim>& elements)
    {
      if(factory.comm().rank() == 0)
      {
        // Insert uniformly spaced vertices
        std::array<unsigned int,dim> vertices = elements;
        for (std::size_t i=0; i<vertices.size(); i++)
          vertices[i]++;

        insertVertices(factory, lowerLeft, upperRight, vertices);

        // Compute the index offsets needed to move to the adjacent
        // vertices in the different coordinate directions
        std::array<unsigned int, dim> unitOffsets =
          computeUnitOffsets(vertices);

        // Loop over all "cubes", and split up each cube into dim!
        // (factorial) simplices
        FactoryUtilities::MultiIndex<dim> elementsIndex(elements);
        size_t cycle = elementsIndex.cycle();

        for (size_t i=0; i<cycle; ++elementsIndex, i++) {

          // 'base' is the index of the lower left element corner
          unsigned int base = 0;
          for (int j=0; j<dim; j++)
            base += elementsIndex[j] * unitOffsets[j];

          // each permutation of the unit vectors gives a simplex.
          std::vector<unsigned int> permutation(dim);
          for (int j=0; j<dim; j++)
            permutation[j] = j;

          do {

            // Make a simplex
            std::vector<unsigned int> corners(dim+1);
            corners[0] = base;

            for (int j=0; j<dim; j++)
              corners[j+1] =
                corners[j] + unitOffsets[permutation[j]];

            factory.insertElement(GeometryTypes::simplex(dim), corners);

          } while (std::next_permutation(permutation.begin(),
                                         permutation.end()));

        }

      }       // if(rank == 0)
    }

    /** \brief Create a structured simplex grid

        This works in all dimensions.  The Coxeter-Freudenthal-Kuhn triangulation is
        used, which splits each cube into dim! (i.e., dim faculty) simplices.  See Allgower and Georg,
        'Numerical Path Following' for a description.

        If the grid dimension is less than the world dimension, the coefficients (dim+1,...,dimworld) in
        the vertex coordinates are set to the corresponding values of the lowerLeft input argument.

        \param lowerLeft Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements Number of elements in each coordinate direction
     */
    static std::unique_ptr<GridType> createSimplexGrid(
      const FieldVector<ctype,dimworld>& lowerLeft,
      const FieldVector<ctype,dimworld>& upperRight,
      const std::array<unsigned int,dim>& elements)
    {
      GridFactory<GridType> factory;
      createSimplexGrid(factory, lowerLeft, upperRight, elements);
      return std::unique_ptr<GridType>(factory.createGrid());
    }

  };

}  // namespace Dune

#endif
