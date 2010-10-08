// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_STRUCTURED_GRID_FACTORY_HH
#define DUNE_STRUCTURED_GRID_FACTORY_HH

/** \file
    \brief A class to construct structured cube and simplex grids using the grid factory
 */

#include <cstddef>

#include <dune/common/array.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/grid/common/gridfactory.hh>

namespace Dune {

  /** \brief Construct structured cube and simplex grids in unstructured grid managers
   */
  template <class GridType>
  class StructuredGridFactory
  {

    typedef typename GridType::ctype ctype;

    static const int dim = GridType::dimension;

    static const int dimworld = GridType::dimensionworld;

    /** \brief dim-dimensional multi-index.  The range for each component can be set individually
     */
    class MultiIndex
      : public array<unsigned int,dim>
    {

      // The range of each component
      array<unsigned int,dim> limits_;

    public:
      /** \brief Constructor with a given range for each digit */
      MultiIndex(const array<unsigned int,dim>& limits)
        : limits_(limits)
      {
        std::fill(this->begin(), this->end(), 0);
      }

      /** \brief Increment the MultiIndex */
      MultiIndex& operator++() {

        for (int i=0; i<dim; i++) {

          // Augment digit
          (*this)[i]++;

          // If there is no carry-over we can stop here
          if ((*this)[i]<limits_[i])
            break;

          (*this)[i] = 0;

        }
        return *this;
      }

      /** \brief Compute how many times you can call operator++ before getting to (0,...,0) again */
      size_t cycle() const {
        size_t result = 1;
        for (int i=0; i<dim; i++)
          result *= limits_[i];
        return result;
      }

    };

    /** \brief Insert a structured set of vertices into the factory */
    static void insertVertices(GridFactory<GridType>& factory,
                               const FieldVector<ctype,dim>& lowerLeft,
                               const FieldVector<ctype,dim>& upperRight,
                               const array<unsigned int,dim>& vertices)
    {

      MultiIndex index(vertices);

      // Compute the total number of vertices to be created
      int numVertices = index.cycle();

      // Create vertices
      for (int i=0; i<numVertices; i++, ++index) {

        // scale the multiindex to obtain a world position
        FieldVector<double,dimworld> pos(0);
        for (int j=0; j<dim; j++)
          pos[j] = index[j] * (upperRight[j]-lowerLeft[j])/(vertices[j]-1);

        factory.insertVertex(pos);

      }

    }

    // Compute the index offsets needed to move to the adjacent vertices
    // in the different coordinate directions
    static array<unsigned int, dim> computeUnitOffsets(const array<unsigned int,dim>& vertices)
    {
      array<unsigned int, dim> unitOffsets;
      if (dim>0)        // paranoia
        unitOffsets[0] = 1;

      for (int i=1; i<dim; i++)
        unitOffsets[i] = unitOffsets[i-1] * vertices[i-1];

      return unitOffsets;
    }

  public:

    /** \brief Create a structured cube grid
        \param lowerLeft Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements Number of elements in each coordinate direction
     */
    static shared_ptr<GridType> createCubeGrid(const FieldVector<ctype,dimworld>& lowerLeft,
                                               const FieldVector<ctype,dimworld>& upperRight,
                                               const array<unsigned int,dim>& elements)
    {
      // The grid factory
      GridFactory<GridType> factory;

      // Insert uniformly spaced vertices
      array<unsigned int,dim> vertices = elements;
      for (int i=0; i<vertices.size(); i++)
        vertices[i]++;

      // Insert vertices for structured grid into the factory
      insertVertices(factory, lowerLeft, upperRight, vertices);

      // Compute the index offsets needed to move to the adjacent vertices
      // in the different coordinate directions
      array<unsigned int, dim> unitOffsets = computeUnitOffsets(vertices);

      // Compute an element template (the cube at (0,...,0).  All other
      // cubes are constructed by moving this template around
      unsigned int nCorners = 1<<dim;

      std::vector<unsigned int> cornersTemplate(nCorners,0);

      for (size_t i=0; i<nCorners; i++)
        for (int j=0; j<dim; j++)
          if ( i & (1<<j) )
            cornersTemplate[i] += unitOffsets[j];

      // Insert elements
      MultiIndex index(elements);

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

        factory.insertElement(GeometryType(GeometryType::cube, dim), corners);

      }

      // Create the grid and hand it to the calling method
      return shared_ptr<GridType>(factory.createGrid());

    }

    /** \brief Create a structured simplex grid

        This works in all dimensions.  The Coxeter-Freudenthal-Kuhn triangulation is
        used, which splits each cube into dim! simplices.  See Allgower and Georg,
        'Numerical Path Following' for a description.
     */
    static shared_ptr<GridType> createSimplexGrid(const FieldVector<ctype,dimworld>& lowerLeft,
                                                  const FieldVector<ctype,dimworld>& upperRight,
                                                  const array<unsigned int,dim>& elements)
    {
      // The grid factory
      GridFactory<GridType> factory;

      // Insert uniformly spaced vertices
      array<unsigned int,dim> vertices = elements;
      for (std::size_t i=0; i<vertices.size(); i++)
        vertices[i]++;

      insertVertices(factory, lowerLeft, upperRight, vertices);

      // Compute the index offsets needed to move to the adjacent vertices
      // in the different coordinate directions
      array<unsigned int, dim> unitOffsets = computeUnitOffsets(vertices);

      // Insert the elements
      std::vector<unsigned int> corners(dim+1);

      // Loop over all "cubes", and split up each cube into dim! (factorial) simplices
      MultiIndex elementsIndex(elements);
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
            corners[j+1] = corners[j] + unitOffsets[permutation[j]];

          factory.insertElement(GeometryType(GeometryType::simplex, dim), corners);

        } while (std::next_permutation(permutation.begin(), permutation.end()));

      }

      // Create the grid and hand it to the calling method
      return shared_ptr<GridType>(factory.createGrid());
    }

  };

}  // namespace Dune

#endif
