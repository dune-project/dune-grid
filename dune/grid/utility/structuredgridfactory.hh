// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_STRUCTURED_GRID_FACTORY_HH
#define DUNE_STRUCTURED_GRID_FACTORY_HH

/** \file
    \brief A class to construct structured cube and simplex grids using the grid factory
 */

#include <algorithm>
#include <cstddef>
#include <utility>

#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
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

    /** \brief parse grid creation options from a parameter tree

        \param params   ParameterTree specifying how to create the grid.
                        The options noted below are interpreted, any other
                        options are ignored.  Options usually have
                        defaults, but if an option is present in the
                        ParameterTree is must parse successfully, otherwise
                        there will be exceptions.
        \param bbox     Parsed from the option \c bbox.lower and \c
                        bbox.upper.
        \param elements Parsed from the option \c elements.

        Options in params:
        \li \c bbox.lower Lower left corner of the grid (default: all 0)
        \li \c bbox.upper Upper right corner of the grid (default: all 1)
        \li \c elements   Number of elements in each coordinate direction
                          (default: all 1)
     */
    static void parseOptions(const ParameterTree &params,
                             std::pair<FieldVector<ctype,dimworld>,
                                 FieldVector<ctype,dimworld> > &bbox,
                             array<unsigned, dim> &elements)
    {
      std::fill(bbox.first.begin(), bbox.first.end(), 0);
      std::fill(bbox.second.begin(), bbox.second.end(), 1);
      std::fill(elements.begin(), elements.end(), 1);

      bbox.first = params.get("bbox.lower", bbox.first);
      bbox.second = params.get("bbox.upper", bbox.second);
      elements = params.get("elements", elements);
    }

  public:

    /** \brief Create a structured cube grid

        \param params ParameterTree specifying how to create the grid.
                      The options noted below are interpreted, and the
                      ParameterTree is then passed down to the
                      GridFactory.  In particular, the interpreted options
                      are not removed from the tree before passing it
                      down.

        Options in params:
        \li \c bbox.lower Lower left corner of the grid (default: all 0)
        \li \c bbox.upper Upper right corner of the grid (default: all 1)
        \li \c elements   Number of elements in each coordinate direction
                          (default: all 1)
     */
    static shared_ptr<GridType> createCubeGrid(const ParameterTree &params)
    {
      std::pair<FieldVector<ctype,dimworld>,
          FieldVector<ctype,dimworld> > bbox;
      array<unsigned, dim> elements;
      parseOptions(params, bbox, elements);
      return createCubeGrid(bbox.first, bbox.second, elements, params);
    }

    /** \brief Create a structured cube grid
        \param lowerLeft Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements Number of elements in each coordinate direction
        \param params ParameterTree with grid specific options.  This is
                      passed directly to the underlying GridFactory.  In
                      particular, this functions does not interpret the \c
                      bbox.lower, \c bbox.upper, or \c elements options
                      from the ParameterTree, use createCubeGrid(const
                      ParameterTree&) for that.
     */
    static shared_ptr<GridType>
    createCubeGrid(const FieldVector<ctype,dimworld>& lowerLeft,
                   const FieldVector<ctype,dimworld>& upperRight,
                   const array<unsigned int,dim>& elements,
                   const ParameterTree &params = ParameterTree())
    {
      // The grid factory
      GridFactory<GridType> factory(params);

      // Insert uniformly spaced vertices
      array<unsigned int,dim> vertices = elements;
      for( size_t i = 0; i < vertices.size(); ++i )
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

        \param params ParameterTree specifying how to create the grid.
                      The options noted below are interpreted, and the
                      ParameterTree is then passed down to the
                      GridFactory.  In particular, the interpreted options
                      are not removed from the tree before passing it
                      down.

        Options in params:
        \li \c bbox.lower Lower left corner of the grid (default: all 0)
        \li \c bbox.upper Upper right corner of the grid (default: all 1)
        \li \c elements   Number of elements in each coordinate direction
                          (default: all 1)
     */
    static shared_ptr<GridType>
    createSimplexGrid(const ParameterTree &params) {
      std::pair<FieldVector<ctype,dimworld>,
          FieldVector<ctype,dimworld> > bbox;
      array<unsigned, dim> elements;
      parseOptions(params, bbox, elements);
      return createSimplexGrid(bbox.first, bbox.second, elements,
                               params);
    }

    /** \brief Create a structured simplex grid
        \param lowerLeft Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements Number of elements in each coordinate direction
        \param params ParameterTree with grid specific options.  This is
                      passed directly to the underlying GridFactory.  In
                      particular, this functions does not interpret the \c
                      bbox.lower, \c bbox.upper, or \c elements options
                      from the ParameterTree, use createCubeGrid(const
                      ParameterTree&) for that.

        This works in all dimensions.  The Coxeter-Freudenthal-Kuhn triangulation is
        used, which splits each cube into dim! simplices.  See Allgower and Georg,
        'Numerical Path Following' for a description.
     */
    static shared_ptr<GridType>
    createSimplexGrid(const FieldVector<ctype,dimworld>& lowerLeft,
                      const FieldVector<ctype,dimworld>& upperRight,
                      const array<unsigned int,dim>& elements,
                      const ParameterTree &params = ParameterTree())
    {
      // The grid factory
      GridFactory<GridType> factory(params);

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
