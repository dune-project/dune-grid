// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ONEDGRID_FACTORY_HH
#define DUNE_ONEDGRID_FACTORY_HH

/** \file
    \brief The specialization of the generic GridFactory for OneDGrid
    \author Oliver Sander
 */

#include <vector>
#include <map>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/onedgrid.hh>

namespace Dune {

  /** \brief Specialization of the generic GridFactory for OneDGrid

   */
  template <>
  class GridFactory<OneDGrid> : public GridFactoryInterface<OneDGrid> {

    /** \brief Type used by the grid for coordinates */
    typedef OneDGrid::ctype ctype;

    typedef std::map<FieldVector<ctype,1>, unsigned int >::iterator VertexIterator;


  public:

    /** \brief Default constructor, optionally with grid-specific options
     *
     * \note There are no grid-specific options for OneDGrid at the
     *       moment, so any parameters specified have no effect.
     */
    explicit GridFactory(const ParameterTree &params = ParameterTree());

    /** \brief Constructor for a given grid object

       If you already have your grid object constructed you can
       hand it over using this constructor.

       If you construct your factory class using this constructor
       the pointer handed over to you by the method createGrid() is
       the one you supplied here.
     */
    GridFactory(OneDGrid* grid);

    /** \brief Destructor */
    ~GridFactory();

    /** \brief Insert a vertex into the coarse grid */
    virtual void insertVertex(const FieldVector<ctype,1>& pos);

    /** \brief Insert an element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
     */
    virtual void insertElement(const GeometryType& type,
                               const std::vector<unsigned int>& vertices);


    /** \brief Insert a boundary segment (== a point).
        This influences the ordering of the boundary segments
     */
    virtual void insertBoundarySegment(const std::vector<unsigned int>& vertices);

    /** \brief Finalize grid creation and hand over the grid

       The receiver takes responsibility of the memory allocated for the grid
     */
    virtual OneDGrid* createGrid();

  private:

    // Initialize the grid structure in UG
    void createBegin();

    // Pointer to the grid being built
    OneDGrid* grid_;

    // True if the factory allocated the grid itself, false if the
    // grid was handed over from the outside
    bool factoryOwnsGrid_;

    /** \brief While inserting the elements this array records the vertices
        of the elements. */
    std::vector<Dune::array<unsigned int, 2> > elements_;

    /** \brief Buffer the vertices until createGrid() is called */
    std::map<FieldVector<ctype,1>, unsigned int > vertexPositions_;

    /** \brief Counter that creates the vertex indices */
    unsigned int vertexIndex_;

    /** \brief Store the explicitly given boundary segments until createGrid() is called */
    std::vector<unsigned int> boundarySegments_;

  };

}

#endif
