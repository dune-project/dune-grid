// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ONEDGRID_FACTORY_HH
#define DUNE_ONEDGRID_FACTORY_HH

/** \file
    \brief The specialization of the generic GridFactory for OneDGrid
    \author Oliver Sander
 */

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/onedgrid.hh>

namespace Dune {

  /** \brief Specialization of the generic GridFactory for OneDGrid

   */
  template <>
  class GridFactory<OneDGrid> : public GridFactoryInterface<OneDGrid> {

    typedef GridFactoryInterface<OneDGrid> Base;

    /** \brief Type used by the grid for coordinates */
    typedef OneDGrid::ctype ctype;

  public:

    // use default implementation from base class
    using Base::insertionIndex;

    /** \brief Default constructor */
    GridFactory();

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

    /** \brief Insert a boundary segment (== a point) and the boundary segment geometry
     *
        This influences the ordering of the boundary segments.
        The BoundarySegment object does not actually have any effect.
     */
    virtual void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                                       const std::shared_ptr<BoundarySegment<1> >& boundarySegment);

    /** \brief Return true if the intersection has been explictily insterted into the factory */
    virtual bool wasInserted(const typename OneDGrid::LeafIntersection& intersection) const;

    /** \brief Return the number of the intersection in the order of insertion into the factory */
    virtual unsigned int insertionIndex(const typename OneDGrid::LeafIntersection& intersection) const;

    /** \brief Finalize grid creation and hand over the grid
        The receiver takes responsibility of the memory allocated for the grid
     */
    virtual std::unique_ptr<OneDGrid> createGrid();

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
    std::vector<std::array<unsigned int, 2> > elements_;

    /** \brief Buffer the vertices until createGrid() is called */
    std::map<FieldVector<ctype,1>, unsigned int > vertexPositions_;

    /** \brief Counter that creates the vertex indices */
    unsigned int vertexIndex_;

    /** \brief Store the explicitly given boundary segments until createGrid() is called */
    std::vector<unsigned int> boundarySegments_;

    /** \brief Store vertex positions sorted by index */
    std::vector<ctype> vertexPositionsByIndex_;
  };

}

#endif
