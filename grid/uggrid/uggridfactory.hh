// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_FACTORY_HH
#define DUNE_UGGRID_FACTORY_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/uggrid.hh>

namespace Dune {

  template <int dimworld>
  class GridFactory<UGGrid<dimworld> > : public GridFactoryInterface<UGGrid<dimworld> > {

    /** \brief Type used by the grid for coordinates */
    typedef typename UGGrid<dimworld>::ctype ctype;

    // UGGrid only in 2d and 3d
    dune_static_assert(dimworld==2 || dimworld || 3, "UGGrid only in 2d and 3d");

  public:

    /** \brief Default constructor */
    GridFactory();

    /** \brief Destructor */
    ~GridFactory();

    /** \brief Insert a vertex into the coarse grid */
    virtual void insertVertex(const FieldVector<ctype,dimworld>& pos);

    /** \brief Insert an element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
     */
    virtual void insertElement(const GeometryType& type,
                               const std::vector<unsigned int>& vertices);

    /** \brief Method to insert an arbitrarily shaped boundary segment into a coarse grid
        \param vertices The indices of the vertices of the segment
        \param boundarySegment Class implementing the geometry of the boundary segment.
        The grid object takes control of this object and deallocates it when destructing itself.
     */
    void insertBoundarySegment(const std::vector<unsigned int> vertices,
                               const BoundarySegment<dimworld>* boundarySegment);


    /** \brief Finalize grid creation and hand over the grid

       The receiver takes responsibility of the memory allocated for the grid
     */
    virtual UGGrid<dimworld>* createGrid();

  private:
    UGGrid<dimworld>* grid_;

  };

}

#endif
