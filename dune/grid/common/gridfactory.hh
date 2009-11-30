// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_FACTORY_HH
#define DUNE_GRID_FACTORY_HH

/** \file
    \brief Provide a generic factory class for unstructured grids.
 */

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/boundarysegment.hh>

namespace Dune
{

  /** \brief Provide a generic factory class for unstructured grids.

      \ingroup GridFactory

      This base class declares the interface.
   */
  template <class GridType>
  class GridFactoryInterface
  {

  protected:
    /** \brief dimension of the grid */
    static const int dimension = GridType::dimension;

    /** \brief The grid world dimension */
    enum {dimworld = GridType::dimensionworld};

    /** \brief Type used by the grid for coordinates */
    typedef typename GridType::ctype ctype;

  public:
    template< int codim >
    struct Codim
    {
      typedef typename GridType::template Codim< codim >::Entity Entity;
    };

    /** \brief Default constructor */
    GridFactoryInterface()
    {}

    /** \brief Constructor for a given grid object

       If you already have your grid object constructed you can
       hand it over using this constructor.  A reason may be that
       you have need a UGGrid object with a non-default heap size.

       If you construct your factory class using this constructor
       the pointer handed over to you by the method createGrid() is
       the one you supplied here.
     */
    GridFactoryInterface(GridType* grid)
    {}

    virtual ~GridFactoryInterface ()
    {}

    /** \brief Insert a vertex into the coarse grid */
    virtual void insertVertex(const FieldVector<ctype,dimworld>& pos) = 0;

    /** \brief Insert an element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
     */
    virtual void insertElement(const GeometryType& type,
                               const std::vector<unsigned int>& vertices) = 0;

    /** \brief Method to insert an arbitrarily shaped boundary segment into a coarse grid
        \param vertices The indices of the vertices of the segment
        \param boundarySegment Class implementing the geometry of the boundary segment.
        The grid object takes control of this object and deallocates it when destructing itself.
     */
    virtual void insertBoundarySegment(const std::vector<unsigned int> vertices,
                                       const BoundarySegment<dimension,dimworld>* boundarySegment)
    {
      DUNE_THROW(GridError, "This grid does not support parametrized boundary segments!");
    }

    /** \brief Finalize grid creation and hand over the grid

       The receiver takes responsibility of the memory allocated for the grid
     */
    virtual GridType* createGrid() = 0;

    /** \brief obtain an element's insertion index
     *
     *  Data can be associated to the created macro grid using the insertion
     *  index of each entity that has been inserted during the grid creation
     *  process.
     *
     *  Between grid construction (createGrid) and the first grid
     *  modification, this method allows to obtain this insertion index from
     *  the grid factory. This way, data can be stored using the index maps
     *  provided by the grid.
     *
     *  \param[in]  entity  entity whose insertion index is requested
     *
     *  \returns insertion index of the entity
     */
    virtual unsigned int
    insertionIndex ( const typename Codim< 0 >::Entity &entity ) const
    {
      DUNE_THROW( NotImplemented, "insertion indices have not yet been implemented." );
    }

    /** \brief obtain a vertex' insertion index
     *
     *  Data can be associated to the created macro grid using the insertion
     *  index of each entity that has been inserted during the grid creation
     *  process.
     *
     *  Between grid construction (createGrid) and the first grid
     *  modification, this method allows to obtain this insertion index from
     *  the grid factory. This way, data can be stored using the index maps
     *  provided by the grid.
     *
     *  \param[in]  entity  entity whose insertion index is requested
     *
     *  \returns insertion index of the entity
     */
    virtual unsigned int
    insertionIndex ( const typename Codim< dimension >::Entity &entity ) const
    {
      DUNE_THROW( NotImplemented, "insertion indices have not yet been implemented." );
    }

    /** \brief obtain a boundary's insertion index
     *
     *  Data can be associated to the created macro grid using the insertion
     *  index of each entity that has been inserted during the grid creation
     *  process.
     *
     *  Between grid construction (createGrid) and the first grid
     *  modification, this method allows to obtain this insertion index from
     *  the grid factory. This way, data can be stored using the index maps
     *  provided by the grid.
     *
     *  \param[in]  intersection  intersection whose insertion index is requested
     *
     *  \returns insertion index of the intersection
     *
     *  \note The insertion index can only be obtained for boundary
     *        intersections that were actually inserted
     *        (see also wasInserted).
     */
    virtual unsigned int
    insertionIndex ( const typename GridType::LeafIntersection &intersection ) const
    {
      DUNE_THROW( NotImplemented, "insertion indices have not yet been implemented." );
    }


    /** \brief determine whether an intersection was inserted
     *
     *  This method allows checking wheter an intersection was actually
     *  inserted into the grid factory.
     *
     *  \note Not all boundary segments need to be inserted into the grid
     *        factory.
     *  \note This method returns \b false for all interior intersections
     *
     *  \param[in]  intersection  intersection in question
     *
     *  \returns \b true, if the intersection was inserted
     */
    virtual bool
    wasInserted ( const typename GridType::LeafIntersection &intersection ) const
    {
      DUNE_THROW( NotImplemented, "insertion indices have not yet been implemented." );
    }

  };


  /** \brief Provide a generic factory class for unstructured grids.

      \ingroup GridFactory

      This is the unspecialized class, which does nothing.  All work is
      done in the specializations for the different grid types.
   */
  template <class GridType>
  class GridFactory : public GridFactoryInterface<GridType> {

    /** \brief The grid world dimension */
    enum {dimworld = GridType::dimensionworld};

    /** \brief Type used by the grid for coordinates */
    typedef typename GridType::ctype ctype;

  public:

    /** \brief Default constructor */
    GridFactory();

    /** \brief Constructor for a given grid object

       If you already have your grid object constructed you can
       hand it over using this constructor.  A reason may be that
       you have need a UGGrid object with a non-default heap size.

       If you construct your factory class using this constructor
       the pointer handed over to you by the method createGrid() is
       the one you supplied here.
     */
    GridFactory(GridType* grid)
    {
      DUNE_THROW(GridError, "There is no grid factory for this grid type!");
    }


    /** \brief Insert a vertex into the coarse grid */
    virtual void insertVertex(const FieldVector<ctype,dimworld>& pos) {
      DUNE_THROW(GridError, "There is no grid factory for this grid type!");
    }

    /** \brief Insert an element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
     */
    virtual void insertElement(const GeometryType& type,
                               const std::vector<unsigned int>& vertices) {
      DUNE_THROW(GridError, "There is no grid factory for this grid type!");
    }

    /** \brief Finalize grid creation and hand over the grid

       The receiver takes responsibility of the memory allocated for the grid
     */
    virtual GridType* createGrid() {
      DUNE_THROW(GridError, "There is no grid factory for this grid type!");
    }

  };

}

#endif
