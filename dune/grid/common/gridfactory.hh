// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_FACTORY_HH
#define DUNE_GRID_FACTORY_HH

/** \file
    \brief Provide a generic factory class for unstructured grids.
 */

#include <memory>
#include <vector>

#include <dune/common/function.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/grid.hh>

namespace Dune
{

  /** \brief Provide a generic factory class for unstructured grids.
   *
   * \ingroup GridFactory
   *
   * This base class declares the interface.
   *
   * Example use: create a grid consisting of a cube triangulated into 6
   * tetrahedra:
   * \image html cube-to-tet-6.png "Left: cube with vertex numbers.  Middle: cube triangulated into six tetrahedra.  Right: exploded version of the middle figure, with number for the tetrahedra."
     \code
     Dune::GridFactory<Grid> gf;

     Dune::FieldVector<typename Grid::ctype, 3> pos;

     pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
     pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
     pos[0] = 0; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
     pos[0] = 1; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
     pos[0] = 0; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
     pos[0] = 1; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
     pos[0] = 0; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);
     pos[0] = 1; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);

     Dune::GeometryType type;
     type.makeTetrahedron();
     std::vector<unsigned int> vid(4);

     vid[0] = 0; vid[1] = 1; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
     vid[0] = 0; vid[1] = 5; vid[2] = 1; vid[3] = 7; gf.insertElement(type, vid);
     vid[0] = 0; vid[1] = 4; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
     vid[0] = 0; vid[1] = 6; vid[2] = 4; vid[3] = 7; gf.insertElement(type, vid);
     vid[0] = 0; vid[1] = 2; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);
     vid[0] = 0; vid[1] = 3; vid[2] = 2; vid[3] = 7; gf.insertElement(type, vid);

     std::shared_ptr<Grid> gridp(gf.createGrid());
     \endcode
   * Make sure that the inserted elements are not inverted, since not all
   * grids support that.  For instance, in the following code snippet the
   * elements 1, 3 and 5 are inverted while elements 0, 2 and 4 are not.
     \code
     vid[0] = 0; vid[1] = 1; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
     vid[0] = 0; vid[1] = 1; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
     vid[0] = 0; vid[1] = 4; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
     vid[0] = 0; vid[1] = 4; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);
     vid[0] = 0; vid[1] = 2; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);
     vid[0] = 0; vid[1] = 2; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
     \endcode
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

    /** \brief virtual destructor */
    virtual ~GridFactoryInterface ()
    {}

    /** \brief Insert a vertex into the coarse grid */
    virtual void insertVertex(const FieldVector<ctype,dimworld>& pos) = 0;

    /** \brief Insert an element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering

        Make sure the inserted element is not inverted (this holds even
        for simplices).  There are grids that can't handle inverted elements.
     */
    virtual void insertElement(const GeometryType& type,
                               const std::vector<unsigned int>& vertices) = 0;

    /** \brief Insert a parametrized element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \param elementParametrization A function prescribing the shape of this element

        Make sure the inserted element is not inverted (this holds even
        for simplices).  There are grids that can't handle inverted elements.
     */
    virtual void insertElement(const GeometryType& type,
                               const std::vector<unsigned int>& vertices,
                               const std::shared_ptr<VirtualFunction<FieldVector<ctype,dimension>,FieldVector<ctype,dimworld> > >& elementParametrization)
    {
      DUNE_THROW(GridError, "This grid does not support parametrized elements!");
    }

    /** \brief insert a boundary segment
     *
     *  This method inserts a boundary segment into the coarse grid. Using
     *  this method has two advantages over not using it:
     *  - The boundary segment gets an insertion index.
     *  - The grid factory can verify that this is actually a boundary segment
     *  .
     *
     *  \note You are not forced to insert all boundary segments. The grid
     *        factory will find the remaining boundary segments itself.
     *
     *  \param[in]  vertices  the indices of the vertices of the segment
     */
    virtual void insertBoundarySegment(const std::vector<unsigned int>& vertices) = 0;

    /** \brief insert an arbitrarily shaped boundary segment
     *
     *  This method inserts a boundary segment into the coarse grid.
     *
     *  \param[in]  vertices         the indices of the vertices of the segment
     *  \param[in]  boundarySegment  user defined implementation of the boundary segment's geometry
     */
    virtual void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                                       const std::shared_ptr<BoundarySegment<dimension,dimworld> >& boundarySegment)
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
     *  This method allows checking whether an intersection was actually
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

      See GridFactoryInterface for an example how to use this class.
   */
  template <class GridType>
  class GridFactory : public GridFactoryInterface<GridType> {

    /** \brief The grid world dimension */
    enum {dimworld = GridType::dimensionworld};

    /** \brief Type used by the grid for coordinates */
    typedef typename GridType::ctype ctype;

  public:

    /** \brief Default constructor */
    GridFactory() {
      DUNE_THROW(GridError, "There is no grid factory for this grid type!");
    }

    /** \brief Insert a vertex into the coarse grid */
    virtual void insertVertex(const FieldVector<ctype,dimworld>& pos) {
      DUNE_THROW(GridError, "There is no grid factory for this grid type!");
    }

    /** \brief Insert an element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering

        Make sure the inserted element is not inverted (this holds even
        for simplices).  There are grids that can't handle inverted tets.
     */
    virtual void insertElement(const GeometryType& type,
                               const std::vector<unsigned int>& vertices) {
      DUNE_THROW(GridError, "There is no grid factory for this grid type!");
    }

    /** \brief insert a boundary segment
     *
     *  This method inserts a boundary segment into the coarse grid. Using
     *  this method has two advantages over not using it:
     *  - The boundary segment gets an insertion index.
     *  - The grid factory can verify that this is actually a boundary segment
     *  .
     *
     *  \note You are not forced to insert all boundary segments. The grid
     *        factory will find the remaining boundary segments itself.
     *
     *  \param[in]  vertices  the indices of the vertices of the segment
     */
    virtual void insertBoundarySegment(const std::vector<unsigned int>& vertices) {
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
