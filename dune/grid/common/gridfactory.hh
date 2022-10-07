// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_GRIDFACTORY_HH
#define DUNE_GRID_COMMON_GRIDFACTORY_HH

/** \file
    \brief Provide a generic factory class for unstructured grids.
 */

#include <memory>
#include <vector>

#include <dune/common/deprecated.hh>
#define DUNE_FUNCTION_HH_SILENCE_DEPRECATION
#include <dune/common/function.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

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
     Dune::GridFactory<Grid> factory;

     factory.insertVertex({0, 0, 0});
     factory.insertVertex({1, 0, 0});
     factory.insertVertex({0, 1, 0});
     factory.insertVertex({1, 1, 0});
     factory.insertVertex({0, 0, 1});
     factory.insertVertex({1, 0, 1});
     factory.insertVertex({0, 1, 1});
     factory.insertVertex({1, 1, 1});

     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 1, 3, 7});
     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 5, 1, 7});
     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 4, 5, 7});
     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 6, 4, 7});
     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 2, 6, 7});
     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 3, 2, 7});

     std::shared_ptr<Grid> gridp(factory.createGrid());
     \endcode
   * Make sure that the inserted elements are not inverted, since not all
   * grids support that.  For instance, in the following code snippet the
   * elements 1, 3 and 5 are inverted while elements 0, 2 and 4 are not.
     \code
     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 1, 3, 7});
     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 1, 5, 7});
     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 4, 5, 7});
     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 4, 6, 7});
     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 2, 6, 7});
     factory.insertElement(Dune::GeometryTypes::tetrahedron, {0, 2, 3, 7});
     \endcode
   */
  template <class GridType>
  class GridFactoryInterface
  {

  protected:
    /** \brief dimension of the grid */
    static const int dimension = GridType::dimension;

    /** \brief The grid world dimension */
    constexpr static int dimworld = GridType::dimensionworld;

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

    DUNE_NO_DEPRECATED_BEGIN
    /** \brief Insert a parametrized element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \param elementParametrization A function prescribing the shape of this element

        Make sure the inserted element is not inverted (this holds even
        for simplices).  There are grids that can't handle inverted elements.

        \deprecated [After Dune 2.7] VirtualFunction is deprecated, use the
                    overload taking a std::function instead
     */
    [[deprecated("[After Dune 2.7]: VirtualFunction is deprecated, use the "
                 "overload taking a std::function instead")]]
    virtual void
    insertElement([[maybe_unused]] const GeometryType& type,
                  [[maybe_unused]] const std::vector<unsigned int>& vertices,
                  [[maybe_unused]] const std::shared_ptr<VirtualFunction<
                                         FieldVector<ctype,dimension>,
                                         FieldVector<ctype,dimworld>
                         > >& elementParametrization)
    {
      DUNE_THROW(GridError, "This grid does not support parametrized elements!");
    }
    DUNE_NO_DEPRECATED_END

    /** \brief Insert a parametrized element into the coarse grid

        \param type                   The GeometryType of the new element
        \param vertices               The vertices of the new element, using
                                      the DUNE numbering
        \param elementParametrization A function prescribing the shape of this
                                      element

        Make sure the inserted element is not inverted (this holds even for
        simplices).  There are grids that can't handle inverted elements.
     */
    virtual void
    insertElement(const GeometryType& type,
                  const std::vector<unsigned int>& vertices,
                  std::function<FieldVector<ctype,dimworld>
                                  (FieldVector<ctype,dimension>)>
                       elementParametrization)
    {
      // note: this forward to the overload taking a Virtual function during
      // the deprecation period, once that is over it should the throwing of
      // the exception should be moved here directly
      using Domain = FieldVector<ctype,dimension>;
      using Range = FieldVector<ctype,dimworld>;
      DUNE_NO_DEPRECATED_BEGIN
      auto f =
        makeVirtualFunction<Domain, Range>(std::move(elementParametrization));
      insertElement(type, vertices,
                    std::make_unique<decltype(f)>(std::move(f)));
      DUNE_NO_DEPRECATED_END
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
    virtual void insertBoundarySegment([[maybe_unused]] const std::vector<unsigned int>& vertices,
                                       [[maybe_unused]] const std::shared_ptr<BoundarySegment<dimension,dimworld> >& boundarySegment)
    {
      DUNE_THROW(GridError, "This grid does not support parametrized boundary segments!");
    }

    /** \brief Finalize grid creation and hand over the grid

       The receiver takes responsibility of the memory allocated for the grid
     */
    virtual std::unique_ptr<GridType> createGrid() = 0;

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
    insertionIndex ( [[maybe_unused]] const typename Codim< 0 >::Entity &entity ) const
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
    insertionIndex ( [[maybe_unused]] const typename Codim< dimension >::Entity &entity ) const
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
    insertionIndex ( [[maybe_unused]] const typename GridType::LeafIntersection &intersection ) const
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
    wasInserted ( [[maybe_unused]] const typename GridType::LeafIntersection &intersection ) const
    {
      DUNE_THROW( NotImplemented, "insertion indices have not yet been implemented." );
    }

    using Communication = Dune::Communication<typename MPIHelper::MPICommunicator>;

    /** \brief Return the Communication used by the grid factory
     *
     * Defaults to the Communication induced by the process-local communicator.
     */
    Communication comm() const
    {
      return Communication(MPIHelper::getLocalCommunicator());
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

    typedef GridFactoryInterface<GridType> Base;

    /** \brief The grid world dimension */
    constexpr static int dimworld = GridType::dimensionworld;

    /** \brief Type used by the grid for coordinates */
    typedef typename GridType::ctype ctype;

  public:

    // use default implementation from base class
    using Base::insertBoundarySegment;

    /** \brief Default constructor */
    GridFactory() {
      DUNE_THROW(GridError, "There is no grid factory for this grid type!");
    }

    /** \brief Insert a vertex into the coarse grid */
    virtual void insertVertex([[maybe_unused]] const FieldVector<ctype,dimworld>& pos) {
      DUNE_THROW(GridError, "There is no grid factory for this grid type!");
    }

    /** \brief Insert an element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering

        Make sure the inserted element is not inverted (this holds even
        for simplices).  There are grids that can't handle inverted tets.
     */
    virtual void insertElement([[maybe_unused]] const GeometryType& type,
                               [[maybe_unused]] const std::vector<unsigned int>& vertices) {
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
    virtual void insertBoundarySegment([[maybe_unused]] const std::vector<unsigned int>& vertices) {
      DUNE_THROW(GridError, "There is no grid factory for this grid type!");
    }

    /** \brief Finalize grid creation and hand over the grid

       The receiver takes responsibility of the memory allocated for the grid
     */
    virtual std::unique_ptr<GridType> createGrid() {
      DUNE_THROW(GridError, "There is no grid factory for this grid type!");
    }

  };

}

#endif
