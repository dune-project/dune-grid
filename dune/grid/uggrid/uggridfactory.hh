// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_UGGRID_FACTORY_HH
#define DUNE_UGGRID_FACTORY_HH

/** \file
    \brief The specialization of the generic GridFactory for UGGrid
    \author Oliver Sander
 */

#include <array>
#include <memory>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/uggrid.hh>

namespace Dune {


  /** \brief Specialization of the generic GridFactory for UGGrid

      \ingroup GridFactory

     <p>
     If you want to write a routine that reads a grid from some
     file into a Dune UGGrid object you have to know how to use the UGGrid
     grid factory.  In the following we
     assume that you have a grid in some file format and an
     empty UGGrid object, created by one of its constructors.
     Hence, your file reader method signature may look like this:
     </p>

     <pre>
     UGGrid&lt;3&gt;* readMyFileFormat(const std::string& filename)
     </pre>

     Now, in order to create a valid UGGrid object do the
     following steps:

     <h2> 1) Create a GridFactory Object </h2>
     <p>
     Get a new GridFactory object by calling
     <pre>
     GridFactory<UGGrid<dim> > factory;
     </pre>

     <h2> 2)  Enter the Vertices </h2>

     <p>
     Insert the grid vertices by calling
     </p>

     <pre>
     factory.insertVertex(const FieldVector&lt;double,dimworld&gt;& position);
     </pre>

     <p>
     for each vertex.  The order of insertion determines the level- and leaf indices
     of your level 0 vertices.
     </p>


     <h2> 3) Enter the elements </h2>

     <p>
     For each element call
     </p>

     <pre>
     factory.insertElement(Dune::GeometryType type, const std::vector&lt;int&gt;& vertices);
     </pre>

     <p>
     The parameters are
     </p>

     <ul>
     <li> <b>type</b> - The element type.  UG supports the types <i>simplex</i> and
     <i>cube</i> in 2d, and <i>simplex, cube, prism</i>, and <i>pyramid</i> in 3d.
     <li> <b>vertices</b> - The Ids of the vertices of this element.</li>
     </ul>

     <p>
     The numbering of the vertices of each element is expected to follow the DUNE conventions.
     Refer to the page on reference elements for the details.

     <h2> 4) Parametrized Domains </h2>

     <p>
     UGGrid supports parametrized domains.  That means that you can provide a
     smooth description of your grid boundary.  The actual grid will always
     be piecewise linear; however, as you refine, the grid will approach your
     prescribed boundary.  You don't have to do this.  If your
     coarse grid boundary describes your domain well read on at Section 5.
     </p>

     <p>
     In order to create curved boundary segments, for each segment you have to write
     a class which implements the correct geometry.  These classes are then handed
     over to the UGGrid object.  Boundary segment implementations must be derived
     from
     <pre>
     template &lt;int dimworld&gt; Dune::BoundarySegment
     </pre>
     This is an abstract base class which requires you to overload the method

     <pre>
     virtual FieldVector&lt; double, dimworld &gt; operator() (const FieldVector&lt; double, dimworld-1 &gt; &local)
     </pre>

     <p>
     This methods must compute the world coordinates from the local ones on the
     boundary segment.  Give these classes to your grid factory by calling
     </p>
     <pre>
     factory.insertBoundarySegment(const std::vector&lt;int&gt;& vertices,
                             const BoundarySegment&lt;dimworld&gt; *boundarySegment);
     </pre>

     <p>
     Control over the allocated objects is taken from you, and the grid object
     will take care of their destruction.
     </p>

     <h2> 5) Finish construction </h2>

     <p>
     To finish off the construction of the UGGrid object call
     </p>

     <pre>
     std::unique_ptr<UGGrid<dim> > grid = factory.createGrid();
     </pre>

     <p>
     This time it is you who gets full responsibility for the allocated object.
     </p>

     <h2> Loading a Grid on a Parallel Machine </h2>
     <p>
     If you're working on a parallel machine, and you want to set up a
     parallel grid, proceed as described only on the rank-0 process.
     On the other processes just create a GridFactory and call createGrid()
     to obtain the grid object. This will create the grid on the master process
     and set up UG correctly on all other process.  Call <tt>loadBalance()</tt>
     to actually distribute the grid.
     </p>

     <p>\warning To use a parametrized boundary on a parallel machine you need
     to hand over the boundary segments to the grid factory on <b>all</b> processes.
     This behavior violates the Dune grid interface specification and will be
     corrected in the future.
     </p>
   */
  template <int dimworld>
  class GridFactory<UGGrid<dimworld> > : public GridFactoryInterface<UGGrid<dimworld> > {

    /** \brief Type used by the grid for coordinates */
    typedef typename UGGrid<dimworld>::ctype ctype;

    // UGGrid only in 2d and 3d
    static_assert(dimworld==2 || dimworld == 3, "UGGrid only in 2d and 3d");

  public:

    /** \brief Default constructor */
    GridFactory();

    /** \brief Constructor for a given grid object

       If you already have your grid object constructed you can
       hand it over using this constructor.

       If you construct your factory class using this constructor
       the pointer handed over to you by the method createGrid() is
       the one you supplied here.
     */
    GridFactory(UGGrid<dimworld>* grid);

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

    /** \brief Method to insert a boundary segment into a coarse grid

       Using this method is optional.  It only influences the ordering of the segments

        \param vertices The indices of the vertices of the segment
     */
    void insertBoundarySegment(const std::vector<unsigned int>& vertices);

    /** \brief Method to insert an arbitrarily shaped boundary segment into a coarse grid
        \param vertices The indices of the vertices of the segment
        \param boundarySegment Class implementing the geometry of the boundary segment.
     */
    void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                               const std::shared_ptr<BoundarySegment<dimworld> > &boundarySegment);


    /** \brief Finalize grid creation and hand over the grid

       The receiver takes responsibility of the memory allocated for the grid
     */
    virtual std::unique_ptr<UGGrid<dimworld>> createGrid();

    static const int dimension = UGGrid<dimworld>::dimension;

    template< int codim >
    struct Codim
    {
      typedef typename UGGrid<dimworld>::template Codim< codim >::Entity Entity;
    };

    /** \brief Return the number of the element in the order of insertion into the factory
     *
     * For UGGrid elements this number is the same as the element level index
     */
    virtual unsigned int
    insertionIndex ( const typename Codim< 0 >::Entity &entity ) const
    {
      return UG_NS<dimension>::levelIndex(entity.impl().target_);
    }

    /** \brief Return the number of the vertex in the order of insertion into the factory
     *
     * For UGGrid vertices this number is the same as the vertex level index
     */
    virtual unsigned int
    insertionIndex ( const typename Codim< dimension >::Entity &entity ) const
    {
      return UG_NS<dimension>::levelIndex(entity.impl().target_);
    }

    /** \brief Return the number of the intersection in the order of insertion into the factory
     *
     * For UGGrid intersections this number is the same as the boundary segment index
     */
    virtual unsigned int
    insertionIndex ( const typename UGGrid<dimworld>::LeafIntersection &intersection ) const
    {
      return intersection.boundarySegmentIndex();
    }

    /** \brief Return true if the intersection has been explictily insterted into the factory */
    virtual bool
    wasInserted ( const typename UGGrid<dimworld>::LeafIntersection &intersection ) const
    {
      return (insertionIndex( intersection ) < boundarySegmentVertices_.size());
    }

    using Communication = typename UGGrid<dimworld>::Communication;

    /** \brief Return the Communication used by the grid factory
     *
     * Use the Communication available from the grid.
     */
    Communication comm() const
    {
        return grid_->comm();
    }

  private:

    // Initialize the grid structure in UG
    void createBegin();

    // Pointer to the grid being built
    UGGrid<dimworld>* grid_;

    // True if the factory allocated the grid itself, false if the
    // grid was handed over from the outside
    bool factoryOwnsGrid_;

    /** \brief Buffer for the vertices of each explicitly given boundary segment */
    std::vector<std::array<int, dimworld*2-2> > boundarySegmentVertices_;

    /** \brief While inserting the elements this array records the number of
        vertices of each element. */
    std::vector<unsigned char> elementTypes_;

    /** \brief While inserting the elements this array records the vertices
        of the elements. */
    std::vector<unsigned int> elementVertices_;

    /** \brief Buffer the vertices until createend() is called */
    std::vector<FieldVector<double, dimworld> > vertexPositions_;

  };

}

#endif
