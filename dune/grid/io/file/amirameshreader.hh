// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_AMIRAMESH_READER_HH
#define DUNE_AMIRAMESH_READER_HH

#warning Support for AmiraMesh is deprecated and will be removed after Dune 2.8.

#include <memory>
#include <string>

#include <dune/common/to_unique_ptr.hh>
#include <dune/grid/common/gridfactory.hh>

#if HAVE_PSURFACE
#include <dune/grid/io/file/amiramesh/psurfaceboundary.hh>

#if HAVE_AMIRAMESH
#include <amiramesh/AmiraMesh.h>
#else
// forward declaration so we can at least compile the header without libamiramesh
class AmiraMesh;
#endif

namespace Dune {

  /** @ingroup AmiraMesh
   * \brief Provides file reading facilities in the AmiraMesh format.
   *
   */
  template<class GridType>
  class AmiraMeshReader {

    using ctype = typename GridType::ctype;

    /** \brief Dimension of the grid */
    enum {dim = GridType::dimension};

    /** \brief Create the boundary description from an explicitly given psurface object */
    static void createDomain(GridFactory<GridType>& factory, const std::shared_ptr<PSurfaceBoundary<dim-1, ctype> >& boundary);

    /** \brief Create the actual grid */
    static void buildGrid(GridFactory<GridType>& factory, AmiraMesh* am);

    /** \brief Create the actual grid */
    static void build2dGrid(GridFactory<GridType>& factory, AmiraMesh* am);

  public:

    /** \brief The method that does the reading.
     *
     * @param filename The filename
     *
     * \return The return type is a special pointer type that casts into GridType*,
     *    std::unique_ptr<GridType>, and std::shared_ptr<GridType>.  It is scheduled
     *    to be replaced by std::unique_ptr<GridType> eventually.
     */
    static ToUniquePtr<GridType> read(const std::string& filename);

    /** \brief Read a grid from file into a given grid object
     *
     * @param grid The grid objects that is to be read
     * @param filename The filename
     */
    static void read(GridType& grid,
                     const std::string& filename);

    /** \brief Read a grid with a parametrized boundary

       Several grid managers support parametrized boundary segments, which carry
       functions describing the true shape of the boundary segment.
       This information will the be considered when refining the grid.

       @param filename The name of the grid file
       @param boundary Pointer to an object holding the description of the grid domain boundary

       \return The return type is a special pointer type that casts into GridType*,
          std::unique_ptr<GridType>, and std::shared_ptr<GridType>.  It is scheduled
          to be replaced by std::unique_ptr<GridType> eventually.
     */
    static ToUniquePtr<GridType> read(const std::string& filename,
                          const std::shared_ptr<PSurfaceBoundary<dim-1, ctype> >& boundary);

  private:
    /** \brief Read a grid with a parametrized boundary into a given grid object

       \note Reading files into existing grid objects is not officially supported,
       and this method is here only for a transition period.

       @param grid The grid objects that is to be read
       @param filename The name of the grid file
       @param boundary The PSurface boundary object
     */
    static void read(GridType& grid,
                     const std::string& filename,
                     const std::shared_ptr<PSurfaceBoundary<dim-1, ctype> >& boundary);

  public:
    /** \brief Read a block vector from an AmiraMesh file
     *
     * The data may either be given on the nodes (P1-Functions) or the elements (P0-Functions).
     *
     * \param f The vector to read into.  Implicitly assumed to be an ISTL vector
     * \param filename Name of the AmiraMesh file
     */
    template<class DiscFuncType>
    static void readFunction(DiscFuncType& f, const std::string& filename);

  };

}

#if HAVE_AMIRAMESH
#include "amiramesh/amirameshreader.cc"
#endif

#endif // #if HAVE_PSURFACE
#endif
