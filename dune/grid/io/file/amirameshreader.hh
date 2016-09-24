// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_AMIRAMESH_READER_HH
#define DUNE_AMIRAMESH_READER_HH

#include <memory>
#include <string>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/gridreader.hh>

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
  class AmiraMeshReader
      : public GridReader<GridType, AmiraMeshReader<GridType>>
  {
    /** \brief Dimension of the grid */
    enum {dim = GridType::dimension};

    /** \brief Create the boundary description from an explicitly given psurface object */
    static void createDomain(GridFactory<GridType>& factory, const std::shared_ptr<PSurfaceBoundary<dim-1> >& boundary);

    /** \brief Create the actual grid */
    static void buildGrid(GridFactory<GridType>& factory, AmiraMesh* am);

    /** \brief Create the actual grid */
    static void build2dGrid(GridFactory<GridType>& factory, AmiraMesh* am);

  protected:

    /** \brief The method that does the reading.
     *
     * @param factory  A grid-factory that does the final construction of the grid
     * @param filename The filename
     *
     * Implementation of \ref GridReader::read.
     */
    static void readFactoryImp(GridFactory<GridType>& factory, const std::string& filename);


    /** \brief Read a grid with a parametrized boundary

       Several grid managers support parametrized boundary segments, which carry
       functions describing the true shape of the boundary segment.
       This information will the be considered when refining the grid.

       @param factory  A grid-factory that does the final construction of the grid
       @param filename The name of the grid file
       @param boundary Pointer to an object holding the description of the grid domain boundary

       Implementation of \ref GridReader::read.
     */
    static void readFactoryImp(GridFactory<GridType>& factory, const std::string& filename,
                               const std::shared_ptr<PSurfaceBoundary<dim-1> >& boundary);

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
