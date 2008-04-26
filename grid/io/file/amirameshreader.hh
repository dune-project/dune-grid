// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_AMIRAMESH_READER_HH
#define DUNE_AMIRAMESH_READER_HH

#include <string>

#include <amiramesh/AmiraMesh.h>

namespace Dune {

  /** @ingroup AmiraMesh
   * \brief Provides file reading facilities in the AmiraMesh format.
   *
   */
  template<class GridType>
  class AmiraMeshReader {

    /** \brief Create the boundary description from an explicitly given psurface file */
    static void createDomain(GridType& grid, const std::string& filename);

    /** \brief Create the actual grid */
    static void buildGrid(GridType& grid, AmiraMesh* am);
  public:

    /** \brief The method that does the reading.
     *
     * @param grid The grid objects that is to be read
     * @param filename The filename
     */
    static void read(GridType& grid,
                     const std::string& filename);

    /** \brief Read a grid with a parametrized boundary

       Several grid managers support parametrized boundary segment which carry
       function describing the true shape of the boundary segment.
       This information will the be considered when refining the grid.

       In
       <em>'Krause, Sander, Automatic Construction of Boundary Parametrizations
       for Geometric Multigrid Solvers, CVS, 2005'</em>,
       the authors describe a
       way to automatically build such boundary descriptions.  Their
       file format can be read by this routine.

       To use this feature you
       have to have the psurface library and build Dune with --with-psurface.
       Ask Oliver sander@mi.fu-berlin.de for help.

       @param grid The grid objects that is to be read
       @param filename The name of the grid file
       @param domainFilename The name of the psurface boundary file
     */
    static void read(GridType& grid,
                     const std::string& filename,
                     const std::string& domainFilename);

    /** \brief Read a block vector from an AmiraMesh file
     *
     * \param f The vector to read into.  Implicitly assumed to be an ISTL vector
     * \param filename Name of the AmiraMesh file
     */
    template<class DiscFuncType>
    static void readFunction(DiscFuncType& f, const std::string& filename);

  };

}

#include "amiramesh/amirameshreader.cc"

#endif
