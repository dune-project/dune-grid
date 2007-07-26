// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_AMIRAMESH_READER_HH
#define DUNE_AMIRAMESH_READER_HH

#include <string>
#include <dune/common/array.hh>
#include <amiramesh/AmiraMesh.h>

namespace Dune {

  /** @ingroup AmiraMesh
   * \brief Provides file reading facilities in the AmiraMesh format.
   *
   */
  template<class GridType>
  class AmiraMeshReader {

  public:

    /** \brief The method that does the reading.
     *
     * @param grid The grid objects that is to be written
     * @param filename The filename
     */
    static void read(GridType& grid,
                     const std::string& filename);

    /** \brief Read a block vector from an AmiraMesh file
     *
     * \param f The vector to read into.  Implicitly assumed to be an ISTL vector
     * \param filename Name of the AmiraMesh file
     */
    template<class DiscFuncType>
    static void readFunction(DiscFuncType& f, const std::string& filename);

    /** \brief Dummy constructor */
    AmiraMeshReader() {}

  };

}

// Default implementation
template<class GridType>
void Dune::AmiraMeshReader<GridType>::read(GridType& grid,
                                           const std::string& filename)
{
  DUNE_THROW(IOError, "No AmiraMesh reading has been implemented for " << grid.name() << "!");
}

#include "amiramesh/amirameshreader.cc"

// the amiramesh reader for UGGrid
#ifdef HAVE_UG
#include "amiramesh/amuggridreader.hh"
#endif

#endif
