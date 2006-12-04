// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_AMIRAMESH_WRITER_HH
#define DUNE_AMIRAMESH_WRITER_HH

#include <string>

namespace Dune {

  /** @ingroup IO
   * \brief Provides file writing facilities in the AmiraMesh format.
   *
   */
  template<class GridType, class IndexSetType>
  class AmiraMeshWriter {

    enum {dim = GridType::dimension};

  public:

    /** \brief Write a grid in AmiraMesh format
     *
     * @param grid The grid objects that is to be written
     * @param filename The filename
     * @param indexSet The index set encodes whether the leaf grid or a level grid is written
     */
    static void writeGrid(const GridType& grid,
                          const std::string& filename,
                          const IndexSetType& indexSet);

    /** \brief Writes an ISTL block vector in AmiraMesh format
        @param filename The filename
        @param indexSet The index set encodes whether the vector lives of the
        leaf grid or on a level grid
     */
    template <class VectorType>
    static void writeBlockVector(const VectorType& f,
                                 const std::string& filename,
                                 const IndexSetType& indexSet);

  };

  /** @ingroup IO
   * \brief Provides file writing facilities in the AmiraMesh format for level grids.
   *
   */
  template<class GridType>
  class LevelAmiraMeshWriter {

    enum {dim = GridType::dimension};

    typedef AmiraMeshWriter<GridType,typename GridType::Traits::LevelIndexSet> WriterType;

  public:

    /** \brief Write a grid in AmiraMesh format
     *
     * @param grid The grid objects that is to be written
     * @param filename The filename
     * @param level The level to be written
     */
    static void writeGrid(const GridType& grid,
                          const std::string& filename,
                          int level) {
      WriterType::writeGrid(grid,
                            filename,
                            grid.levelIndexSet(level));
    }

    /** \brief Writes an ISTL block vector in AmiraMesh format
        @param grid The grid objects that the vector lives on
        @param f The vector to be written.  Has to comply with the ISTL conventions
        @param filename The filename
        @param level The level of the grid that the vector lives on
     */
    template <class VectorType>
    static void writeBlockVector(const GridType& grid,
                                 const VectorType& f,
                                 const std::string& filename,
                                 int level) {

      WriterType::template writeBlockVector<VectorType>(f, filename, grid.levelIndexSet(level));
    }

  };

  /** @ingroup IO
   * \brief Provides file writing facilities in the AmiraMesh format for leaf grids.
   *
   */
  template<class GridType>
  class LeafAmiraMeshWriter {

    enum {dim = GridType::dimension};

    typedef AmiraMeshWriter<GridType,typename GridType::Traits::LeafIndexSet> WriterType;

  public:

    /** \brief Write a grid in AmiraMesh format
     *
     * @param grid The grid objects that is to be written
     * @param filename The filename
     */
    static void writeGrid(const GridType& grid,
                          const std::string& filename) {
      WriterType::writeGrid(grid, filename, grid.leafIndexSet());
    }

    /** \brief Writes an ISTL block vector in AmiraMesh format
        @param grid The grid objects that the vector lives on
        @param f The vector to be written.  Has to comply with the ISTL conventions
        @param filename The filename
     */
    template <class VectorType>
    static void writeBlockVector(const GridType& grid,
                                 const VectorType& f,
                                 const std::string& filename) {
      WriterType::template writeBlockVector<VectorType>(f, filename, grid.leafIndexSet());
    }

  };

}

// The default implementation
#include "amiramesh/amirameshwriter.cc"

// the amiramesh writer for SGrid
//#include "amiramesh/amsgridwriter.cc"

#endif
