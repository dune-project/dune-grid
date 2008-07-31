// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_AMIRAMESH_WRITER_HH
#define DUNE_AMIRAMESH_WRITER_HH

#include <string>

#include <dune/common/fixedarray.hh>
#include <amiramesh/AmiraMesh.h>


namespace Dune {

  /** @ingroup AmiraMesh
   * \brief Provides file writing facilities in the AmiraMesh format.
   *
   */
  template<class GridType, class IndexSetType>
  class AmiraMeshWriter {

    enum {dim = GridType::dimension};

  public:

    /** \brief Add a grid view to the file
        \param grid GridView to be written
        \param splitQuads Amira doesn't support 2d quad grids.  If this is set any
        quadrilateral will be split in two triangles.  If not, the file is not
        readable by standard Amira.
     */
    template <class GridView>
    void addGrid(const GridView& gridView, bool splitQuads=false);

    /** \brief Add level grid
        \param grid Grid to be written
        \param indexSet Level of the level grid that is to be written
        \param splitQuads Amira doesn't support 2d quad grids.  If this is set any
        quadrilateral will be split in two triangles.  If not, the file is not
        readable by standard Amira.
     */
    template <class GridType2>
    void addLevelGrid(const GridType2& grid, int level, bool splitQuads=false);

    /** \brief Add leaf grid
        \param grid Grid to be written
        \param splitQuads Amira doesn't support 2d quad grids.  If this is set any
        quadrilateral will be split in two triangles.  If not, the file is not
        readable by standard Amira.
     */
    template <class GridType2>
    void addLeafGrid(const GridType2& grid, bool splitQuads=false);

    /** \brief Add cell data */
    template <class GridType2, class DataContainer>
    void addCellData(const DataContainer& data, const GridType2& grid);

    /** \brief Add vertex data
        \param An ISTL compliant vector type
        \param IndexSet of the grid that the data belongs to
     */
    template <class DataContainer>
    void addVertexData(const DataContainer& data, const IndexSetType& indexSet);

    /** \brief Write AmiraMesh object to disk
        \param filename Name of the file to write to
        \param ascii Set this if you want an ascii AmiraMesh file
     */
    void write(const std::string& filename, bool ascii=false) const;

    /** \brief Write data on a uniform grid into an AmiraMesh file
     */
    template <class GridType2, class DataContainer>
    void addUniformData(const GridType2& grid,
                        const array<unsigned int, GridType2::dimension>& n,
                        const DataContainer& data);

  protected:

    AmiraMesh amiramesh_;
  };

  /** @ingroup AmiraMesh
   * \brief Provides file writing facilities in the AmiraMesh format for level grids.
   *
   */
  template<class GridType>
  class LevelAmiraMeshWriter
    : public AmiraMeshWriter<GridType,typename GridType::Traits::LevelIndexSet>
  {

    enum {dim = GridType::dimension};

    typedef AmiraMeshWriter<GridType,typename GridType::Traits::LevelIndexSet> WriterType;

  public:

    /** \brief Default constructor */
    LevelAmiraMeshWriter() {}

    /** \brief Constructor which initializes the AmiraMesh object with a given level grid */
    LevelAmiraMeshWriter(const GridType& grid, int level) {
      this->addLevelGrid(grid, level);
    }

    /** \brief Write a grid in AmiraMesh format
     *
     * @param grid The grid objects that is to be written
     * @param filename The filename
     * @param level The level to be written
     */
    static void writeGrid(const GridType& grid,
                          const std::string& filename,
                          int level) {
      LevelAmiraMeshWriter amiramesh(grid, level);
      amiramesh.write(filename);
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

      LevelAmiraMeshWriter amiramesh;
      amiramesh.addVertexData(f, grid.levelIndexSet(level));
      amiramesh.write(filename);
    }

  };

  /** @ingroup AmiraMesh
   * \brief Provides file writing facilities in the AmiraMesh format for leaf grids.
   *
   */
  template<class GridType>
  class LeafAmiraMeshWriter
    : public AmiraMeshWriter<GridType,typename GridType::Traits::LeafIndexSet>
  {

    enum {dim = GridType::dimension};

    typedef AmiraMeshWriter<GridType,typename GridType::Traits::LeafIndexSet> WriterType;

  public:

    /** \brief Default constructor */
    LeafAmiraMeshWriter() {}

    /** \brief Constructor which initializes the AmiraMesh object with a given leaf grid */
    LeafAmiraMeshWriter(const GridType& grid) {
      this->addLeafGrid(grid);
    }

    /** \brief Write a grid in AmiraMesh format
     *
     * @param grid The grid objects that is to be written
     * @param filename The filename
     */
    static void writeGrid(const GridType& grid,
                          const std::string& filename) {
      LeafAmiraMeshWriter amiramesh(grid);
      amiramesh.write(filename);
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
      LeafAmiraMeshWriter amiramesh;
      amiramesh.addVertexData(f, grid.leafIndexSet());
      amiramesh.write(filename);
    }

  };

}

// The default implementation
#include "amiramesh/amirameshwriter.cc"

#endif
