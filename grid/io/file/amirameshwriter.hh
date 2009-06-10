// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_AMIRAMESH_WRITER_HH
#define DUNE_AMIRAMESH_WRITER_HH

#include <string>

#include <dune/common/array.hh>
#include <amiramesh/AmiraMesh.h>


namespace Dune {

  /** @ingroup AmiraMesh
   * \brief Provides file writing facilities in the AmiraMesh format.
   *
   */
  template<class GridView>
  class AmiraMeshWriter {

    enum {dim = GridView::dimension};

  public:

    /** \brief Add a grid view to the file
        \param grid GridView to be written
        \param splitQuads Amira doesn't support 2d quad grids.  If this is set any
        quadrilateral will be split in two triangles.  If not, the file is not
        readable by standard Amira.
     */
    void addGrid(const GridView& gridView, bool splitAll=false);

    /** \brief Add level grid
        \param grid Grid to be written
        \param indexSet Level of the level grid that is to be written
        \param splitQuads Amira doesn't support 2d quad grids.  If this is set any
        quadrilateral will be split in two triangles.  If not, the file is not
        readable by standard Amira.
     */
    template <class GridType2>
    void addLevelGrid(const GridType2& grid, int level, bool splitAll=false);

    /** \brief Add leaf grid
        \param grid Grid to be written
        \param splitQuads Amira doesn't support 2d quad grids.  If this is set any
        quadrilateral will be split in two triangles.  If not, the file is not
        readable by standard Amira.
     */
    template <class GridType2>
    void addLeafGrid(const GridType2& grid, bool splitAll=false);

    /** \brief Add cell data
        \param An ISTL compliant vector type
        \param Grid view that the data belongs to
     */
    template <class DataContainer>
    void addCellData(const DataContainer& data, const GridView& gridView);

    /** \brief Add vertex data
        \param An ISTL compliant vector type
        \param Grid view that the data belongs to
     */
    template <class DataContainer>
    void addVertexData(const DataContainer& data, const GridView& gridView, bool GridSplitUp=false);

    /** \brief Write AmiraMesh object to disk
        \param filename Name of the file to write to
        \param ascii Set this if you want an ascii AmiraMesh file
     */
    void write(const std::string& filename, bool ascii=false) const;

    /** \brief Write data on a uniform grid into an AmiraMesh file
     */
    template <class DataContainer>
    void addUniformData(const GridView& gridView,
                        const array<unsigned int, dim>& n,
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
    : public AmiraMeshWriter<typename GridType::LevelGridView>
  {

  public:

    /** \brief Default constructor */
    LevelAmiraMeshWriter() {}

    /** \brief Constructor which initializes the AmiraMesh object with a given level grid */
    LevelAmiraMeshWriter(const GridType& grid, int level) {
      this->addGrid(grid.levelView(level));
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
                                 int level,
                                 bool GridSplitUp=false) {

      LevelAmiraMeshWriter amiramesh;
      amiramesh.addVertexData(f, grid.levelView(level),GridSplitUp);
      amiramesh.write(filename);
    }

  };

  /** @ingroup AmiraMesh
   * \brief Provides file writing facilities in the AmiraMesh format for leaf grids.
   *
   */
  template<class GridType>
  class LeafAmiraMeshWriter
    : public AmiraMeshWriter<typename GridType::LeafGridView>
  {

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
                                 const std::string& filename,
                                 bool GridSplitUp = false) {
      LeafAmiraMeshWriter amiramesh;
      amiramesh.addVertexData(f, grid.leafView(),GridSplitUp);
      amiramesh.write(filename);
    }

  };

}

// implementation
#include "amiramesh/amirameshwriter.cc"

#endif
