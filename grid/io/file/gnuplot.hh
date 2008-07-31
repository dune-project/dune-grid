// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IO_GNUPLOT_HH
#define DUNE_IO_GNUPLOT_HH

/** @file
    @author Christian Engwer
    @brief Provides gnuplot output for 1D Grids
 */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/helpertemplates.hh>

#include <dune/grid/common/grid.hh>

namespace Dune {

  /** \brief Writer for 1D grids in gnuplot format
      \ingroup Gnuplot
      \tparam GridType the grid
      \tparam GridView Level- or LeafGridView
   */
  template<class GridType, class GridView>
  class GnuplotWriter {

    typedef typename GridType::ctype ctype;

    enum {dimworld = GridType::dimensionworld};

  public:
    GnuplotWriter (const GridType & g, const GridView & gv) : _grid(g), _is(gv.indexSet()), _gv(gv)
    {
      dune_static_assert(dimworld==1 || dimworld==2, "GnuPlot export only works for worlddim==1 and worlddim==2");
      // allocate _data buffer
      _data.resize(_is.size(0)*2);
    }

    /** \brief Add cell data
        \param An ISTL compliant vector type
        \param name associated with the data
     */
    template <class DataContainer>
    void addCellData(const DataContainer& data, const std::string & name)
    {
      if (dimworld!=1)
        DUNE_THROW(IOError, "Gnuplot cell data writing is only supported for grids in a 1d world!");
      addData(cellData, data, name);
    }

    /** \brief Add vertex data
        \param An ISTL compliant vector type
        \param name associated with the data
     */
    template <class DataContainer>
    void addVertexData(const DataContainer& data, const std::string & name)
    {
      addData(vertexData, data, name);
    }

    /** \brief Write Gnuplot file to disk
        \param filename Name of the file to write to
     */
    void write(const std::string& filename) const;

  private:
    enum DataType { vertexData, cellData };
    const GridType & _grid;
    const typename GridView::IndexSet & _is;
    const GridView _gv;
    std::vector< std::vector< float > > _data;
    std::vector< std::string > _names;

    template <class DataContainer>
    void addData(DataType t, const DataContainer& data, const std::string & name);

    void writeRow(std::ostream & file,
                  const FieldVector<ctype, dimworld>& position,
                  const std::vector<float> & data) const;
  };

  /** \brief GnuplotWriter on the leaf grid
      \ingroup Gnuplot
   */
  template<class G>
  class LeafGnuplotWriter : public GnuplotWriter<G,typename G::template Codim<0>::LeafIndexSet>
  {
  public:
    /** \brief Construct a Gnuplot writer for the leaf level of a given grid */
    LeafGnuplotWriter (const G& grid)
      : GnuplotWriter<G,typename G::template Codim<0>::LeafIndexSet>(grid,grid.leafIndexSet())
    {}
  };

  /** \brief GnuplotWriter on a given level grid
      \ingroup Gnuplot
   */
  template<class G>
  class LevelGnuplotWriter : public GnuplotWriter<G, typename G::template Codim<0>::LevelIndexSet>
  {
  public:
    /** \brief Construct a Gnuplot writer for a certain level of a given grid */
    LevelGnuplotWriter (const G& grid, int level)
      : GnuplotWriter<G,typename G::template Codim<0>::LevelIndexSet>(grid,grid.levelIndexSet(level))
    {}
  };

}

#include "gnuplot/gnuplot.cc"

#endif // DUNE_IO_GNUPLOT_HH
