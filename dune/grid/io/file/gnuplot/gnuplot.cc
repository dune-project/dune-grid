// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** @file
    @author Christian Engwer
    @brief Implementation of gnuplot output for 1D and 2D grids
 */

#include "../gnuplot.hh"

namespace Dune {

  /** \brief Write Gnuplot file to disk
      \param filename Name of the file to write to
   */
  template<class GridView>
  void
  GnuplotWriter<GridView>::write(const std::string& filename) const
  {
    // open file
    std::ofstream file(filename.c_str());
    // write all column names
    file << "# coord\t";
    for (size_t i=0; i<_names.size(); i++)
      file << _names[i] << "\t";
    file << "\n";

    if (dimworld==1) {
#if !NDEBUG
      int counter = 0;
#endif
      typedef typename GridView::template Codim<0>::Iterator CellIterator;
      CellIterator it = _gv.template begin<0>();
      CellIterator end = _gv.template end<0>();
      for (; it != end; ++it)
      {
        int i = _is.index(*it);
        // check that the elements are numbered consecutively
        assert (i == counter++);
        // calc positions
        assert(it->geometry().corners() == 2);
        const FieldVector<ctype,dimworld>& left = it->geometry().corner(0);
        const FieldVector<ctype,dimworld>& right = it->geometry().corner(1);
        assert(left[0] < right[0]);
        // write gnuplot rows for left & right vertex
        writeRow(file, left, _data[2*i]);
        writeRow(file, right, _data[2*i+1]);
      }

    } else {

      typedef typename GridView::template Codim<dimworld>::Iterator VertexIterator;
      VertexIterator it  = _gv.template begin<dimworld>();
      VertexIterator end = _gv.template end<dimworld>();
      for (; it != end; ++it) {

        // write gnuplot rows for vertex
        writeRow(file, it->geometry().corner(0), _data[_is.index(*it)]);

      }

    }

  }

  template<class GridView>
  void
  GnuplotWriter<GridView>::writeRow(std::ostream & file,
                                    const FieldVector<ctype,dimworld>& position,
                                    const std::vector<float> & data) const
  {
    assert (data.size() == _names.size());
    // write position
    file << position << "\t";
    // write all data columns
    for (size_t j=0; j<data.size(); j++)
      file << data[j] << "\t";
    file << "\n";
  }

  /** \brief Add data (internal)
      \param t type of data (cellData/vertexData)
      \param data An ISTL compliant vector type
      \param name associated with the data
   */
  template<class GridView>
  template<class DataContainer>
  void
  GnuplotWriter<GridView>::addData(DataType t, const DataContainer& data, const std::string & name)
  {
    assert((t == cellData && _is.size(0) == data.size())
           || (t == vertexData && _is.size(GridView::dimension) == data.size()) );
    _names.push_back(name);

    // copy data to new container

    if (dimworld==1) {

      // data is transformed to nonconforming vertex data
      int c = 0;
      int shift = (t==vertexData ? 1 : 0);
      for (size_t i=0; i<_is.size(0); i++)
      {
        _data[c++].push_back(data[i]);
        _data[c++].push_back(data[i+shift]);
      };

    } else {

      // 2d: only vertex data is allowed
      for (size_t i=0; i<_is.size(dimworld); i++)
        _data[i].push_back(data[i]);

    }

  }

}
