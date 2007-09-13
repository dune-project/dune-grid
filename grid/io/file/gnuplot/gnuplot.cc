// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** @file
    @author Christian Engwer
    @brief Implementation of gnuplot output for 1D Grids
 */

#include "../gnuplot.hh"

namespace Dune {

  /** \brief Write Gnuplot file to disk
      \param filename Name of the file to write to
   */
  template<class GridType, class IndexSet>
  void
  GnuplotWriter<GridType,IndexSet>::write(const std::string& filename) const
  {
    // open file
    std::ofstream file(filename.c_str());
    // write all column names
    file << "coord\t";
    for (size_t i=0; i<_names.size(); i++)
      file << _names[i] << "\t";
    file << "\n";

    int counter = 0;
    typedef typename IndexSet::template Codim<0>::template Partition<InteriorBorder_Partition>::Iterator CellIterator;
    CellIterator it = _is.template begin<0,InteriorBorder_Partition>();
    CellIterator end = _is.template end<0,InteriorBorder_Partition>();
    for (; it != end; ++it)
    {
      int i = _is.index(*it);
      // check that the elements are numbered consecutively
      assert (i == counter++);
      // calc positions
      assert(it->geometry().corners() == 2);
      float left = it->geometry()[0];
      float right = it->geometry()[1];
      assert(left < right);
      // write gnuplot rows for left & right vertex
      writeRow(file, left, _data[2*i]);
      writeRow(file, right, _data[2*i+1]);
    }
  }

  template<class GridType, class IndexSet>
  void
  GnuplotWriter<GridType,IndexSet>::writeRow(std::ostream & file, float position, const std::vector<float> & data) const
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
  template<class GridType, class IndexSet>
  template<class DataContainer>
  void
  GnuplotWriter<GridType,IndexSet>::addData(DataType t, const DataContainer& data, const std::string & name)
  {
    assert((t == cellData && (size_t) _is.size(0) == data.size())
           || (t == vertexData && (size_t) _is.size(0)+1 == data.size()) );
    _names.push_back(name);

    // copy data to new container
    // data is transformed to nonconforming vertex data
    int c = 0;
    int shift = (t==vertexData ? 1 : 0);
    for (int i=0; i<_is.size(0); i++)
    {
      _data[c++].push_back(data[i]);
      _data[c++].push_back(data[i+shift]);
    };
  }

}
