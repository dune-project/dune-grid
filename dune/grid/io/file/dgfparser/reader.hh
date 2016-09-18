// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_READER_HH
#define DUNE_DGF_READER_HH

#include <memory>
#include <string>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/gridreader.hh>
#include <dune/grid/io/file/dgfparser/dgfgridfactory.hh>

namespace Dune
{
  template <typename GridType>
  class DgfReader
      : private GridReader<GridType, DgfReader<GridType>>
  {
    typedef GridReader<GridType, DgfReader<GridType>> Base;

  public:
    typedef GridType Grid;

    static std::unique_ptr<Grid> read (const std::string& filename,
                                       MPIHelper::MPICommunicator comm = MPIHelper::getCommunicator())
    {
      DGFGridFactory<Grid> factory(filename, comm);

      return std::unique_ptr<Grid>{ factory.grid() };
    }

    static void read (Dune::GridFactory<Grid>& /*factory*/, const std::string& /*filename*/)
    {
      DUNE_THROW( DGFException, "DGF reader does not support to read into a GridFactory directly." );
    }
  };


} // end namespace Dune

#endif // END: INCLUDE-GUARD
