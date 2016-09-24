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
  /** \ingroup DuneGridFormatParser
   *
   *  \brief Read DGF mesh file
   *
   *  Read a .dgf file and construct a grid using the dgf-grid factory interface.
   *  Reader is implemented conforming the \ref GridReader interface. The exception is
   *  that the `read` method accepting a \ref GridFactory is not supported.
   */
  template <typename GridType>
  class DgfReader
      : public GridReader<GridType, DgfReader<GridType>>
  {
  protected:
    typedef GridType Grid;

    //! \brief Read .dgf file and return a unique_ptr to the created grid.
    //! Implementation of \ref GridReader::read.
    static std::unique_ptr<Grid> readImp (const std::string& filename,
                                          MPIHelper::MPICommunicator comm = MPIHelper::getCommunicator())
    {
      DGFGridFactory<Grid> factory(filename, comm);

      return std::unique_ptr<Grid>{ factory.grid() };
    }
  };


} // end namespace Dune

#endif // END: INCLUDE-GUARD
