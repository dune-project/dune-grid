// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_READER_HH
#define DUNE_GRID_READER_HH

#include <memory>
#include <string>
#include <utility>

#include <dune/grid/common/gridfactory.hh>

namespace Dune {

  template< class Grid, class GridReaderImp >
  class GridReader
  {
  protected:
    // type of underlying implementation, for internal use only
    typedef GridReaderImp Implementation;

  public:

    template <class... Args>
    static std::unique_ptr<Grid> read(const std::string &filename, Args&&... args)
    {
      GridFactory<Grid> factory;
      read(factory, filename, std::forward<Args>(args)...);

      return std::unique_ptr<Grid>{ factory.createGrid() };
    }

    template <class... Args>
    static void read(GridFactory<Grid> &factory, const std::string &filename, Args&&... args)
    {
      Implementation::read(factory, filename, std::forward<Args>(args)...);
    }

  };

} // end namespace Dune

#endif // END: INCLUDE-GUARD
