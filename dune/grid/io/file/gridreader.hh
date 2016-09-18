// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_READER_HH
#define DUNE_GRID_READER_HH

#include <memory>
#include <string>
#include <utility>

#include <dune/grid/common/gridfactory.hh>

namespace Dune {

  /** \brief Grid reader abstract base class
   *
   *  Interface class for a file-reader for grids. The reader interface defines two
   *  methods to read grids from file: 1) Given a filename it returns a unique_ptr to the
   *  created grid. 2) Given a \ref GridFactory and filename it fills the factory and the
   *  user can create the grid from the factory by calling \ref createGrid() manually.
   *
   *  The interface provides a default implementation for the `read(filename)` method, by
   *  creating a \ref GridFactory and calling the second read method on this factory instance.
   *
   *  The factory `read` method redirects to the concrete implementation provided as template
   *  parameter `GridReaderImp`.
   *
   *  Both `read` methods forward additional arguments to the concrete implementation of the
   *  derived class.
   *
   *  Example of usage:
   *  ```
   *  typedef Dune::UGGrid<3> Grid;
   *  typedef Dune::AlbertaReader<Grid> Reader;
   *  auto grid = Reader::read("docs/grids/amc/grid-3-3.amc");
   *  ```
   */
  template< class Grid, class GridReaderImp >
  class GridReader
  {
  protected:
    // type of underlying implementation, for internal use only
    typedef GridReaderImp Implementation;

  public:

    //! Read the grid from a file with filename and return a unique_ptr to the created grid
    template <class... Args>
    static std::unique_ptr<Grid> read(const std::string &filename, Args&&... args)
    {
      GridFactory<Grid> factory;
      read(factory, filename, std::forward<Args>(args)...);

      return std::unique_ptr<Grid>{ factory.createGrid() };
    }

    //! Read the grid from a file with filename into a grid-factory. Must be implemented by derived class.
    template <class... Args>
    static void read(GridFactory<Grid> &factory, const std::string &filename, Args&&... args)
    {
      Implementation::read(factory, filename, std::forward<Args>(args)...);
    }

  };

} // end namespace Dune

#endif // END: INCLUDE-GUARD
