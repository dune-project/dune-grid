// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_UTILITY_GRIDINFO_GMSH_MAIN_HH
#define DUNE_GRID_UTILITY_GRIDINFO_GMSH_MAIN_HH

#include <exception>
#include <iostream>
#include <ostream>
#include <string>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/utility/gridinfo.hh>

/** \file
 *  \author Jö Fahlke <jorrit@jorrit.de>
 *  \date 2011
 *
 * \brief Generic main() function for printing information about a mesh read
 *        from a .msh-file.
 *
 * This header contains a generic main() function.  To use it for your grid,
 * write a .cc file like this:
 * \code
 *#ifdef HAVE_CONFIG_H
 *#include "config.h"
 *#endif

 *#include <string>

 *#include <dune/grid/mygrid.hh>

   const std::string programName = "dune-gridinfo-gmsh-mygrid-3d";
   typedef Dune::MyGrid<3> Grid;

 *#include <dune/grid/utility/gridinfo-gmsh-main.hh>
 * \endcode
 * Write an automake target for your program as usual.  No special libraries
 * are needed for you program, beyond what is needed for the grid in question.
 */

#ifdef HEADERCHECK
// define so headercheck will run
const std::string programName = "headercheck";
#endif // HEADERCHECK

#ifndef DOXYGEN
namespace {
  // anonymous namespace so we don't freakishly conflict with another usage()
  // function that may be linked in from another compilation unit.
  void usage(std::ostream &stream) {
    stream << "USAGE:\n"
           << "  " << programName << " GRIDFILE\n"
           << "\n"
           << "PARAMTERS:\n"
           << "  GRIDFILE Name of the .msh file to read the grid from.\n"
           << std::flush;
  }
}

#ifndef HEADERCHECK
int main(int argc, char **argv) {
  try {
    const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // check that we are not run through mpirun
    if(mpiHelper.size() > 1) {
      if(mpiHelper.rank() == 0)
        std::cerr << programName << ": Sorry, this program works only in "
                  << "serial." << std::endl;
      return 1;
    }

    // check command line arguments
    if(argc < 2) {
      std::cerr << "Need name of a .msh file to read.\n"
                << std::endl;
      usage(std::cerr);
      return 1;
    }
    if(argc > 2) {
      std::cerr << "Too many arguments.\n"
                << std::endl;
      usage(std::cerr);
      return 1;
    }
    std::string gridFileName = argv[1];
    if(gridFileName == "-h" || gridFileName == "-?" ||
       gridFileName == "--help" || gridFileName == "-help")
    {
      usage(std::cout);
      return 0;
    }

    // read grid
    typedef Dune::GmshReader<Grid> Reader;
    Dune::shared_ptr<Grid> gridp(Reader::read(gridFileName));

    // collect information
    Dune::GridViewInfo<Grid::ctype> gridViewInfo;
    Dune::fillGridViewInfoSerial(gridp->leafView(), gridViewInfo);

    // print it
    std::cout << gridViewInfo << std::flush;
  }
  catch(const std::exception &e) {
    std::cerr << "Caught exception of type " << Dune::className(e)
              << std::endl
              << "e.what(): " << e.what() << std::endl;
    throw;
  }
  catch(const Dune::Exception &e) {
    std::cerr << "Caught exception of type " << Dune::className(e)
              << std::endl
              << "Exception message: " << e << std::endl;
    throw;
  }
  catch(const std::string &s) {
    std::cerr << "Caught exception of type " << Dune::className(s)
              << std::endl
              << "Exception message: " << s << std::endl;
    throw;
  }
  catch(...) {
    std::cerr << "Caught exception of unknown type" << std::endl;
    throw;
  }
}
#endif // !HEADERCHECK
#endif // !DOXYGEN

#endif // DUNE_GRID_UTILITY_GRIDINFO_GMSH_MAIN_HH
