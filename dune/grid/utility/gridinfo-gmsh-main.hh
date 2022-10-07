// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_UTILITY_GRIDINFO_GMSH_MAIN_HH
#define DUNE_GRID_UTILITY_GRIDINFO_GMSH_MAIN_HH

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

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
   \code
   #ifdef HAVE_CONFIG_H
   #include "config.h"
   #endif

   #include <string>

   #include <dune/grid/mygrid.hh>

   const std::string programName = "dune-gridinfo-gmsh-mygrid-3d";
   typedef Dune::MyGrid<3> Grid;

   #include <dune/grid/utility/gridinfo-gmsh-main.hh>
   \endcode
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
           << "  " << programName << " [-R REFINES] GRIDFILE\n"
           << "\n"
           << "PARAMETERS:\n"
           << "  -R REFINES How many global refines to do after reading\n"
           << "    (default: 0)\n"
           << "  GRIDFILE Name of the .msh file to read the grid from.\n"
           << std::flush;
  }

  bool prefix_match(const std::string &prefix, const std::string &str)
  {
    return str.compare(0,prefix.size(), prefix) == 0;
  }

  void error_argument_required(const std::string &opt) {
    std::cerr << "Error: option " << opt << " requires argument\n";
    usage(std::cerr);
    std::exit(1);
  }

  void error_unknown_option(const std::string &opt) {
    std::cerr << "Error: unknown option: " << opt << "\n";
    usage(std::cerr);
    std::exit(1);
  }

  void error_parsing_optarg(const std::string &opt, const std::string &error) {
    std::cerr << "Error: option " << opt << ": " << error << "\n";
    usage(std::cerr);
    std::exit(1);
  }

  template<class T>
  void parse(const std::string &arg, T &val) {
    std::istringstream s(arg);
    s >> val;
    bool good = !s.fail();
    if(good) {
      char dummy;
      s >> dummy;
      good = s.fail() && s.eof();
    }
    if(!good) {
      std::ostringstream s;
      s << "Can't parse \"" << arg << "\" as a " << Dune::className(val);
      throw std::runtime_error(s.str());
    }
  }

  std::size_t refines = 0;
  std::string gridFileName = "";

  void parseOptions(int argc, char **argv) {
    std::vector<std::string> params;
    for(++argv; *argv; ++argv) {
      std::string arg = *argv;
      if(prefix_match("-", arg)) {
        std::string opt = arg;
        if(opt == "--") {
          for(++argv; *argv; ++argv)
            params.push_back(*argv);
          break;
        }
        else if(prefix_match("-h", opt) || prefix_match("-?", opt) ||
                opt == "--help")
        {
          usage(std::cout);
          std::exit(0);
        }
        else if(opt == "-R" || opt == "--global-refines") {
          ++argv;
          if(!*argv) error_argument_required(opt);
          try { parse(*argv, refines); }
          catch(const std::runtime_error &e)
          { error_parsing_optarg(opt, e.what()); }
        }
        else if(prefix_match("-R", opt)) {
          try { parse(*argv+std::strlen("-R"), refines); }
          catch(const std::runtime_error &e)
          { error_parsing_optarg(opt, e.what()); }
        }
        else if(prefix_match("--global-refines=", opt)) {
          try { parse(*argv+std::strlen("--global-refines="), refines); }
          catch(const std::runtime_error &e)
          { error_parsing_optarg(opt, e.what()); }
        }
        else
          error_unknown_option(opt);
      }
      else
        params.push_back(arg);
    }
    // check command line arguments
    if(params.size() < 1) {
      std::cerr << "Need name of a .msh file to read.\n"
                << std::endl;
      usage(std::cerr);
      std::exit(1);
    }
    if(params.size() > 1) {
      std::cerr << "Too many arguments.\n"
                << std::endl;
      usage(std::cerr);
      std::exit(1);
    }
    gridFileName = params[0];
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

    parseOptions(argc, argv);

    // read grid
    typedef Dune::GmshReader<Grid> Reader;
    std::shared_ptr<Grid> gridp(Reader::read(gridFileName));
    gridp->globalRefine(refines);

    // collect information
    Dune::GridViewInfo<Grid::ctype> gridViewInfo;
    Dune::fillGridViewInfoSerial(gridp->leafGridView(), gridViewInfo);

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
