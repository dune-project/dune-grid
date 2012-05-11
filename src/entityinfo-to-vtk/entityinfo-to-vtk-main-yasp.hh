// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_SRC_ENTITIYINFO_TO_VTK_ENTITIYINFO_TO_VTK_MAIN_YASP_HH
#define DUNE_GRID_SRC_ENTITIYINFO_TO_VTK_ENTITIYINFO_TO_VTK_MAIN_YASP_HH

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/yaspgrid.hh>

#include "entityinfo-to-vtk.hh"

/** \file
 *  \author Jö Fahlke <jorrit@jorrit.de>
 *  \date 2012
 *
 * \brief Generic main() function writing information about a YaspGrid's
 *       entities to VTK files,
 *
 * This header contains a generic main() function.  To use it for your grid,
 * write a .cc file like this:
 * \code
 *#ifdef HAVE_CONFIG_H
 *#include "config.h"
 *#endif

 *#include <string>

   const std::string programName = "dune-entityinfo-to-vtk-yasp-3d";
   static const int dim = 3;

 *#include "entityinfo-to-vtk-main-yasp.hh"
 * \endcode
 * Write an automake target for your program as usual.  No special libraries
 * are needed for you program.
 */

#ifdef HEADERCHECK
// define so headercheck will run
const std::string programName = "headercheck";
static const int dim = 3;
#endif // HEADERCHECK

#ifndef DOXYGEN
namespace {
  typedef Dune::YaspGrid<dim> Grid;
  typedef Grid::ctype ctype;

  // anonymous namespace so we don't freakishly conflict with another usage()
  // function that may be linked in from another compilation unit.
  void usage(std::ostream &stream) {
    stream <<
    "USAGE:\n"
    "  " << programName << " [-E ELEMENTS] [-O OVERLAP] [-R REFINES] [-S SIZE]\n"
    "        [-T OUTPUTTYPE] VTKNAME\n"
    "\n"
    "PARAMTERS:\n"
    "  -E ELEMENTS Whitespace seperated vector of number of elements for each\n"
    "    coordinate axis (default: 1 for each dimension)\n"
    "  -O OVERLAP Size of the overlap (default: 1)\n"
    "  -R REFINES How many global refines to do after creating the grid\n"
    "    (default: 0)\n"
    "  -S SIZE Whitespace seperated vector of domain sizes for each coordinate\n"
    "    axis (default: 1 for each dimension)\n"
    "  -T OUTPUTTYPE Mode to use when writing VTK files.  One of ascii, base64,\n"
    "    appendedraw (default), appendedbase64.\n"
    "  VTKNAME Name for the VTK output files.  This can include a directory\n"
    "    component.  Each VTK file will have a prefix inserted between directory\n"
    "    and basename and a suffix appended.\n"
    << std::flush;
  }

  bool prefix_match(const std::string &prefix, const std::string &str)
  {
    return str.compare(0,prefix.size(), prefix) == 0;
  }

  void error_argument_required(const std::string &opt, bool silent) {
    if(!silent) {
      std::cerr << "Error: option " << opt << " requires argument\n\n";
      usage(std::cerr);
    }
    std::exit(1);
  }

  void error_unknown_option(const std::string &opt, bool silent) {
    if(!silent) {
      std::cerr << "Error: unknown option: " << opt << "\n\n";
      usage(std::cerr);
    }
    std::exit(1);
  }

  void error_parsing_optarg(const std::string &opt, const std::string &error,
                            bool silent)
  {
    if(!silent) {
      std::cerr << "Error: option " << opt << ": " << error << "\n\n";
      usage(std::cerr);
    }
    std::exit(1);
  }

  bool opt_with_arg(const std::string &prefix, char **&argv, std::string &rest)
  {
    if(prefix_match(prefix, *argv)) {
      rest = *argv + prefix.size();
      ++argv;
      return true;
    }
    else
      return false;
  }

  bool opt_and_arg(const std::string &opt, char **&argv, std::string &optarg,
                   bool silent)
  {
    bool match = opt == *argv;
    if(match) {
      ++argv;
      if(!*argv)
        error_argument_required(opt, silent);
      optarg = *argv;
      ++argv;
    }
    return match;
  }

  bool optarg(const std::string &shortopt, const std::string &longopt,
              char **&argv, std::string &optarg, bool silent)
  {
    return opt_and_arg(shortopt, argv, optarg, silent) ||
           opt_with_arg(shortopt, argv, optarg) ||
           opt_and_arg(longopt, argv, optarg, silent) ||
           opt_with_arg(longopt+"=", argv, optarg);
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

  template<class T, int size>
  void parse(const std::string &arg, Dune::FieldVector<T, size> &val) {
    std::istringstream s(arg);
    for(int i = 0; i < size; ++i)
      s >> val[i];
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

  void parse(const std::string &arg, Dune::VTK::OutputType &val) {
    std::istringstream s(arg);
    std::string tmp;
    s >> tmp;
    bool good = !s.fail();
    if(good) {
      char dummy;
      s >> dummy;
      good = s.fail() && s.eof();
    }
    if(good) {
      if(tmp == "ascii")
        val = Dune::VTK::ascii;
      else if(tmp == "base64")
        val = Dune::VTK::base64;
      else if(tmp == "appendedraw")
        val = Dune::VTK::appendedraw;
      else if(tmp == "appendedbase64")
        val = Dune::VTK::appendedbase64;
      else
        good = false;
    }
    if(!good) {
      std::ostringstream s;
      s << "Can't parse \"" << arg << "\" as a " << Dune::className(val) << " "
        << "(valid values are: ascii, base64, appendedraw, appendedbase64)";
      throw std::runtime_error(s.str());
    }
  }

  Dune::FieldVector<int, dim> elements(1);
  int overlap = 1;
  std::size_t refines = 0;
  Dune::FieldVector<ctype, dim> size(1);
  Dune::VTK::OutputType outputtype = Dune::VTK::appendedraw;
  std::string vtkName = "";

  void parseOptions(char **argv, bool silent) {
    std::vector<std::string> params;
    std::string arg;
    ++argv;
    while(*argv)
      try {
        arg = *argv;
        if(arg == "-h" || arg == "-?" || arg == "--help") {
          if(!silent)
            usage(std::cout);
          std::exit(0);
        }
        else if(optarg("-E", "--elements", argv, arg, silent))
          parse(arg, elements);
        else if(optarg("-O", "--overlap", argv, arg, silent))
          parse(arg, overlap);
        else if(optarg("-R", "--global-refines", argv, arg, silent))
          parse(arg, refines);
        else if(optarg("-S", "--size", argv, arg, silent))
          parse(arg, size);
        else if(optarg("-T", "--output-type", argv, arg, silent))
          parse(arg, outputtype);
        else if(arg == "--")
          for(++argv; *argv; ++argv)
            params.push_back(*argv);
        else if(prefix_match("-", arg))
          error_unknown_option(arg, silent);
        else {
          params.push_back(arg);
          ++argv;
        }
      }
      catch(const std::runtime_error &e)
      { error_parsing_optarg(arg, e.what(), silent); }
    // check command line arguments
    if(params.size() < 1) {
      if(!silent) {
        std::cerr << "Need name for VTK files.\n\n";
        usage(std::cerr);
      }
      std::exit(1);
    }
    if(params.size() > 1) {
      if(!silent) {
        std::cerr << "Too many arguments.\n\n";
        usage(std::cerr);
      }
      std::exit(1);
    }
    vtkName = params[0];
  }
}

#ifndef HEADERCHECK
int main(int argc, char **argv) {
  try {
    const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    parseOptions(argv, mpiHelper.rank() != 0);

    Grid grid(mpiHelper.getCommunicator(), size, elements,
              Dune::FieldVector<bool, dim>(false), overlap);
    grid.globalRefine(refines);

    entityinfoToVTK(grid.leafView(), vtkName, outputtype);
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

#endif // DUNE_GRID_SRC_ENTITIYINFO_TO_VTK_ENTITIYINFO_TO_VTK_MAIN_YASP_HH
