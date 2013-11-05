// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <algorithm>
#include <bitset>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/array.hh>
#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include "to-vtk.hh"

const std::size_t dim = DIMENSION;

typedef Dune::YaspGrid<dim> Grid;

namespace {
  // anonymous namespace so we don't freakishly conflict with another help()
  // function that may be linked in from another compilation unit.
  void help(std::ostream &stream) {
    stream <<
"Convert a ALU-macrogridfile to VTK\n"
"\n"
"SYNOPSIS:\n"
"  gridnumbering-yasp-" << dim << "d [-s SIZE] [-e ELEMENTS] VTK_PREFIX\n"
"\n"
"PARAMETERS:\n"
"  -s SIZE Size of grid.  Either a space seperated list (one value per\n"
"    dimension) or a single value that is used for all dimensions.\n"
"  -e ELEMENTS Number of elements in each direction.  Either a space\n"
"    seperated list (one value per dimension) or a single value that is used\n"
"    for all dimensions.\n"
"  VTK_FILENAME Filename-prefix for VTK to write to.\n"
        << std::flush;
  }

  bool prefix_match(const std::string &prefix, const std::string &str)
  {
    return str.compare(0,prefix.size(), prefix) == 0;
  }

  bool equal_match(const std::string &a, const std::string &b)
  {
    return a == b;
  }

  void error_argument_required(const std::string &opt) {
    std::cerr << "Error: option " << opt << " requires argument\n";
    help(std::cerr);
    std::exit(1);
  }

  void error_unknown_option(const std::string &opt) {
    std::cerr << "Error: unknown option: " << opt << "\n";
    help(std::cerr);
    std::exit(1);
  }

  void error_parsing_optarg(const std::string &opt, const std::string &error) {
    std::cerr << "Error: option " << opt << ": " << error << "\n";
    help(std::cerr);
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

  template<class T>
  void parsevec(const std::string &arg, T &vec) {
    try {
      parse(arg, vec[0]);
      std::fill(vec.begin()+1, vec.end(), vec[0]);
    }
    catch(const std::runtime_error &) {
      std::istringstream s(arg);
      for(auto &val : vec)
        s >> val;
      bool good = !s.fail();
      if(good) {
        char dummy;
        s >> dummy;
        good = s.fail() && s.eof();
      }
      if(!good) {
        std::ostringstream s;
        s << "Can't parse \"" << arg << "\" as a " << Dune::className(vec);
        throw std::runtime_error(s.str());
      }
    }
  }

  bool getOptArgReq(char **&argv, const std::string &match,
                    std::string &opt, std::string &arg)
  {
    opt = match;
    if(*argv == opt) {
      ++argv;
      if(!*argv) error_argument_required(opt);
      arg = *argv;
      ++argv;
      return true;
    }
    if(prefix_match(opt, *argv)) {
      arg = *argv + opt.size();
      ++argv;
      return true;
    }
    return false;
  }

  bool getLongOptArgReq(char **&argv, const std::string &match,
                        std::string &opt, std::string &arg)
  {
    opt = match;
    if(*argv == opt) {
      ++argv;
      if(!*argv) error_argument_required(opt);
      arg = *argv;
      ++argv;
      return true;
    }
    if(prefix_match(opt+"=", *argv)) {
      arg = *argv + opt.size()+1;
      ++argv;
      return true;
    }
    return false;
  }

  Dune::FieldVector<Grid::ctype, dim> size(1);
  Dune::array<int, dim> elems;
  std::string vtkPrefix;

  void parseOptions(int argc, char **argv) {
    std::fill(elems.begin(), elems.end(), 10);

    std::vector<std::string> params;
    ++argv;
    while(*argv) {
      if(prefix_match("-", *argv)) {
        std::string opt, arg;
        if(equal_match("--", *argv)) {
          for(++argv; *argv; ++argv)
            params.push_back(*argv);
          break;
        }
        else if(prefix_match("-h", *argv) || prefix_match("-?", *argv) ||
                equal_match("--help", *argv))
        {
          help(std::cout);
          std::exit(0);
        }
        else if(getOptArgReq(argv, "-s", opt, arg) ||
                getLongOptArgReq(argv, "--size", opt, arg))
        {
          try { parsevec(*argv, size); }
          catch(const std::runtime_error &e)
          { error_parsing_optarg(opt, e.what()); }
        }
        else if(getOptArgReq(argv, "-e", opt, arg) ||
                getLongOptArgReq(argv, "--elems", opt, arg))
        {
          try { parsevec(*argv, elems); }
          catch(const std::runtime_error &e)
          { error_parsing_optarg(opt, e.what()); }
        }
        else
          error_unknown_option(*argv);
      }
      else {
        params.push_back(*argv);
        ++argv;
      }
    }
    // check command line arguments
    if(params.size() < 1) {
      std::cerr << "Need a prefix to name output vtk files.\n"
                << std::endl;
      help(std::cerr);
      std::exit(1);
    }
    if(params.size() > 1) {
      std::cerr << "Too many arguments.\n"
                << std::endl;
      help(std::cerr);
      std::exit(1);
    }
    vtkPrefix = params[0];
  }
}

int main(int argc, char **argv) {

  try {
    Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    parseOptions(argc, argv);

    Grid grid(mpiHelper.getCommunicator(), size, elems, std::bitset<dim>(), 0);

    numberingToVTK(grid.leafView(), vtkPrefix);
  }
  catch(const Dune::Exception &e) {
    std::cerr << "Caught Dune exception: " << e << std::endl;
    throw;
  }
  catch(const std::exception &e) {
    std::cerr << "Caught std::exception: " << e.what() << std::endl;
    throw;
  }
  catch(...) {
    std::cerr << "Caught unknown exception" << std::endl;
    throw;
  }
}
