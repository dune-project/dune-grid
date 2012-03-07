// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_UTILITY_STRUCTURED_TO_ALU_MAIN_HH
#define DUNE_GRID_UTILITY_STRUCTURED_TO_ALU_MAIN_HH

/** @file
 *
 * This file is meant to be included directly by a .cc file.  The including
 * file must setup the following names in the global namespace prior to the
 * include directive:
 *
 * - \c Grid        Type of the grid to be used.
 * - \c programName A std::string initialized to the name of the program for
 *                  help messages etc.
 * - \c useCube     A boolean constant suitable for use as a template
 *                  argument.  Whether to use cubes or simplices.
 *
 * The including file is also responsible for including any headers necessary
 * to setup the names above -- this includes "config.h".
 */

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>

#include <dune/grid/utility/structuredgridfactory.hh>

void help(std::ostream &s) {
  s <<
  "Generate a structured with the StructuredGridFactory, distribute it and\n"
  "write it out as a partitioned ALU-macrogridfile\n"
  "\n"
  "SYNOPSIS:\n"
  "  mpirun -np SIZE " << programName << " [--] ALU_PREFIX LOWER UPPER ELEMS\n"
  "  " << programName << " {-h|--help}\n"
  "\n"
  "PARAMETERS:\n"
  "  -h, --help Print this help.\n"
  "  SIZE Number of ranks to distribute to.\n"
  "  ALU_PREFIX Filename-prefix for ALUGrid to write the macrogrid to.  Each\n"
  "    process appends its rank to this filename to form something like\n"
  "    ALUPREFIX.RANK.\n"
  "  LOWER, UPPER Lower and upper world coordinates of the box to create the\n"
  "    structured mesh in.  Seperate individual components with spaces, but\n"
  "    make sure to escape them from the shell.\n"
  "  ELEMS Space seperated list of size dimgrid.  Each list entry gives the\n"
  "    number elements to split the mesh into for the corresponding coordinate\n"
  "    direction.\n"
  << std::flush;
}

template<class OutputIterator>
bool parseRangeHelper(OutputIterator it, const OutputIterator &end,
                      const std::string &arg)
{
  std::istringstream s(arg);
  for(; it != end; ++it)
    s >> *it;
  if(s.fail() || s.bad())
    return false;
  char dummy;
  s >> dummy;
  if(!s.fail())
    return false;
  return true;
}

template<class Range>
void parseRange(int rank, Range &range, const std::string &arg,
                const std::string &argname)
{
  if(!parseRangeHelper(range.begin(), range.end(), arg)) {
    if(rank == 0) {
      std::cerr << programName << ": Can't parse commandline argument \""
                << arg << "\" for parameter " << argname << ".\n\n";
      help(std::cerr);
    }
    std::exit(1);
  }
}

template<class Coord, class Elems>
void parseOptions(int rank, int argc, char **argv, std::string &aluName,
                  Coord &lower, Coord &upper, Elems &elems)
{
  std::vector<const char *> params;
  bool noMoreSwitches = false;
  for(char const* const* argp = argv+1; *argp; ++argp) {
    if(!noMoreSwitches && **argp == '-') {
      std::string opt = *argp;
      if(opt == "--")
        noMoreSwitches = true;
      else if(opt == "-h" || opt == "--help") {
        if(rank == 0)
          help(std::cout);
        std::exit(0);
      }
      else {
        if(rank == 0) {
          std::cerr << programName << ": Unknown option: " << opt << "\n\n";
          help(std::cerr);
        }
        std::exit(1);
      }
    }
    else
      params.push_back(*argp);
  }
  if(params.size() < 4) {
    if(rank == 0) {
      std::cerr << programName << ": Too few commandline arguments\n\n";
      help(std::cerr);
    }
    std::exit(1);
  }
  if(params.size() > 4) {
    if(rank == 0) {
      std::cerr << programName << ": Too many commandline arguments\n\n";
      help(std::cerr);
    }
    std::exit(1);
  }

  aluName = params[0];
  parseRange(rank, lower, params[1], "LOWER");
  parseRange(rank, upper, params[2], "UPPER");
  parseRange(rank, elems, params[3], "ELEMS");
}

template<bool = useCube>
struct GridCreator;

template<>
struct GridCreator<false> {
  typedef Grid::ctype ctype;
  static const unsigned dim = Grid::dimension;
  static const unsigned dimw = Grid::dimensionworld;

  static Dune::shared_ptr<Grid>
  create(const Dune::FieldVector<ctype, dimw> &lower,
         const Dune::FieldVector<ctype, dimw> &upper,
         const Dune::array<unsigned, dim> &elems)
  {
    return Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower, upper,
                                                                elems);
  }
};

template<>
struct GridCreator<true> {
  typedef Grid::ctype ctype;
  static const unsigned dim = Grid::dimension;
  static const unsigned dimw = Grid::dimensionworld;

  static Dune::shared_ptr<Grid>
  create(const Dune::FieldVector<ctype, dimw> &lower,
         const Dune::FieldVector<ctype, dimw> &upper,
         const Dune::array<unsigned, dim> &elems)
  {
    return Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper,
                                                             elems);
  }
};

int main(int argc, char **argv) {

  try {

    typedef Grid::ctype ctype;
    static const unsigned dim = Grid::dimension;
    static const unsigned dimw = Grid::dimensionworld;

    Dune::Timer timer;

    const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    //////////////////////////////////////////////////////////////////////
    //
    //  parse commandline
    //

    std::string aluName;
    Dune::FieldVector<ctype, dimw> lower, upper;
    Dune::array<unsigned, dim> elems;
    parseOptions(mpiHelper.rank(), argc, argv, aluName, lower, upper, elems);

    //////////////////////////////////////////////////////////////////////
    //
    //  create Grid
    //

    Dune::shared_ptr<Grid> gridp(GridCreator<>::create(lower, upper, elems));

    //////////////////////////////////////////////////////////////////////
    //
    //  actual loadbalancing
    //

    if(mpiHelper.rank() == 0)
      std::cout << "[" << timer.elapsed() << "] Load-balancing grid"
                << std::endl;
    gridp->loadBalance();

    //////////////////////////////////////////////////////////////////////
    //
    //  write grid
    //

    {
      if(mpiHelper.rank() == 0)
        std::cout << "[" << timer.elapsed() << "] Writing grid files "
                  << aluName << ".*" << std::endl;

      std::string path;
      std::string basename;

      std::size_t pos = aluName.rfind('/');
      switch(pos) {
      case std::string::npos :
        path = ".";
        basename = aluName;
        break;
      case 0 :
        path = "/";
        basename = aluName.substr(1);
        break;
      default :
        path = aluName.substr(0, pos);
        basename = aluName.substr(pos+1);
        break;
      }
      gridp->writeMacroGrid(path, basename);
    }

    if(mpiHelper.rank() == 0)
      std::cout << "[" << timer.elapsed() << "] Done." << std::endl;
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

#endif // DUNE_GRID_UTILITY_STRUCTURED_TO_ALU_MAIN_HH
