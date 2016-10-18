#include <config.h>

#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/dgfparser.hh>
#include "../dgfoned.hh"
#include "test-dgf.cc"

using namespace Dune;

int main(int argc, char ** argv)
try {
  runDGFTest<OneDGrid>(argc,argv);
  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch (...)
{
  std::cerr << "Generic exception!" << std::endl;
  return 1;
}
