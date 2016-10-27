#include <config.h>

#include <dune/grid/albertagrid.hh>
#include <dune/grid/io/file/dgfparser.hh>
#include "checkdgf.hh"

using namespace Dune;

int main(int argc, char ** argv)
try {
  runDGFTest<AlbertaGrid<2>>(argc,argv);
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
