#include <config.h>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser.hh>
#include "../dgfyasp.hh"
#include "checkdgf.hh"

using namespace Dune;

int main(int argc, char ** argv)
try {
#ifdef TESTCOORDINATES
  using Grid=YaspGrid<2,EquidistantOffsetCoordinates<double,2>>;
#else
  using Grid=YaspGrid<3>;
#endif
  runDGFTest<Grid>(argc,argv);
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
