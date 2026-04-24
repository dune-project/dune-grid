// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <dune/grid/albertagrid.hh>
#include <dune/grid/io/file/dgfparser.hh>
#include "checkdgf.hh"

using namespace Dune;

template< int dim, int dimworld >
struct EnableSubIndexCheck< Dune::AlbertaGrid< dim, dimworld > >
{
  static const bool v = false;
};


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
