// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
 *  \author Matrin Nolte
 *  \brief a small program converting a gmsh file into a DGF file
 *
 *  gmsh2dgf is a small example program for the DGFWriter. It reads a gmsh file
 *  into any grid (selected by gridtype.hh) and writes it back as a DGF file.
 *
 *  The program's usage is as follows:
    \code
    ./gmsh2dgf <gmshfile>
    \endcode
 */

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/dgfparser/dgfwriter.hh>

using namespace Dune;

int main ( int argc, char *argv[] )
try
{
  Dune::MPIHelper::instance( argc, argv );
  typedef Dune::GridSelector::GridType Grid;

  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <gmshfile>" << std::endl;
    return 1;
  }

  const std::string gmshFileName( argv[ 1 ] );
  std::string dgfFileName( gmshFileName );
  dgfFileName.resize( dgfFileName.find_last_of( "." ) );
  dgfFileName += ".dgf";

  Grid *grid = GmshReader< Grid  >::read( gmshFileName );

  typedef Grid::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  DGFWriter< GridView > dgfWriter( gridView );
  dgfWriter.write( dgfFileName );

  const GridView::IndexSet &indexSet = gridView.indexSet();
  std::cerr << "Grid successfully written: "
            << indexSet.size( Grid::dimension ) << " vertices, "
            << indexSet.size( 0 ) << " elements."
            << std::endl;

  return 0;
}
catch( const Exception &exception )
{
  std::cerr << exception << std::endl;
  return 1;
}
