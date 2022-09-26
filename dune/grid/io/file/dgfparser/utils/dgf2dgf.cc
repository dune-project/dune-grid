// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
 *  \author Matrin Nolte
 *  \brief a small program converting a DGF file into a DGF file
 *
 *  dgf2dgf is a small example program for the DGFWriter. It reads a DGF file
 *  into any grid, optionally refining the grid
 *  globally. The leaf grid is then written back as a DGF file.
 *
 *  The program's usage is as follows:
    \code
    ./dgf2dgf <dgffile> [refinement level]
    \endcode
 *
 *  While the program may seem completely useless, it has the following usages:
 *  - Convert an interval block into a simplex or cube grid (depending on the
 *    grid implementation used).
 *  - Resolve the simplex generator block into a vertex and a simplex block, so
 *    that it can be used without triangle or tetgen.
 *  - Construct a refined macro grid, which is very useful when setting up
 *    parallel computations with dune-ALUGrid.
 *  .
 *
 *  The source code of this program also demonstrates the easy use of the DGF
 *  parser and the DGFWriter.
 */

#include <config.h>

#include <iostream>

#include <dune/grid/io/file/dgfparser/dgfwriter.hh>

using namespace Dune;

int main ( int argc, char *argv[] )
try
{
  typedef GridSelector::GridType Grid;

  Dune::MPIHelper::instance( argc, argv );

  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ]
              << " <dgffile> [refinement level]" << std::endl;
    return 1;
  }

  const std::string dgfFileName( argv[ 1 ] );
  const int level = (argc < 3 ? 0 : atoi( argv[ 2 ] ));

  GridPtr< Grid > gridptr( dgfFileName );
  if( level > 0 )
    gridptr->globalRefine( level );

  typedef Grid::LeafGridView GridView;
  GridView gridView = gridptr->leafGridView();

  DGFWriter< GridView > dgfWriter( gridView );
  dgfWriter.write( std::cout );

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
