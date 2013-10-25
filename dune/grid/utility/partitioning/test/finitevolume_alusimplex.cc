// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"               // know what grids are present

#include <string>

#include <dune/grid/alugrid.hh>

typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming> Grid;

static const std::string gridDescription = "ALUSimplex<2,2>";
static const std::string gridPrefix = "alusimplex";

// number of global refines to double the effort
static const int dblref = 1;

#include "finitevolume_main.hh"
