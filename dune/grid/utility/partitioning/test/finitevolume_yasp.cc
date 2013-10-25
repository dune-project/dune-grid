// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"               // know what grids are present

#include <string>

#include <dune/grid/yaspgrid.hh>

typedef Dune::YaspGrid<2> Grid;

static const std::string gridDescription = "Yasp<2>";
static const std::string gridPrefix = "yasp";

// number of global refines to double the effort
static const int dblref = 1;

#include "finitevolume_main.hh"
