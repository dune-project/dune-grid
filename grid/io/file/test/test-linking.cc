// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id: test-sgrid.cc 5487 2009-09-25 06:39:23Z mnolte $

#include <config.h>

#include <iostream>

#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/virtualrefinement.hh>

/*
   #include "gridcheck.cc"
   #include "checkgeometryinfather.cc"
   #include "checkintersectionit.cc"
   #include "checkpartition.cc"
 */

template <class Grid>
void test(Grid& grid)
{}
