// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
    \brief Test whether the DGFAlu parser can be used together with UGGrid

    By the time of writing, this did not work, because dgfalu.hh (indirectly)
    includes X.h, which defines a lot of preprocessor macros.  These (in
    particular the line #define Convex 2, screw up the uggrid headers.
 */

#include "config.h"

#include <dune/grid/io/file/dgfparser/dgfalu.hh>

#include <dune/grid/uggrid.hh>

int main()
{}
