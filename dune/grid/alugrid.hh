// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_HH
#define DUNE_ALUGRID_HH

// only include this code, if HAVE_ALUGRID is true
#if HAVE_ALUGRID || DOXYGEN

#include <dune/grid/alugrid/common/declaration.hh>

#include <dune/grid/alugrid/3d/alugrid.hh>
#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>

// 2d version
#include <dune/grid/alugrid/2d/alugrid.hh>
#include <dune/grid/alugrid/2d/alu2dgridfactory.hh>

#include <dune/grid/alugrid/common/persistentcontainer.hh>

/** @file
    @author Robert Kloefkorn
    @brief Provides base classes for ALUGrid
 **/

//- include declaration of ALUGrid
#include <dune/grid/alugrid/common/declaration.hh>

#endif // #ifdef HAVE_ALUGRID || DOXYGEN

#endif
