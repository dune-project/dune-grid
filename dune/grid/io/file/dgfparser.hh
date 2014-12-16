// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSER_HH
#define DUNE_DGFPARSER_HH
// include dgf parser
#include "dgfparser/dgfparser.hh"
/* include the implementations */
#include <dune/grid/albertagrid/dgfparser.hh>
#include "dgfparser/dgfalu.hh"
#include "dgfparser/dgfug.hh"
#include "dgfparser/dgfoned.hh"
#include "dgfparser/dgfyasp.hh"

#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
// this include is to be removed when SGrid is removed
#define DUNE_AVOID_SGRID_DEPRE_WARNING_BECAUSE_I_KNOW_WHAT_IM_DOING
#include "dgfparser/dgfs.hh"
#undef DUNE_AVOID_SGRID_DEPRE_WARNING_BECAUSE_I_KNOW_WHAT_IM_DOING
#endif

#include "dgfparser/dgfgeogrid.hh"
#include "dgfparser/dgfidentitygrid.hh"
#endif
