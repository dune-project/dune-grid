// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFALLIMPLEMENTATIONS_HH
#define DUNE_DGFALLIMPLEMENTATIONS_HH
/* include the implementations */
#if HAVE_ALBERTA
#include "dgfalberta.hh"
#endif
#if HAVE_ALUGRID
#include "dgfalu.hh"
#endif
#if HAVE_UG
#include "dgfug.hh"
#endif
#include "dgfoned.hh"
#include "dgfyasp.hh"
#include "dgfs.hh"
#endif
