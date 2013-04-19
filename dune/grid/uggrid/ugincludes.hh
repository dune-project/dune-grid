// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGINCLUDES_HH
#define DUNE_UGINCLUDES_HH

/** \file
 * \brief All includes of UG headers in one single spot

   All UG includes have to be made from this file, and from
   this file only!  This is because undefAllMacros.pl takes
   all headers from this file and undefs the macros defined
   therein.
 */

#include <ug/gm.h>
#include <ug/std_domain.h>
#include <ug/initug.h>
#include <ug/commands.h>
#include <ug/formats.h>
#include <ug/elements.h>
#include <ug/shapes.h>
#include <ug/algebra.h>
#include <ug/refine.h>
#include <ug/ugm.h>
#include <ug/rm.h>
#if defined ModelP
#include <ug/parallel.h>
#endif

#endif
