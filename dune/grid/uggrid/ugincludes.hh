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

#include <initug.h>
#include <np/udm/formats.h>
#include <gm/elements.h>
#include <dom/std/std_domain.h>
#include <gm/algebra.h>
#include <gm/gm.h>
#include <gm/refine.h>
#include <gm/shapes.h>
#include <gm/rm.h>
#include <gm/ugm.h>
#include <ui/commands.h>
#if defined ModelP
#include <parallel/dddif/parallel.h>
#endif

#endif
