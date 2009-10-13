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

#include <gm.h>
#ifdef UG_LGMDOMAIN
#include <lgm_domain.h>
#else
#include <std_domain.h>
#endif
#include <initug.h>
#include <commands.h>
#include <formats.h>
#include <elements.h>
#include <shapes.h>
#include <algebra.h>
#include <refine.h>
#include <ugm.h>
#include <rm.h>
#if defined ModelP
#include "parallel.h"
#endif

#endif
