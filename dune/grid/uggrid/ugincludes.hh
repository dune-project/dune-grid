// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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


#include <dune/uggrid/gm/gm.h>
#include <dune/uggrid/domain/std_domain.h>
#include <dune/uggrid/initug.h>
#include <dune/uggrid/gm/elements.h>
#include <dune/uggrid/gm/algebra.h>
#include <dune/uggrid/gm/shapes.h>
#include <dune/uggrid/gm/refine.h>
#include <dune/uggrid/gm/ugm.h>
#include <dune/uggrid/gm/rm.h>
#if defined ModelP
#include <dune/uggrid/parallel/dddif/parallel.h>
#endif

#endif
