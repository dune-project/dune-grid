// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDTYPE_HH
#define DUNE_GRIDTYPE_HH

/**
 * @file
 * @author Andreas Dedner
 *
 * @brief This file can be included directly following config.h to
 *        test if a grid type was correctly selected.
 *
 **/

#ifndef HEADERCHECK

// NOGRID is used to specify that no default was set during configure
// If NOGRID and HAVE_GRIDTYPE are both not set then no grid was selected
// and an error is produced
#if defined NOGRID
  #if ! HAVE_GRIDTYPE
    #error "No grid type selected, use GRIDTYPE=..."
  #endif
#else
  #if ! HAVE_GRIDTYPE
    #error "No grid type selected, typo in GRIDTYPE=...?"
  #endif
#endif

#endif  // HEADERCHECK

#endif  // DUNE_GRIDTYPE_HH
