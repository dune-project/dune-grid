// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAHEADER_HH
#define DUNE_ALBERTAHEADER_HH

#if HAVE_ALBERTA

#if not (ALBERTA_DIM > 0)
  #if HEADERCHECK
    #undef ALBERTA_DIM
    #define ALBERTA_DIM 2
  #else
    #error ALBERTA_DIM should be 1, 2, or 3
  #endif
#endif

// Set ALBERTA's DIM_OF_WORLD preprocessor variable
#ifndef ALBERTA_DIM
#error "ALBERTA_DIM needed to use AlbertaGrid."
#endif
#define DIM_OF_WORLD ALBERTA_DIM

// if we have ALBERTA C++ lib define namespace for ALBERTA
#ifdef __ALBERTApp__
#define ALBERTA Alberta::
#else
#define ALBERTA ::
#endif

#ifndef ALBERTA_DEBUG
#define ALBERTA_DEBUG 0
#endif

// MAX, MIN, and ABS are defined macros of ALBERTA
// if they are not defined elsewhere, they are undefined here
#ifndef MAX
#define _MAX_NOT_DEFINED_
#endif

#ifndef MIN
#define _MIN_NOT_DEFINED_
#endif

#ifndef ABS
#define _ABS_NOT_DEFINED_
#endif

#ifndef DIM_OF_WORLD
#error "DIM_OF_WORLD not defined."
#endif

#ifdef HAVE_CONFIG_H
#define ALBERTASAVE_HAVE_CONFIG_H HAVE_CONFIG_H
#undef HAVE_CONFIG_H
#endif

#include <alberta/alberta.h>

#ifdef ALBERTASAVE_HAVE_CONFIG_H
#define HAVE_CONFIG_H ALBERTASAVE_HAVE_CONFIG_H
#undef ALBERTASAVE_HAVE_CONFIG_H
#endif

#ifndef _ALBERTA_H_
#error "Unable to include alberta.h."
#endif

// Macro nil may be defined by alberta_util.h. If so, undefine it.
#ifdef nil
#undef nil
#endif

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTAHEADER_HH
