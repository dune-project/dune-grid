// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAHEADER_HH
#define DUNE_ALBERTAHEADER_HH

#if HAVE_ALBERTA

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

// we dont use the el->index, its for debugging
#ifndef EL_INDEX
#define EL_INDEX 0
#else
#if EL_INDEX != 0
#warning "EL_INDEX != 0, but not used in interface implementation!\n"
#endif
#endif


#ifndef NEIGH_IN_EL
// neighbours were calculated on walkthrough
#define NEIGH_IN_EL 0
#else
#if NEIGH_IN_EL != 0
#error "NEIGH_IN_EL != 0 is not supported by this implementation!\n"
#endif
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

#include <alberta.h>

#ifndef _ALBERTA_H_
#error "Unable to include alberta.h."
#endif

// Macro nil may be defined by alberta_util.h. If so, undefine it.
#ifdef nil
#undef nil
#endif

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTAHEADER_HH
