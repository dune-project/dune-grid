// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAHEADER_HH
#define DUNE_ALBERTAHEADER_HH

#if HAVE_ALBERTA

#include <dune/grid/albertagrid/griddim.hh>

// if we have ALBERTA C++ lib define namespace for ALBERTA
#ifdef __ALBERTApp__
#define ALBERTA Alberta::
#else
#define ALBERTA ::
#endif

// the keyword ALBERTA stands for ALBERTA routines
#ifndef __ALBERTApp__
extern "C"
{
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

#ifndef DIM
#error "DIM or DIM_OF_WORLD not defined!"
#endif

#include <alberta.h>

// Macro nil may be defined by alberta_util.h. If so, undefine it.
#ifdef nil
#undef nil
#endif

// for version 1.2 thing are different
#if DUNE_ALBERTA_VERSION < 0x200

// face is not defined but should be the value of edge
//#ifndef FACE
//#define FACE EDGE
//#endif

static inline void meshTraverse(MESH *mesh,
                                int level, FLAGS fill_flag,
                                void (*el_fct)(const EL_INFO *))
{
  mesh_traverse(mesh, level, fill_flag, el_fct);
}

#define GET_EL_FROM_LIST(rc_list_el) (rc_list_el).el

//////////////////////////////////////////////////
#else // version 2.0
//////////////////////////////////////////////////

static inline void wrapped_el_fct(const EL_INFO* elinfo, void * data)
{
  ((void (*)(const EL_INFO *))data)(elinfo);
}

static inline void meshTraverse(MESH *mesh,
                                int level, FLAGS fill_flag,
                                void (*el_fct)(const EL_INFO *))
{
  mesh_traverse(mesh, level, fill_flag, wrapped_el_fct, (void*)el_fct );
}

#define GET_EL_FROM_LIST(rc_list_el) (rc_list_el).el_info.el

#endif // end Version 2.0
//////////////////////////////////////////////////////////////////////

#ifndef _ALBERTA_H_
#error "Couldn't find alberta.h for include! "
#endif

#ifndef __ALBERTApp__
} // end extern "C"
#endif

#endif // HAVE_ALBERTA

#endif
