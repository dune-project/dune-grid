// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __GRAPEWRAPPER__
#define __GRAPEWRAPPER__

#ifdef __cplusplus
extern "C" {
#endif

#define class GrapeClass
#define private GrapePrivate
#define friend GrapeFriend
#define explicit GrapeExplicit

#define G_CPP
#include <grape.h>
#undef G_CPP

#undef class
#undef private
#undef friend
#undef explicit

// make cast from const char * to char *
// otherwise wont work with gcc 4.2.x
#define GRAPE_CALL(obj,meth) GRAPE(obj,((char *)meth))

#ifdef __cplusplus
}
#endif

#endif
