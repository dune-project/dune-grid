// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __GRAPEWRAPPER__
#define __GRAPEWRAPPER__

#include <string>

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

// make cast from const char * to char *
// otherwise wont work with gcc 4.2.x
inline void g_newerrorbox (const char * a, const char * b, int c, const char * d)
{
  g_errorbox(((char *)a),((char *)b),c,((char *)d));
}

// define new allert macro that uses g_newerrorbox
#define GRAPE_ALERT(condition,message,error_exit) \
  do {if(!(condition)) { \
        g_newerrorbox(message,__FILE__,__LINE__,# condition); \
        error_exit; \
      }} while(0)

#ifdef __cplusplus
}
#endif

#endif
