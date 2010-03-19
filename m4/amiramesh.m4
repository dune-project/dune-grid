## -*- autoconf -*-
# $Id: dune_amira.m4 5156 2008-04-14 09:28:06Z christi $
# searches for amiramesh-headers and libs

# DUNE_PATH_AMIRAMESH()
# 
# shell variables:
#   with_amiramesh
#     no or yes
#   AMIRAMESHROOT
#   AMIRAMESH_LIB_PATH
#   AMIRAMESH_INCLUDE_PATH
#   AMIRAMESH_CPPFLAGS
#   AMIRAMESH_LDFLAGS
#   AMIRAMESH_LIBS
#   HAVE_AMIRAMESH
#     1 or 0 or undef
#
# substitutions:
#   AMIRAMESH_CPPFLAGS
#   AMIRAMESH_LDFLAGS
#   AMIRAMESH_LIBS
#
# defines:
#   HAVE_AMIRAMESH
#
# conditionals:
#   AMIRAMESH
AC_DEFUN([DUNE_PATH_AMIRAMESH],[
  AC_REQUIRE([AC_PROG_CXX])

  AC_ARG_WITH(amiramesh,
    AC_HELP_STRING([--with-amiramesh=PATH],[directory with AmiraMesh inside]))

# store values
ac_save_LDFLAGS="$LDFLAGS"
ac_save_CPPFLAGS="$CPPFLAGS"
ac_save_LIBS="$LIBS"

# initialize
HAVE_AMIRAMESH=0

## do nothing if --without-amiramesh is used
if test x$with_amiramesh != xno ; then

# is --with-amiramesh=bla used?
if test "x$with_amiramesh" != x ; then
	if ! test -d $with_amiramesh; then
        AC_MSG_WARN([Amiramesh directory $with_amiramesh does not exist])
	else
        # expand tilde / other stuff
		AMIRAMESHROOT=`cd $with_amiramesh && pwd`
	fi
fi
if test "x$AMIRAMESHROOT" = x; then
    # use some default value...
    AMIRAMESHROOT="/usr/local/amiramesh"
fi

AMIRAMESH_LIB_PATH="$AMIRAMESHROOT/lib"
AMIRAMESH_INCLUDE_PATH="$AMIRAMESHROOT/include"

CPPFLAGS="$CPPFLAGS -I$AMIRAMESH_INCLUDE_PATH"

AC_LANG_PUSH([C++])

# check for header
AC_CHECK_HEADER([amiramesh/AmiraMesh.h], 
   [AMIRAMESH_CPPFLAGS="-I$AMIRAMESH_INCLUDE_PATH -DHX_HAS_STDIOSTREAM"
	HAVE_AMIRAMESH="1"],
  AC_MSG_WARN([AmiraMesh.h not found in $AMIRAMESH_INCLUDE_PATH/amiramesh]),
  [#define HX_HAS_STDIOSTREAM])

CPPFLAGS="$ac_save_CPPFLAGS $AMIRAMESH_CPPFLAGS"

# if header is found...
if test x$HAVE_AMIRAMESH = x1 ; then
   LIBS="-L$AMIRAMESH_LIB_PATH -lamiramesh $LIBS"

   AC_LINK_IFELSE(AC_LANG_PROGRAM([#include "amiramesh/AmiraMesh.h"], [AmiraMesh* am = AmiraMesh::read("test");]),
	[AMIRAMESH_LIBS="-L$AMIRAMESH_LIB_PATH -lamiramesh"
         AMIRAMESH_LDFLAGS=""],
	[HAVE_AMIRAMESH="0"
	AC_MSG_WARN(libamiramesh not found!)])
fi

AC_LANG_POP([C++])

## end of amiramesh check (--without wasn't set)
fi

with_amiramesh="no"
# survived all tests?
if test x$HAVE_AMIRAMESH = x1 ; then
  AC_SUBST(AMIRAMESH_LIBS, $AMIRAMESH_LIBS)
  AC_SUBST(AMIRAMESH_LDFLAGS, $AMIRAMESH_LDFLAGS)
  AC_SUBST(AMIRAMESH_CPPFLAGS, $AMIRAMESH_CPPFLAGS)
  AC_DEFINE(HAVE_AMIRAMESH, 1, [Define to 1 if amiramesh-library is found])

  # add to global list
  DUNE_ADD_ALL_PKG([AmiraMesh], [\$(AMIRAMESH_CPPFLAGS)],
                   [\$(AMIRAMESH_LDFLAGS)], [\$(AMIRAMESH_LIBS)]) 

  # set variable for summary
  with_amiramesh="yes"
else
  AC_SUBST(AMIRAMESH_LIBS, "")
  AC_SUBST(AMIRAMESH_LDFLAGS, "")
  AC_SUBST(AMIRAMESH_CPPFLAGS, "")
fi
  
# also tell automake
AM_CONDITIONAL(AMIRAMESH, test x$HAVE_AMIRAMESH = x1)

# reset old values
LIBS="$ac_save_LIBS"
CPPFLAGS="$ac_save_CPPFLAGS"
LDFLAGS="$ac_save_LDFLAGS"

DUNE_ADD_SUMMARY_ENTRY([AmiraMesh],[$with_amiramesh])

])
