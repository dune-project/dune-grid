## -*- autoconf -*-
# searches for psurface-headers and lib

# DUNE_PATH_PSURFACE()
#
# DUNE_PATH_AMIRAMESH must be called previously if psurface was built with
# AmiraMesh support
#
# shell variables:
#   with_psurface
#     yes or no
#   PSURFACEROOT
#   PSURFACE_LIB_PATH
#   PSURFACE_INCLUDE_PATH
#   PSURFACE_CPPFLAGS
#   PSURFACE_LDFLAGS
#   PSURFACE_LIBS
#   HAVE_PSURFACE
#     0 or 1
#
# substitutions:
#   PSURFACE_LIBS
#   PSURFACE_LDFLAGS
#   PSURFACE_CPPFLAGS
#
# defines:
#   HAVE_PSURFACE
#
# conditionals:
#   PSURFACE
AC_DEFUN([DUNE_PATH_PSURFACE],[
  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([DUNE_PATH_AMIRAMESH])

  AC_ARG_WITH(psurface,
    AC_HELP_STRING([--with-psurface=PATH],[directory with the psurface library inside]))

# store values
ac_save_LDFLAGS="$LDFLAGS"
ac_save_CPPFLAGS="$CPPFLAGS"
ac_save_LIBS="$LIBS"

# initialize to sane value
HAVE_PSURFACE=0

## do nothing if --without-psurface is used
if test x$with_psurface != xno ; then

# is --with-psurface=bla used?
if test "x$with_psurface" != x ; then
    if test -d $with_psurface; then
      # expand tilde / other stuff
      PSURFACEROOT=`cd $with_psurface && pwd`
    else
      AC_MSG_ERROR([directory $with_psurface does not exist])
    fi      
fi
if test "x$PSURFACEROOT" = x; then  
    # use some default value...
    PSURFACEROOT="/usr/local/psurface"
fi

PSURFACE_LIB_PATH="$PSURFACEROOT/lib"
PSURFACE_INCLUDE_PATH="$PSURFACEROOT/include"

CPPFLAGS="$CPPFLAGS -I$PSURFACE_INCLUDE_PATH"

AC_LANG_PUSH([C++])

# check for header
AC_CHECK_HEADER([psurface/PSurface.h], 
   [PSURFACE_CPPFLAGS="-I$PSURFACE_INCLUDE_PATH"
	HAVE_PSURFACE="1"],
   [if test "x$with_psurface" != x ; then
    AC_MSG_WARN([psurface/PSurface.h not found in $PSURFACE_INCLUDE_PATH])
    fi
   ])

CPPFLAGS="$CPPFLAGS $PSURFACE_CPPFLAGS"

# if header is found...
if test x$HAVE_PSURFACE = x1 ; then
   AC_MSG_CHECKING([psurface library -lpsurface])

   # Why are the $AMIRAMESH_LDFLAGS $AMIRAMESH_LIBS here?
   # OS: This is a hack.  psurface can be compiled with and without AmiraMesh
   # support.  If it is compiled with AmiraMesh support, then it needs these flags
   # for the test to link (and since these flags must be set properly the AmiraMesh
   # test must have been called successfully before).  If psurface is compiled
   # without AmiraMesh support than the additional flags here do not matter.
   LIBS="-L$PSURFACE_LIB_PATH -lpsurface $AMIRAMESH_LIBS $LIBS"
   LDFLAGS="$LDFLAGS $AMIRAMESH_LDFLAGS"

   AC_LINK_IFELSE([AC_LANG_PROGRAM([#include "psurface/PSurface.h"], [[PSurface<2,double> foo;]])],
	[PSURFACE_LIBS="-L$PSURFACE_LIB_PATH -lpsurface"
         PSURFACE_LDFLAGS=""
         AC_MSG_RESULT(yes)],
	[HAVE_PSURFACE="0"
        AC_MSG_RESULT(no)
	AC_MSG_WARN([psurface header found, but libpsurface missing!])])
fi

AC_LANG_POP([C++])

## end of psurface check (--without wasn't set)
fi

with_psurface="no"
# survived all tests?
if test x$HAVE_PSURFACE = x1 ; then
  AC_SUBST(PSURFACE_LIBS, $PSURFACE_LIBS)
  AC_SUBST(PSURFACE_LDFLAGS, $PSURFACE_LDFLAGS)
  AC_SUBST(PSURFACE_CPPFLAGS, $PSURFACE_CPPFLAGS)
  AC_DEFINE(HAVE_PSURFACE, 1, [Define to 1 if psurface-library is found])

  # add to global list
  DUNE_ADD_ALL_PKG([psurface], [$PSURFACE_CPPFLAGS],
                   [$PSURFACE_LDFLAGS], [$PSURFACE_LIBS])

  # set variable for summary
  with_psurface="yes"
else
  AC_SUBST(PSURFACE_LIBS, "")
  AC_SUBST(PSURFACE_LDFLAGS, "")
  AC_SUBST(PSURFACE_CPPFLAGS, "")
fi
  
# also tell automake
AM_CONDITIONAL(PSURFACE, test x$HAVE_PSURFACE = x1)

# reset old values
LIBS="$ac_save_LIBS"
CPPFLAGS="$ac_save_CPPFLAGS"
LDFLAGS="$ac_save_LDFLAGS"

DUNE_ADD_SUMMARY_ENTRY([psurface],[$with_psurface])

])
