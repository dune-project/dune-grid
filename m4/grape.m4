## -*- autoconf -*-
# searches for albert-headers and libs

# grape.h und libgr.a/libgr.so are located in the same discretory 

# DUNE_PATH_GRAPE()
#
# configure shell/makefile variables:
#   GRAPE_CPPFLAGS
#   GRAPE_LDFLAGS
#   GRAPE_LIBS
#
# preprocessor defines:
#   HAVE_GRAPE ("ENABLE_GRAPE" or undefined)
#
# automake conditionals:
#   GRAPE
AC_DEFUN([DUNE_PATH_GRAPE],[
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PATH_XTRA])
  AC_REQUIRE([DUNE_PATH_OPENGL])
  AC_REQUIRE([AC_PROG_LD_GNU])

  AC_ARG_WITH(grape,
    AC_HELP_STRING([--with-grape=PATH],[directory with Grape inside]))

# store old values
ac_save_LDFLAGS="$LDFLAGS"
ac_save_CPPFLAGS="$CPPFLAGS"
ac_save_LIBS="$LIBS"

# don't even start testing if X wasn't found
if test "x$no_x" != xyes && test x$with_grape != xno ; then

  LIBS="$X_PRE_LIBS $X_LIBS $X_EXTRA_LIBS"

  # is --with-grape=bla used?
  if test x$with_grape != x ; then
    if test -d $with_grape; then
      # expand tilde / other stuff
      GRAPEROOT=`cd $with_grape && pwd`
    else
      AC_MSG_ERROR([directory $with_grape does not exist])
    fi      
  else
    # set some kind of default grape-path...
    GRAPEROOT="/usr/local/grape/"
  fi

  CPPFLAGS="$CPPFLAGS -I$GRAPEROOT"
  LDFLAGS="$LDFLAGS -L$GRAPEROOT"

  # check for header
  # we have to use CC for checking the header!!
  AC_LANG_PUSH([C])
  AC_CHECK_HEADER([grape.h],
    [GRAPE_CPPFLAGS="-I$GRAPEROOT"
     HAVE_GRAPE="1"])
  AC_LANG_POP

  # check for lib if header was found
  if test x$HAVE_GRAPE = x1 ; then
    # if GL was found, add it implicitly...
    #   This is not the best choice, but testing without GL first and
    #   then trying again fails due to caching...
    CPPFLAGS="$GRAPE_CPPFLAGS $GL_CFLAGS -DENABLE_GRAPE"
    LIBS="$LIBS $GL_LIBS -lXext"
    LDFLAGS="$LDFLAGS $GL_LDFLAGS"

    # if we use the gnu linker add the grape path 
    if test x$lt_cv_prog_gnu_ld = xyes ; then 
      GRAPE_LINKER_FLAGS="-Wl,--rpath -Wl,$GRAPEROOT"
    fi  

    AC_CHECK_LIB(gr, grape, 
      [GRAPE_LDFLAGS="$GL_LDFLAGS $GRAPE_LINKER_FLAGS"
       GRAPE_CPPFLAGS="$CPPFLAGS"
       GRAPE_LIBS="-L$GRAPEROOT -lgr $GL_LIBS -lXext"], 
      [HAVE_GRAPE="0"])
  fi

  # did it work?
  if test x$HAVE_GRAPE = x1 ; then
    AC_SUBST(GRAPE_LIBS, $GRAPE_LIBS)
    AC_SUBST(GRAPE_LDFLAGS, $GRAPE_LDFLAGS)
    AC_SUBST(GRAPE_CPPFLAGS, $GRAPE_CPPFLAGS)
    AC_DEFINE(HAVE_GRAPE, ENABLE_GRAPE,
          [This is only true if grape-library was found by configure 
           _and_ if the application uses the GRAPE_CPPFLAGS])

    # add to global list
    DUNE_ADD_ALL_PKG([GRAPE], [$GRAPE_CPPFLAGS], [$GRAPE_LDFLAGS], [$GRAPE_LIBS])
  fi
elif test "x$X_LIBS" = x ; then 
  AC_MSG_WARN([X libraries were not found and therefore not Grape check possible! See ./configure --help for X library options.])
fi


# report to summary
if test x$HAVE_GRAPE = x1 ; then
  with_grape="yes"
else
  with_grape="no"
fi

# also tell automake	
AM_CONDITIONAL(GRAPE, test x$HAVE_GRAPE = x1)

# reset old values
LIBS="$ac_save_LIBS"
CPPFLAGS="$ac_save_CPPFLAGS"
LDFLAGS="$ac_save_LDFLAGS"
  
DUNE_ADD_SUMMARY_ENTRY([Grape],[$with_grape])

])
