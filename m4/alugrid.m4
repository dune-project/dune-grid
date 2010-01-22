dnl -*- mode: autoconf; tab-width: 4; indent-tabs-mode: nil; -*-
# searches for alugrid-headers and libs

# DUNE_PATH_ALUGRID()
#
# shell variables:
#   with_alugrid
#     no or yes
#   ALUGRIDROOT
#   ALUGRID_VERSIONNO
#   ALUGRID_LIB_PATH
#   ALUGRID_INCLUDE_PATH
#   ALUGRID_CPPFLAGS
#   ALUGRID_LDFLAGS
#   ALUGRID_LIBS
#   HAVE_ALUGRID
#     undef or 1 or 0
#
# substitutions:
#   ALUGRID_CPPFLAGS
#   ALUGRID_LDFLAGS
#   ALUGRID_LIBS
#
# defines:
#   HAVE_ALUGRID
#     ENABLE_ALUGRID or undefined
#   ALUGRID_PARALLEL_H
#   ALUGRID_SERIAL_H
#
# conditionals:
#   ALUGRID
AC_DEFUN([DUNE_PATH_ALUGRID],[
  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([DUNE_MPI])

  AC_ARG_WITH(alugrid,
    AC_HELP_STRING([--with-alugrid=PATH],[directory where ALUGrid is installed]))
# do not use alugrid debug lib 

# store old values
ac_save_LDFLAGS="$LDFLAGS"
ac_save_CPPFLAGS="$CPPFLAGS"
ac_save_LIBS="$LIBS"

# initilize to sane value
HAVE_ALUGRID=0

## do nothing if no --with-alugrid was supplied
if test x$with_alugrid != x && test x$with_alugrid != xno ; then

  if test x$with_alugrid = xyes ; then
    # use some default value...
    ALUGRIDROOT="/usr/local/alugrid"
  else
    ALUGRIDROOT="$with_alugrid"
  fi

  # expand tilde / other stuff
  ALUGRIDROOT=`cd "$ALUGRIDROOT" && pwd`

  if ! test -d "$ALUGRIDROOT"; then
    AC_MSG_ERROR([ALUGrid directory $ALUGRIDROOT does not exist or is inaccessible])
  fi

    ALUGRID_VERSIONCHECK=$ALUGRIDROOT/bin/alugridversion
    ## check version number 
    NEEDEDALUGRID_VERSION=1.1

    AC_MSG_CHECKING([ALUGrid version >= $NEEDEDALUGRID_VERSION])
    if test -f $ALUGRID_VERSIONCHECK; then 
      ALUGRID_VERSION=`$ALUGRID_VERSIONCHECK -c $NEEDEDALUGRID_VERSION`
      if test "x$ALUGRID_VERSION" != "x-1"; then 
        ALUGRID_VERSIONNO=`$ALUGRID_VERSIONCHECK -v`
        AC_MSG_RESULT([yes (ALUGrid-$ALUGRID_VERSIONNO)])
      else 
        AC_MSG_ERROR([ALUGrid version is too old!])
      fi 
      
    else 
      AC_MSG_ERROR([Couldn't find ALUGrid version checker! ALUGrid version too old or ALUGrid not installed in $ALUGRIDROOT!])
    fi   

  ALUGRID_LIB_PATH="$ALUGRIDROOT/lib"
  ALUGRID_INCLUDE_PATH="$ALUGRIDROOT/include"

  AC_LANG_PUSH([C++])

  # set variables so that tests can use them
  ALU3D_INC_FLAG="-I$ALUGRID_INCLUDE_PATH -I$ALUGRID_INCLUDE_PATH/serial -I$ALUGRID_INCLUDE_PATH/duneinterface -DENABLE_ALUGRID"
  CPPFLAGS="$ac_save_CPPFLAGS $ALU3D_INC_FLAG"
  # check for header
  AC_CHECK_HEADERS([alugrid_serial.h], 
     [ALUGRID_CPPFLAGS="$ALU3D_INC_FLAG"
      ALUGRID_LDFLAGS=""
      ALUGRID_LIBS="-L$ALUGRID_LIB_PATH -lalugrid"
    HAVE_ALUGRID="1"],
    AC_MSG_WARN([alugrid_serial.h not found in $ALUGRID_INCLUDE_PATH]))
   
  # Yes, we do check whether either alugrid_serial.h or alugrid_parallel.h
  # works.  Dune decides which one to use depending on how the
  # alugrid_defineparallel.h header defines ALU3DGRID_BUILD_FOR_PARALLEL.
  # This could be improved.
  ALU3D_INC_FLAG_PARA="-I$ALUGRID_INCLUDE_PATH/parallel"
  CPPFLAGS="$ac_save_CPPFLAGS $DUNEMPICPPFLAGS $ALU3D_INC_FLAG_PARA $ALU3D_INC_FLAG"
  # check for parallel header 
  AC_CHECK_HEADERS([alugrid_parallel.h], 
     [ALUGRID_CPPFLAGS="\${DUNEMPICPPFLAGS} $ALU3D_INC_FLAG $ALU3D_INC_FLAG_PARA"
      ALUGRID_LDFLAGS="${DUNEMPILDFLAGS}"
      ALUGRID_LIBS="-L$ALUGRID_LIB_PATH -lalugrid \${DUNEMPILIBS}"
      # for use with the later library test
      LDFLAGS="$LDFLAGS $DUNEMPILDFLAGS"
      LIBS="$DUNEMPILIBS $LIBS"
    HAVE_ALUGRID="1"],
    AC_MSG_WARN([alugrid_parallel.h not found in $ALUGRID_INCLUDE_PATH]))

  # We check only whether linking with the library works, not for an actual
  # function from that library.  So we won't need any special stuff in the
  # CPPFLAGS
  CPPFLAGS="$ac_save_CPPFLAGS"
  # This is a kludge to pass the right libpath before the library on the
  # linker command line.  In the result, the -L flag has to go into the LIBS
  # variable.
  LDFLAGS="$LDFLAGS -L$ALUGRID_LIB_PATH"

  # if header is found...
  if test x$HAVE_ALUGRID = x1 ; then
    AC_CHECK_LIB(alugrid,[malloc],
      [: #dumy argument to avoid default action
      ],
	  [HAVE_ALUGRID="0"
	  AC_MSG_WARN(libalugrid not found!)])
  fi

  AC_LANG_POP([C++])

## end of alugrid check (--without wasn't set)
fi

# survived all tests?
if test x$HAVE_ALUGRID = x1 ; then
  AC_SUBST(ALUGRID_LIBS, $ALUGRID_LIBS)
  AC_SUBST(ALUGRID_LDFLAGS, $ALUGRID_LDFLAGS)
  AC_SUBST(ALUGRID_CPPFLAGS, $ALUGRID_CPPFLAGS)
  AC_DEFINE(HAVE_ALUGRID, ENABLE_ALUGRID,
    [This is only true if alugrid-library was found by configure 
     _and_ if the application uses the ALUGRID_CPPFLAGS])

  # add to global list
  DUNE_ADD_ALL_PKG([ALUGrid], [$ALUGRID_CPPFLAGS],
                   [$ALUGRID_LDFLAGS], [$ALUGRID_LIBS])

  # set variable for summary
  with_alugrid="yes"
else
  AC_SUBST(ALUGRID_LIBS, "")
  AC_SUBST(ALUGRID_LDFLAGS, "")
  AC_SUBST(ALUGRID_CPPFLAGS, "")

  # set variable for summary
  with_alugrid="no"
fi
  
# also tell automake
AM_CONDITIONAL(ALUGRID, test x$HAVE_ALUGRID = x1)

# reset old values
LIBS="$ac_save_LIBS"
CPPFLAGS="$ac_save_CPPFLAGS"
LDFLAGS="$ac_save_LDFLAGS"

DUNE_ADD_SUMMARY_ENTRY([ALUGrid],[$with_alugrid])

])
