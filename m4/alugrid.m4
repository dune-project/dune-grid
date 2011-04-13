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

  AC_ARG_WITH([alugrid-libdir],dnl
    AS_HELP_STRING([--with-alugrid-libdir=PATH],dnl
      [Directory where ALUGrid library is installed (deprecated). Note that this will override library path detection, so use this parameter only if default library detection fails and you know exactly where your ALUGrid library is located.]))dnl


# do not use alugrid debug lib 

# store old values
ac_save_LDFLAGS="$LDFLAGS"
ac_save_CPPFLAGS="$CPPFLAGS"
ac_save_LIBS="$LIBS"

# initilize to sane value
HAVE_ALUGRID=0

## do nothing if no --with-alugrid was supplied
if test x$with_alugrid != xno ; then

  # is --with-alugrid=PATH used?
  AS_IF([test "x$with_alugrid" != "x"],[
    AS_IF([test -d $with_alugrid],[
      AC_MSG_NOTICE([searching for ALUGrid in $with_alugrid...])
      ALUGRIDROOT=`cd $with_alugrid && pwd`
    ],[
      AC_MSG_WARN([ALUGrid directory '$with_alugrid' does not exist or is inaccessible])
    ])
  ],[
    # educated guess for alugrid root
    for d in /usr /usr/local /usr/local/alugrid /opt/alugrid; do
      AC_MSG_NOTICE([searching for ALUGrid in $d...])
      AS_IF([test -f $d/lib/pkgconfig/alugrid.pc -o -x $d/bin/alugridversion],[
        ALUGRIDROOT="$d"
        break
      ])
    done
  ])

  REM_PKG_CONFIG_PATH=$PKG_CONFIG_PATH
  PKG_CONFIG_PATH="$ALUGRIDROOT:$ALUGRIDROOT/lib/pkgconfig:$PKG_CONFIG_PATH"

  ## check version number 
  NEEDEDALUGRID_VERSION=1.23

  AC_MSG_CHECKING([ALUGrid version >= $NEEDEDALUGRID_VERSION])
  if $PKG_CONFIG --atleast-version=$NEEDEDALUGRID_VERSION alugrid ; then 
    ALUGRID_VERSION=`$PKG_CONFIG --modversion alugrid`
    AC_MSG_RESULT([yes (ALUGrid-$ALUGRID_VERSION)])
    ALUGRID_INCLUDE_PATH=`$PKG_CONFIG --variable=includedir alugrid`
    ALUGRID_LIB_PATH=`$PKG_CONFIG --variable=libdir alugrid`
    if test x"$with_alugrid_libdir" != x"" && test x"$with_alugrid_libdir" != x"no" ; then 
      AC_MSG_WARN([--with-alugrid-libdir option is obsolete!])
    fi 
  else   
    # lib dir and include path 
    ALUGRID_INCLUDE_PATH="$ALUGRIDROOT/include"

    ALUGRID_LIB_PATH="$ALUGRIDROOT/lib"
    if test x"$with_alugrid_libdir" != x"" && test x"$with_alugrid_libdir" != x"no"
    then
      if ! test -d "$with_alugrid_libdir"
      then
        AC_MSG_ERROR([library directory $with_alugrid_libdir for ALUGrid does not exist or is inaccessible.])dnl
      else
        ALUGRID_LIB_PATH="$with_alugrid_libdir"
      fi
    fi

    # old check version 
    ALUGRID_VERSIONCHECK=$ALUGRIDROOT/bin/alugridversion
    if test -f $ALUGRID_VERSIONCHECK; then 
      ALUGRID_VERSION=`$ALUGRID_VERSIONCHECK -c $NEEDEDALUGRID_VERSION`
      if test "x$ALUGRID_VERSION" != "x-1"; then 
        ALUGRID_VERSION=`$ALUGRID_VERSIONCHECK -v`
        AC_MSG_RESULT([yes (ALUGrid-$ALUGRID_VERSION)])
      else 
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([ALUGrid version is too old!])
      fi
    else 
      AC_MSG_RESULT([no])
      AS_IF([test "x$ALUGRIDROOT" != "x"],[
        AC_MSG_ERROR([$PKG_CONFIG couldn't find alugrid.pc or wrong version! ALUGrid version is too old or ALUGrid is not installed in $ALUGRIDROOT! You need at least ALUGrid-$NEEDEDALUGRID_VERSION!])
      ])
    fi
  fi

  # restore PKG_CONFIG_PATH 
  PKG_CONFIG_PATH=$REM_PKG_CONFIG_PATH

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

  HAVE_ALUGRID_PARALLEL="0"
  # if the serial header was found then also check for the parallel header
  if test x"$HAVE_ALUGRID" != "x0" ; then
    AC_TRY_COMPILE([#include <alugrid_defineparallel.h> 
                    #if ALU3DGRID_BUILD_FOR_PARALLEL == 0 
                    #error
                    #endif
                   ],
                   [],
                   [HAVE_ALUGRID_PARALLEL="1"],
                   [HAVE_ALUGRID_PARALLEL="0"
                   AC_MSG_WARN("ALUGRID was not build for parallel support!")
                  ])

    # only check for parallel header when ALUGrid was build for parallel support
    if test x"$HAVE_ALUGRID_PARALLEL" != "x0" ; then 
      # check for parallel header 
      AC_CHECK_HEADERS([alugrid_parallel.h], 
         [ALUGRID_CPPFLAGS="\${DUNEMPICPPFLAGS} $ALU3D_INC_FLAG $ALU3D_INC_FLAG_PARA"
          ALUGRID_LDFLAGS="\${DUNEMPILDFLAGS}"
          ALUGRID_LIBS="-L$ALUGRID_LIB_PATH -lalugrid \${DUNEMPILIBS}"
          # for use with the later library test
          LDFLAGS="$LDFLAGS $DUNEMPILDFLAGS"
          LIBS="$DUNEMPILIBS $LIBS"
        HAVE_ALUGRID="1"],
        AC_MSG_WARN([alugrid_parallel.h not found in $ALUGRID_INCLUDE_PATH]))
    fi 
  fi

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

  DUNE_DEFINE_GRIDTYPE([ALUGRID_CONFORM],[],[Dune::ALUConformGrid< dimgrid, dimworld >],[dune/grid/alugrid.hh],[dune/grid/io/file/dgfparser/dgfalu.hh])
  DUNE_DEFINE_GRIDTYPE([ALUGRID_CUBE],[],[Dune::ALUCubeGrid< dimgrid, dimworld >],[dune/grid/alugrid.hh],[dune/grid/io/file/dgfparser/dgfalu.hh])
  DUNE_DEFINE_GRIDTYPE([ALUGRID_SIMPLEX],[],[Dune::ALUSimplexGrid< dimgrid, dimworld >],[dune/grid/alugrid.hh],[dune/grid/io/file/dgfparser/dgfalu.hh])

  # set variable for summary
  if test x"$HAVE_ALUGRID_PARALLEL" != "x0" ; then 
    with_alugrid_parallel="(parallel)"
  else 
    with_alugrid_parallel="(serial)"
  fi 
  with_alugrid="version $ALUGRID_VERSION $with_alugrid_parallel"
  with_alugrid_long="$ALUGRIDROOT"
else
  AC_SUBST(ALUGRID_LIBS, "")
  AC_SUBST(ALUGRID_LDFLAGS, "")
  AC_SUBST(ALUGRID_CPPFLAGS, "")

  # set variable for summary
  with_alugrid="no"
  with_alugrid_long=""
fi
  
# also tell automake
AM_CONDITIONAL(ALUGRID, test x$HAVE_ALUGRID = x1)

# reset old values
LIBS="$ac_save_LIBS"
CPPFLAGS="$ac_save_CPPFLAGS"
LDFLAGS="$ac_save_LDFLAGS"

DUNE_ADD_SUMMARY_ENTRY([ALUGrid],[$with_alugrid],[$with_alugrid_long])

])
