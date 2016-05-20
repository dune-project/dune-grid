## -*- autoconf -*-
# searches for alberta-headers and libs

# Substitutes the following make variables:
#   ALBERTA_DIM
#     defaults to $WORLD_DIM
#
#   ALBERTAROOT = /usr/local/alberta
#     Root dir of the alberta installation.  Set from --with-alberta=...
#
#   ALBERTA_BASE_LIBS = $(ALBERTA_LIBPATHFLAGS) -lalberta_util $ALBERTA_EXTRA
#     LIBS that are always required independent of dimension
#
#   ALBERTA_LIBPATHFLAGS = -L$(ALBERTAROOT)/lib
#     Library path required for alberta
#
#   ALBERTA%DIM%D_LIBS = -L$(DUNE_GRID_LIBDIR) -ldunealbertagrid_%DIM%d -ldunegrid \
#              $(ALBERTA_LIBPATHFLAGS) -lalberta_%DIM%d \
#              $(ALBERTA_BASE_LIBS)
#     *OR*       = $(top_builddir)/lib/libdunealbertagrid_%DIM%d.la \
#              $(top_builddir)/lib/libdunegrid.la \
#              $(ALBERTA_LIBPATHFLAGS) -lalberta_%DIM%d \
#              $(ALBERTA_BASE_LIBS)
#     All LIBS required for dimension %DIM% (1, 2, or 3).  The first value is
#     substituted by default and is appropriate for modules depending on
#     dune-grid.  dune-grid itself will overwrite that with the second value
#     however in configure.ac.
#
#   ALBERTA_LIBS = $(ALBERTA($ALBERTA_DIM)D_LIBS)
#     All LIBS required for the configured dimension.  The value of this
#     variable will be empty for dimensions other than 1, 2, or 3.
#
#   ALBERTA_INCLUDE_CPPFLAGS = -I$(ALBERTAROOT)/include/alberta
#     Include path required for Alberta.
#
#   ALBERTA%DIM%D_CPPFLAGS = $(ALBERTA_INCLUDE_CPPFLAGS) \
#              -DALBERTA_DIM=%DIM% -DENABLE_ALBERTA
#     All CPPFLAGS required for dimension %DIM% (1, 2, or 3).
#
#   ALBERTA_CPPFLAGS = $(ALBERTA$(ALBERTA_DIM)D_CPPFLAGS)
#     All CPPFLAGS required for the configured dimension.  The value of this
#     variable will be empty for dimensions other than 1, 2, or 3, thus
#     ENABLE_ALBERTA will not be defined inside the program, disabling alberta
#     support.
#
#   ALBERTA%DIM%D_LDFLAGS =
#     All LDFLAGS required for dimension %DIM% (1, 2, or 3).  These are
#     currently empty and exist just for consistency.
#
#   ALBERTA_LDFLAGS = $(ALBERTA$(ALBERTA_DIM)_LDFLAGS)
#     All LDFLAGS required for the configured dimension.  These are currently
#     empty and exist just for consistency.
#
#   If you want to use the configured dimension, you have to use
#   $(ALBERTA_LIBS), $(ALBERTA_CPPFLAGS) and $(ALBERTA_LDFLAGS).  If the
#   configured dimension is anything other than 1, 2, or 3, these variable
#   will substitute empty values, thus disabling support for alberta in the
#   program.
#
#   If want to use a specific dimension, say 2, you have to use
#   $(ALBERTA2D_LIBS), $(ALBERTA2D_CPPFLAGS) and $(ALBERTA2D_LDFLAGS).
#
# Defines the folling CPP macro
#   ALBERTA_DIM
#     The Alberta dimension this binary will be linked with.
#   DUNE_ALBERTA_VERSION
#     Alberta version found by configure, always 0x300 for 3.0
#   HAVE_ALBERTA
#     This is only true if alberta-library was found by configure 
#     _and_ if the application uses the ALBERTA_CPPFLAGS
#
# Defines the following automake conditional
#    ALBERTA
#
# configure shell variables:
#    HAVE_ALBERTA
#      1 if a working Alberta was found.
AC_DEFUN([DUNE_PATH_ALBERTA],[
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PROG_F77])
  AC_REQUIRE([AC_PATH_XTRA])
  AC_REQUIRE([DUNE_PATH_OPENGL])

  ALBERTA_DIM='$(WORLDDIM)'

  AC_ARG_WITH(alberta,
    AC_HELP_STRING([--with-alberta=PATH],[directory where ALBERTA (ALBERTA
    version 3.0) is installed.  You can pass additional required
    libraries in the ALBERTA_EXTRA environment variable (in a form suitable
    for $LIBS)]))

  AC_ARG_ENABLE([alberta-libcheck],
    AS_HELP_STRING([--disable-alberta-libcheck],
      [Do not try to link against libalberta_Nd.]))

  # do not use alberta debug lib 
  with_alberta_debug=0

  # store old values
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LIBS="$LIBS"
  # LIBS=""

  ## do nothing if no --with-alberta was supplied
  AS_IF([test x$with_alberta != xno],[

    # is --with-alberta=PATH used?
    AS_IF([test "x$with_alberta" != "x"],[
      AS_IF([test -d $with_alberta],[
        AC_MSG_NOTICE([searching for ALBERTA in $with_alberta...])
        ALBERTAROOT=`cd $with_alberta && pwd`
      ],[
        AC_MSG_WARN([ALBERTA directory '$with_alberta' does not exist])
      ])
    ],[
      # educated guess for alberta root
      for d in /usr /usr/local /usr/local/alberta /opt/alberta; do
        AC_MSG_NOTICE([searching for ALBERTA in $d...])
        AS_IF([test -d $d/include/alberta],[
          ALBERTAROOT="$d"
          break
        ])
      done
    ])


    # set variables so that tests can use them
    ALBERTA_INCLUDE_CPPFLAGS="-I$ALBERTAROOT/include -I$ALBERTAROOT/include/alberta"

    # define varaible flags depending on problem and world dim, to change afterwards easily
    ALBERTA_CPPFLAGS='$(ALBERTA$(ALBERTA_DIM)D_CPPFLAGS)'
    ALBERTA_LDFLAGS='$(ALBERTA$(ALBERTA_DIM)D_LDFLAGS)'
    ALBERTA_LIBS='$(ALBERTA$(ALBERTA_DIM)D_LIBS)'

    # initialize dimension dependent CPPFLAGS, LDFLAGS and LIBS to default values
    for N in 1 2 3 4 5 6 7 8 9 ; do
      eval ALBERTA${N}D_CPPFLAGS=
      eval ALBERTA${N}D_LDFLAGS=
      eval ALBERTA${N}D_LIBS=
    done

    # check for header
    CPPFLAGS="$ac_save_CPPFLAGS $ALBERTA_INCLUDE_CPPFLAGS -DDIM_OF_WORLD=3 -DEL_INDEX=0"
    AC_CHECK_HEADER([alberta/alberta.h], [HAVE_ALBERTA="1"],[
      AC_MSG_WARN([alberta/alberta.h not found in $ALBERTA_INCLUDE_CPPFLAGS])
    ])

    if test "x$HAVE_ALBERTA" = "x1" ; then
      AC_CHECK_MEMBER([struct el_info.wall_bound],[ALBERTA_VERSION="3.0"],
        [AC_MSG_WARN([version 3 not found, deactivating Alberta])
          HAVE_ALBERTA="0"],
        [#include <alberta/alberta.h>])
    fi

    CPPFLAGS="$ac_save_CPPFLAGS $ALBERTA_INCLUDE_CPPFLAGS"

    # TODO: check if static flag exists 
    # link_static_flag defines the flag for the linker to link only static
    # didn't work, with $link_static_flag, so quick hack here

    # check for libalberta_util...
    AS_IF([test "x$HAVE_ALBERTA" = "x1"],[
      AC_CACHE_CHECK([ALBERTA utilities library],[dune_grid_cv_lib_alberta_utils],[
        dune_grid_cv_lib_alberta_utils=no
        for alberta_lib_dir in lib lib64 ; do
          for lib_alberta_utils in alberta_utilities alberta_util ; do
            ALBERTA_LIBPATHFLAGS="-L$ALBERTAROOT/$alberta_lib_dir"
            DUNEALBERTA_LIBPATHFLAGS="-L$top_builddir/$alberta_lib_dir"
            LDFLAGS="$LDFLAGS -L$ALBERTAROOT/$alberta_lib_dir"
            LIBS="-l$lib_alberta_utils $ALBERTA_EXTRA $ac_save_LIBS"
            AC_LINK_IFELSE(
              [AC_LANG_CALL([], [alberta_calloc])],
              [dune_grid_cv_lib_alberta_utils=$lib_alberta_utils; 
                break],
              [])
          done
          if test "x$dune_grid_cv_lib_alberta_utils" != "xno"; then break; fi
        done
      ])
      AS_IF([test "x$dune_grid_cv_lib_alberta_utils" = "xno"],[HAVE_ALBERTA=0])
    ])

    # check for ALBERTA grid library...
    AS_IF([test "x$HAVE_ALBERTA" = "x1"],[
      AS_IF([test x$enable_alberta_libcheck = xno],[
        AC_MSG_WARN([Disabled checking whether libalberta_Nd can be linked.])
      ],[
        AC_CACHE_CHECK([ALBERTA world dimensions],[dune_grid_cv_alberta_world_dims],[
          dune_grid_cv_alberta_world_dims=
          for N in 1 2 3 4 5 6 7 8 9; do
            LIBS="-lalberta_${N}d -l$dune_grid_cv_lib_alberta_utils $ALBERTA_EXTRA $ac_save_LIBS"
            AC_LINK_IFELSE(
              [AC_LANG_CALL([], [mesh_traverse])],
              [dune_grid_cv_alberta_world_dims="$dune_grid_cv_alberta_world_dims $N"],
              [LIBS="-lalberta_${N}d -l$dune_grid_cv_lib_alberta_utils $ALBERTA_EXTRA -lltdl $ac_save_LIBS"
                AC_LINK_IFELSE(
                  [AC_LANG_CALL([], [mesh_traverse])],
                  [dune_grid_cv_alberta_world_dims="$dune_grid_cv_alberta_world_dims $N"
                    ALBERTA_EXTRA="$ALBERTA_EXTRA -lltdl"],
                  [])
              ])
          done
        ])
        ALBERTA_WORLD_DIMS="$dune_grid_cv_alberta_world_dims"
      ])

      AS_IF([test "x$ALBERTA_WORLD_DIMS" != "x"],[
        ALBERTA_BASE_LIBS="\$(ALBERTA_LIBPATHFLAGS) -l$dune_grid_cv_lib_alberta_utils $ALBERTA_EXTRA"

        # define library variables for all found libraries
        for N in $ALBERTA_WORLD_DIMS ; do
          eval ALBERTA${N}D_CPPFLAGS="'\$(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=${N} -DENABLE_ALBERTA'"
          eval ALBERTA${N}D_LDFLAGS=
          # Dune was installed into directory given by with-dunegrid
          eval ALBERTA${N}D_LIBS="'-L\$(DUNE_GRID_LIBDIR) -ldunealbertagrid_${N}d -ldunegrid \$(ALBERTA_LIBPATHFLAGS) -lalberta_${N}d \$(ALBERTA_BASE_LIBS)'"
        done
      ],[
        HAVE_ALBERTA=0
      ])

    ])

  ])  # end of alberta check (--without wasn't set)

  # survived all tests?
  AS_IF([test "x$HAVE_ALBERTA" = "x1"],[

    AC_DEFINE(HAVE_ALBERTA, ENABLE_ALBERTA,
      [This is only true if alberta-library was found by configure 
       _and_ if the application uses the ALBERTA_CPPFLAGS])

    if test "$ALBERTA_VERSION" = "3.0" ; then
      AC_DEFINE([DUNE_ALBERTA_VERSION], [0x300], [Alberta version found by configure, should be 0x300 for 3.0])
    else
      AC_MSG_ERROR([Internal Inconsistency: Invalid Alberta version reported: $ALBERTA_VERSION.])
    fi

    # add to global list
    DUNE_ADD_ALL_PKG([Alberta], [\${ALBERTA_CPPFLAGS}],
                     [\${ALBERTA_LDFLAGS}], [\${ALBERTA_LIBS}])

    DUNE_DEFINE_GRIDTYPE([ALBERTAGRID],[WORLDDIM == ALBERTA_DIM],[Dune::AlbertaGrid< dimgrid >],[dune/grid/albertagrid.hh],[dune/grid/albertagrid/dgfparser.hh])

    # set variable for summary
    with_alberta="version $ALBERTA_VERSION"
    with_alberta_long="$ALBERTAROOT ; world dims $ALBERTA_WORLD_DIMS"

  ],[

    # clear all variables
    ALBERTA_DIM= 
    ALBERTAROOT= 
    ALBERTA_BASE_LIBS=
    DUNEALBERTA_LIBPATHFLAGS=
    ALBERTA_LIBPATHFLAGS=
    ALBERTA_INCLUDE_CPPFLAGS=
    ALBERTA_DIM_CPPFLAGS=
    ALBERTA_CPPFLAGS=
    ALBERTA_LDFLAGS=
    ALBERTA_LIBS=

    ALBERTA_WORLD_DIMS=

    for N in 1 2 3 4 5 6 7 8 9 ; do
      eval ALBERTA${N}D_CPPFLAGS=
      eval ALBERTA${N}D_LDFLAGS=
      eval ALBERTA${N}D_LIBS=
    done

    # set variable for summary
    with_alberta="no"
    with_alberta_long=""

  ])
    
  AC_SUBST([ALBERTA_DIM]) 
  AC_SUBST([ALBERTAROOT]) 
  AC_SUBST([ALBERTA_BASE_LIBS])
  AC_SUBST([DUNEALBERTA_LIBPATHFLAGS])
  AC_SUBST([ALBERTA_LIBPATHFLAGS])
  AC_SUBST([ALBERTA_INCLUDE_CPPFLAGS])
  AC_SUBST([ALBERTA_DIM_CPPFLAGS])
  AC_SUBST([ALBERTA_CPPFLAGS])
  AC_SUBST([ALBERTA_LDFLAGS])
  AC_SUBST([ALBERTA_LIBS])

  AC_SUBST([ALBERTA1D_CPPFLAGS])
  AC_SUBST([ALBERTA2D_CPPFLAGS])
  AC_SUBST([ALBERTA3D_CPPFLAGS])
  AC_SUBST([ALBERTA4D_CPPFLAGS])
  AC_SUBST([ALBERTA5D_CPPFLAGS])
  AC_SUBST([ALBERTA6D_CPPFLAGS])
  AC_SUBST([ALBERTA7D_CPPFLAGS])
  AC_SUBST([ALBERTA8D_CPPFLAGS])
  AC_SUBST([ALBERTA9D_CPPFLAGS])

  AC_SUBST([ALBERTA1D_LDFLAGS])
  AC_SUBST([ALBERTA2D_LDFLAGS])
  AC_SUBST([ALBERTA3D_LDFLAGS])
  AC_SUBST([ALBERTA4D_LDFLAGS])
  AC_SUBST([ALBERTA5D_LDFLAGS])
  AC_SUBST([ALBERTA6D_LDFLAGS])
  AC_SUBST([ALBERTA7D_LDFLAGS])
  AC_SUBST([ALBERTA8D_LDFLAGS])
  AC_SUBST([ALBERTA9D_LDFLAGS])

  AC_SUBST([ALBERTA1D_LIBS])
  AC_SUBST([ALBERTA2D_LIBS])
  AC_SUBST([ALBERTA3D_LIBS])
  AC_SUBST([ALBERTA4D_LIBS])
  AC_SUBST([ALBERTA5D_LIBS])
  AC_SUBST([ALBERTA6D_LIBS])
  AC_SUBST([ALBERTA7D_LIBS])
  AC_SUBST([ALBERTA8D_LIBS])
  AC_SUBST([ALBERTA9D_LIBS])

  # also tell automake
  AM_CONDITIONAL(ALBERTA, test x$HAVE_ALBERTA = x1)

  AM_CONDITIONAL(ALBERTA_1D,[test ! -z "$(echo $ALBERTA_WORLD_DIMS | grep 1)"])
  AM_CONDITIONAL(ALBERTA_2D,[test ! -z "$(echo $ALBERTA_WORLD_DIMS | grep 2)"])
  AM_CONDITIONAL(ALBERTA_3D,[test ! -z "$(echo $ALBERTA_WORLD_DIMS | grep 3)"])
  AM_CONDITIONAL(ALBERTA_4D,[test ! -z "$(echo $ALBERTA_WORLD_DIMS | grep 4)"])
  AM_CONDITIONAL(ALBERTA_5D,[test ! -z "$(echo $ALBERTA_WORLD_DIMS | grep 5)"])
  AM_CONDITIONAL(ALBERTA_6D,[test ! -z "$(echo $ALBERTA_WORLD_DIMS | grep 6)"])
  AM_CONDITIONAL(ALBERTA_7D,[test ! -z "$(echo $ALBERTA_WORLD_DIMS | grep 7)"])
  AM_CONDITIONAL(ALBERTA_8D,[test ! -z "$(echo $ALBERTA_WORLD_DIMS | grep 8)"])
  AM_CONDITIONAL(ALBERTA_9D,[test ! -z "$(echo $ALBERTA_WORLD_DIMS | grep 9)"])

  # reset old values
  LIBS="$ac_save_LIBS"
  CPPFLAGS="$ac_save_CPPFLAGS"
  LDFLAGS="$ac_save_LDFLAGS"

  DUNE_ADD_SUMMARY_ENTRY([ALBERTA],[$with_alberta],[$with_alberta_long])
])
