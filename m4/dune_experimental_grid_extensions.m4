dnl This macro introduces a configure flag --enable-experimental-grid-extensions
dnl that is used to publicly grant access to the implementation of the dune-grid
dnl facades (e.g., Entity, Geometry, etc.).

AC_DEFUN([DUNE_EXPERIMENTAL_GRID_EXTENSIONS],[
  AC_ARG_ENABLE(experimental-grid-extensions,
    AS_HELP_STRING([--enable-experimental-grid-extensions],[If this is set, public access to the implementation of facades like Entity, Geometry, etc. is granted.]))

  AS_IF([test "x$enable_experimental_grid_extensions" = "xyes"],
    AC_DEFINE(DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS, 1, [If this is set, public access to the implementation of facades like Entity, Geometry, etc. is granted.]))
])
