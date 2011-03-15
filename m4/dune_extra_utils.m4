## -*- autoconf -*-
AC_DEFUN([DUNE_EXTRA_UTILS], [
  AC_ARG_ENABLE([extra-utils],
    [AS_HELP_STRING(
        [--enable-extra-utils],
        [Enable compilation and installation of extra utilities from the
          "src" subdirectory.])],
    [], [enable_extra_utils=no])
  case "$enable_extra_utils" in
  yes|no)
    DUNE_ADD_SUMMARY_ENTRY([Extra utils enabled], [$enable_extra_utils]);;
  *) AC_MSG_ERROR([invalid argument for --enable-extra-utils: $enable_extra_utils]);;
  esac
  AM_CONDITIONAL([EXTRA_UTILS], [test "x$enable_extra_utils" = xyes])
])

