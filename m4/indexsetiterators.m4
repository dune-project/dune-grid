AC_DEFUN([DUNE_INDEXSET_ITERATORS],[
  AC_ARG_ENABLE([indexset-iterators],
    AS_HELP_STRING([--enable-indexset-iterators],
                   [enable deprecated iterators on index sets])
  )
  AS_IF([test x$enable_indexset_iterators = xyes],[
    AC_DEFINE([DUNE_ENABLE_INDEXSET_ITERATORS],[1],[Provide deprecated iterators on index sets])
  ])
])
