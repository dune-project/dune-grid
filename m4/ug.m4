## -*- autoconf -*-
# $Id: ug.m4 5156 2008-04-14 09:28:06Z christi $
# searches for UG headers and libs

# DUNE_PATH_UG()
#
# configure shell variables:
#   UGROOT
#   UG_LIB_PATH
#   UG_INCLUDE_PATH
#   UG_CPPFLAGS, UG_LDFLAGS, UG_LIBS
#       flags and libs with indirect references for the makefiles, for
#       instance the literal string '${DUNEMPICPPFLAGS}
#   direct_UG_CPPFLAGS, direct_UG_LDFLAGS, direct_UG_LIBS
#       flags and libs with direct values for use in configure, for instance
#       the value of DUNEMPICPPFLAGS
#   HAVE_UG
#       1 or 0 or undefined
#   with_ug
#       "no" or "yes" with stuff appended
#   enable_ug_lgmdomain
#       yes or no
#
# configure substitutions/makefile variables:
#   UG_CPPFLAGS
#   UG_LDFLAGS
#   UG_LIBS
#
# preprocessor defines:
#   HAVE_UG
#     undefined or ENABLE_UG
#   UG_LGMDOMAIN
#     undefined or 1
#
# automake conditionals:
#   UG
#   UG_LGMDOMAIN
AC_DEFUN([DUNE_PATH_UG],[
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PATH_XTRA])
  AC_REQUIRE([DUNE_MPI])

  AC_ARG_WITH(ug,
    AC_HELP_STRING([--with-ug=PATH],[directory where UG is installed]))

  # store old values
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LIBS="$LIBS"
  
  # initialize
  HAVE_UG=0

  ## do nothing if --without-ug is used
  if test x$with_ug != xno ; then
      
      # is --with-ug=bla used?
      if test "x$with_ug" != x ; then
          if ! test -d $with_ug; then
              AC_MSG_WARN([UG directory $with_ug does not exist!])
          else
              # expand tilde / other stuff
              UGROOT=`cd $with_ug && pwd`
          fi
      fi
      if test "x$UGROOT" = x; then
          # use some default value...
          UGROOT="/usr/local/ug"
      fi
      
      # intermediate variables
      UG_LIB_PATH="$UGROOT/lib"
      UG_INCLUDE_PATH="$UGROOT/include/ug/"
      
      UG_LDFLAGS=""

      # set variables so that tests can use them
      CPPFLAGS="$CPPFLAGS -I$UG_INCLUDE_PATH -DENABLE_UG"
      # hack around limitation of AC_CHECK_LIBS: -L really belong into LIBS,
      # but it has to be in front of the library that AC_CHECK_LIBS inserts on
      # the linker command line
      LDFLAGS="$LDFLAGS -L$UG_LIB_PATH"

      AC_ARG_ENABLE(ug-lgmdomain,
        AC_HELP_STRING([--enable-ug-lgmdomain],[use UG LGM domain (default is standard domain)]))
      if test x"$enable_ug_lgmdomain" = xyes ; then
        UG_LIBS="-L$UG_LIB_PATH -lugL2 -lugL3 -ldevS"
        direct_UG_LIBS="-L$UG_LIB_PATH -lugL2 -lugL3 -ldevS"
      else
        UG_LIBS="-L$UG_LIB_PATH -lugS2 -lugS3 -ldevS"
        direct_UG_LIBS="-L$UG_LIB_PATH -lugS2 -lugS3 -ldevS"
      fi

      REM_PKG_CONFIG_PATH=$PKG_CONFIG_PATH
      PKG_CONFIG_PATH="$UGROOT/lib/pkgconfig:$PKG_CONFIG_PATH"

      # Check whether UG is installed at all and has a suitable version
      if $PKG_CONFIG --atleast-version=3.9.2 libug; then
		UG_CPPFLAGS="-I$UG_INCLUDE_PATH -DENABLE_UG"
		direct_UG_CPPFLAGS="-I$UG_INCLUDE_PATH -DENABLE_UG"
	    HAVE_UG="1"
	  else
		HAVE_UG="0"
		AC_MSG_WARN([UG not found in $UGROOT])
      fi

      # pre-set variable for summary
      with_ug="no"
      
      if test x$HAVE_UG = x1; then

        # Okay.  We have found a UG installation.  But has it been built with --enable-dune?
        if test x$HAVE_UG = x1 ; then
              
          AC_MSG_CHECKING([whether UG has been built with --enable-dune])

          if test x`$PKG_CONFIG --variable=fordune libug` == xyes; then
              AC_MSG_RESULT(yes)
          else
              AC_MSG_RESULT(no)
              AC_MSG_WARN([UG has not been built with --enable-dune!])
              HAVE_UG="0"
              with_ug="no"
          fi
            
        fi
        
      fi

      if test x$HAVE_UG = x1; then

		if test x`$PKG_CONFIG --variable=parallel libug` == xyes; then
			
          # Add additional flags needed for parallel UG  
		  UG_LDFLAGS="\${DUNEMPILDFLAGS} $UG_LDFLAGS"
          direct_UG_LDFLAGS="$DUNEMPILDFLAGS $direct_UG_LDFLAGS"
          UG_CPPFLAGS="\${DUNEMPICPPFLAGS} $UG_CPPFLAGS -DModelP"
          direct_UG_CPPFLAGS="$DUNEMPICPPFLAGS $direct_UG_CPPFLAGS -DModelP"
          UG_LIBS="$UG_LIBS \${DUNEMPILIBS}"
          direct_UG_LIBS="$direct_UG_LIBS $DUNEMPILIBS"
          with_ug="yes (parallel)"
		   
		else
			
          with_ug="yes (sequential)"
				   
	    fi

      fi

      # restore PKG_CONFIG_PATH 
      PKG_CONFIG_PATH=$REM_PKG_CONFIG_PATH
  
  # end of "no --without-ug"
  fi

  # did it work?
  if test x$HAVE_UG = x0 ; then
      # reset flags, so they do not appear in makefiles
      UG_CPPFLAGS=""
      direct_UG_CPPFLAGS=""
      UG_LDFLAGS=""
      direct_UG_LDFLAGS=""
      UG_LIBS=""
      direct_UG_LIBS=""
  fi

  AC_SUBST([UG_LDFLAGS])
  AC_SUBST([UG_LIBS])
  AC_SUBST([UG_CPPFLAGS])

  # add to global list
  DUNE_ADD_ALL_PKG([UG], [\${UG_CPPFLAGS}], [\${UG_LDFLAGS}], [\${UG_LIBS}])

  if test x$HAVE_UG = x1 ; then

      # add support for GRIDTYPE=UGGRID to config.h
      DUNE_DEFINE_GRIDTYPE([UGGRID],[GRIDDIM == WORLDDIM],[Dune::UGGrid< dimgrid >],[dune/grid/uggrid.hh],[dune/grid/io/file/dgfparser/dgfug.hh])

      AC_DEFINE(HAVE_UG, ENABLE_UG, 
        [This is only true if UG was found by configure 
         _and_ if the application uses the UG_CPPFLAGS])
      if test x"$enable_ug_lgmdomain" = xyes ; then
        AC_DEFINE(UG_LGMDOMAIN, 1, [use UG LGM domain])
      fi
  fi 
      
  # tell automake   
  AM_CONDITIONAL(UG, test x$HAVE_UG = x1)
  AM_CONDITIONAL(UG_LGMDOMAIN, test x$HAVE_UG = x1 && test x$UG_LGMDOMAIN = x1)
  
  # restore variables
  LDFLAGS="$ac_save_LDFLAGS"
  CPPFLAGS="$ac_save_CPPFLAGS"
  LIBS="$ac_save_LIBS"

  DUNE_ADD_SUMMARY_ENTRY([UG],[$with_ug])
  
])
