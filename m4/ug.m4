## -*- autoconf -*-
# $Id: ug.m4 5156 2008-04-14 09:28:06Z christi $
# searches for UG headers and libs

# DUNE_PATH_UG()
#
# configure shell variables:
#   UGROOT
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
#
# configure substitutions/makefile variables:
#   UG_CPPFLAGS
#   UG_LDFLAGS
#   UG_LIBS
#
# preprocessor defines:
#   HAVE_UG
#     undefined or ENABLE_UG
#
# automake conditionals:
#   UG
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

      # If an explicit path has been provided it needs to be appended
      # temporarily to PKG_CONFIG_PATH
      REM_PKG_CONFIG_PATH=$PKG_CONFIG_PATH
      # The first additional path is for uninstalled UG, the second one for the installed UG
	  UGLIBPKCONFIG=`echo $UGROOT/lib*/pkgconfig | sed -e 's/\s\+/:/g'`
      PKG_CONFIG_PATH="$UGROOT:$UGLIBPKCONFIG:$PKG_CONFIG_PATH"

      UG_LDFLAGS=""

      UG_LIBS="`PKG_CONFIG_PATH=$PKG_CONFIG_PATH $PKG_CONFIG --libs-only-L libug` -lugS2 -lugS3 -ldevS"
      direct_UG_LIBS="`PKG_CONFIG_PATH=$PKG_CONFIG_PATH $PKG_CONFIG --libs-only-L libug` -lugS2 -lugS3 -ldevS"
      
      AC_MSG_CHECKING([for UG])

      # Check whether UG is installed at all
      if PKG_CONFIG_PATH=$PKG_CONFIG_PATH $PKG_CONFIG --exists libug; then
	    HAVE_UG="1"
        AC_MSG_RESULT(yes)
	  else
		HAVE_UG="0"
        AC_MSG_RESULT(no)
		AC_MSG_WARN([UG not found])
      fi

      ## check version number 
      NEEDEDUG_VERSION=3.9.1-patch7

      if test x$HAVE_UG = x1; then
          
          AC_MSG_CHECKING([whether UG version is recent enough])

          # Does it have a suitable version?
          if PKG_CONFIG_PATH=$PKG_CONFIG_PATH $PKG_CONFIG --atleast-version=$NEEDEDUG_VERSION libug; then
              AC_MSG_RESULT(yes)
          else
              HAVE_UG="0"
              AC_MSG_RESULT(no)
              AC_MSG_WARN([UG version is too old (you need at least $NEEDEDUG_VERSION)])
          fi
      fi

      # pre-set variable for summary
      with_ug="no"
   
      if test x$HAVE_UG = x1; then

        # Okay.  We have found a UG installation.  But has it been built with --enable-dune?
        if test x$HAVE_UG = x1 ; then
              
          AC_MSG_CHECKING([whether UG has been built with --enable-dune])

          if test x`PKG_CONFIG_PATH=$PKG_CONFIG_PATH $PKG_CONFIG --variable=fordune libug` == xyes; then
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

        # Set the compiler flags
		UG_CPPFLAGS="`PKG_CONFIG_PATH=$PKG_CONFIG_PATH $PKG_CONFIG --cflags-only-I libug` -DENABLE_UG"
        direct_UG_CPPFLAGS="`PKG_CONFIG_PATH=$PKG_CONFIG_PATH $PKG_CONFIG --cflags-only-I libug` -DENABLE_UG"

          if test x`PKG_CONFIG_PATH=$PKG_CONFIG_PATH $PKG_CONFIG --variable=parallel libug` == xyes; then
			
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
  fi 
      
  # tell automake   
  AM_CONDITIONAL(UG, test x$HAVE_UG = x1)
  
  # restore variables
  LDFLAGS="$ac_save_LDFLAGS"
  CPPFLAGS="$ac_save_CPPFLAGS"
  LIBS="$ac_save_LIBS"

  DUNE_ADD_SUMMARY_ENTRY([UG],[$with_ug])
  
])
