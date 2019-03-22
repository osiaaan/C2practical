# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT
AM_INIT_AUTOMAKE([foreign 1.9 tar-pax])
AC_CONFIG_SRCDIR([src/1/exercise1a.c])
AM_CONFIG_HEADER([config.h])


# this module depends on dune-common dune-grid dune-fem
# this implies checking for [dune-common], [dune-grid], [dune-istl], [dune-fem]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS $BLAS_LIBS"

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/common/Makefile
                 src/problem/Makefile
                 src/1/Makefile
                 src/2/Makefile
                 src/3/Makefile
                 src/4/Makefile
                 src/5/Makefile
                 src/6/Makefile
                 src/7/Makefile
                 src/8/Makefile
                 src/9/Makefile
                 src/10/Makefile
                 doc/Makefile
                 doc/doxygen/Makefile
                 doc/doxygen/pages/Makefile
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
