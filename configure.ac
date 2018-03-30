#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.62)
AC_INIT([libactpol], [1.2.0], [Mike Nolta <mike@nolta.net>])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/astro.c])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
LT_INIT

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CC_C99
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([m], [floor])
AC_CHECK_LIB([slalib], [slaDat])
AC_CHECK_LIB([slarefro], [slaf_refro])
AC_CHECK_LIB([sofa_c], [iauUtctai])
AC_CHECK_LIB([wcs], [wcsini])
AC_CHECK_LIB([cfitsio], [ffopen])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_OPENMP

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor modf pow sqrt])

# Turn off dirfile support
AC_ARG_ENABLE([dirfile], AS_HELP_STRING([--disable-dirfile], [Disable dirfile support]))
AM_CONDITIONAL([DIRFILE], [test "x$enable_dirfile" != "xno"])
# Do not compile OLDACT code
AC_ARG_ENABLE([oldact], AS_HELP_STRING([--disable-oldact], [Disable oldact code]))
AM_CONDITIONAL([OLDACT], [test "x$enable_oldact" != "xno"])

# Do not use full SLALIB
AC_ARG_ENABLE([slalib], AS_HELP_STRING([--disable-slalib], [Disable slalib]))
AM_CONDITIONAL([SLALIB], [test "x$enable_slalib" != "xno"])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
