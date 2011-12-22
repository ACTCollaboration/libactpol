#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.62)
AC_INIT([libactpol], [0.1], [Mike Nolta <mike@nolta.net>])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/astro.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CC_C99
AC_PROG_INSTALL
AC_PROG_RANLIB

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [floor])
# FIXME: Replace `main' with a function in `-lslalib':
AC_CHECK_LIB([slalib], [slaDat])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor modf pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
