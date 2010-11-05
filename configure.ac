#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([Open Bicycle Dynamics], [0.01], [hazelnusse@gmail.com])
AC_CONFIG_SRCDIR([src/whipple.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Check for libraries
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

# Checks for header files.
AC_CHECK_HEADERS([getopt.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt])

# Process Makefile.in and src/Makefile.in
AC_CONFIG_FILES([Makefile src/Makefile])

# Generate configure script
AC_OUTPUT
