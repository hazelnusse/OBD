AC_PREREQ([2.67])
AC_INIT([Open Bicycle Dynamics], [0.02], [hazelnusse@gmail.com])
AM_INIT_AUTOMAKE([dist-bzip2 -Wall])
#AM_SILENT_RULES([yes])

AM_PATH_PYTHON([2.6])
AC_CONFIG_SRCDIR([whipple/whippleeig.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB

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

# All files that will be generated from templates
AC_CONFIG_FILES([Makefile
                 whipple/common/Makefile
                 whipple/tests/Makefile
                 whipple/Makefile])

# Generate configure script
AC_OUTPUT
