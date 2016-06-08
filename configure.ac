# Copyright 2016, Alistair Boyle, 3-clause BSD License
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([zedhat],[0.0.1],[alistair.js.boyle@gmail.com])

AC_COPYRIGHT([Copyright 2016, Alistair Boyle, 3-clause BSD License])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])

AC_CONFIG_SRCDIR([src/zedhat.c])
AC_CONFIG_HEADERS([src/config.h])
AC_USE_SYSTEM_EXTENSIONS

# debugging
CFLAGS="$CFLAGS -Wall -Werror"
AC_MSG_CHECKING(whether to enable debugging)
default="no"
DEBUG=
AC_ARG_ENABLE(debug, [  --enable-debug=[no/yes] turn on debugging
                       [default=$default]],, enable_debug=$default)
if test "x$enable_debug" = "xyes"; then
  CFLAGS="$CFLAGS -g -DDEBUG"
  DEBUG="${DEBUG} -g";
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi
# profiling
AC_MSG_CHECKING(whether to enable profiling)
default="no"
AC_ARG_ENABLE(debug, [  --enable-profile=[no/yes] turn on profiling
                       [default=$default]],, enable_profile=$default)
if test "x$enable_profile" = "xyes"; then
  CFLAGS="${CFLAGS} -g -DDEBUG -pg"
  CXXFLAGS="${CXXFLAGS} -g -DDEBUG -pg"
  DEBUG="${DEBUG} -pg";
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi
# .. and export so Makefile will use them
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(DEBUG)


# makefile
AM_INIT_AUTOMAKE([1.11 subdir-objects]) # need at least v1.11 of automake for SILENT_RULES
AM_SILENT_RULES([yes]) # enable kernel like cleaner build
AC_CONFIG_FILES([Makefile])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O

# Checks for libraries.
#MC_WITH_LIB([UMFPACK], [umfpack_di_symbolic], [umfpack])
AX_BLAS

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h sys/time.h cblas.h])
AC_HEADER_ASSERT

# Checks for typedefs, structures, and compiler characteristics.
# Checks for library functions.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_TYPE_UINT32_T
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRNLEN
AC_CHECK_FUNCS([gettimeofday strndup memset strcasecmp strchr strdup strerror])
AC_CHECK_HEADERS([malloc.h stddef.h])

# Setup the test suite
AC_CONFIG_TESTDIR([.],[.:tests])
AC_CONFIG_FILES([atlocal])

AC_MSG_CHECKING(for beer in the fridge)
  AC_MSG_RESULT(all systems go)

AC_OUTPUT


AC_MSG_RESULT([])
AC_MSG_RESULT([   zedhat Configuration Summary       ])
AC_MSG_RESULT([======================================])
AC_MSG_RESULT([ C Compiler: $CC])
AC_MSG_RESULT([     CFLAGS: $CFLAGS])
dnl AC_MSG_RESULT([     Shared Libraries: $enable_shared])
dnl AC_MSG_RESULT([     Static Libraries: $enable_static])
AC_MSG_RESULT([Solvers  -----------------------------])
AS_IF([test "x$ax_blas_ok" = xyes -a "x$BLAS_LIBS" != ""],
      AC_MSG_RESULT([       BLAS: $BLAS_LIBS]),
      AC_MSG_RESULT([       BLAS: $ax_blas_ok]))

AC_MSG_RESULT([])
