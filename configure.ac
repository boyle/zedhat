# Copyright 2016, Alistair Boyle, 3-clause BSD License
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([zedhat],[0.0.2],[alistair.js.boyle@gmail.com])

AC_COPYRIGHT([Copyright 2016, Alistair Boyle, 3-clause BSD License])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])

AC_CONFIG_SRCDIR([src/zedhat.c])
AC_CONFIG_HEADERS([src/config.h])
AC_USE_SYSTEM_EXTENSIONS
PKG_PROG_PKG_CONFIG

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
AM_CONDITIONAL([CC_IS_CLANG], [test x"${CC}" == x"clang"])
AM_PROG_CC_C_O
AC_PROG_CC_STDC
AC_C_INLINE
AC_TYPE_SIZE_T
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_TYPE_UINT32_T
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRNLEN
AC_CHECK_HEADER_STDBOOL # matrix.h
AC_CHECK_HEADERS([float.h stdlib.h string.h sys/time.h malloc.h stddef.h signal.h],,AC_MSG_ERROR([missing header]))
AC_CHECK_FUNCS([gettimeofday strndup memset strcasecmp strchr strdup strerror signal])
AC_HEADER_ASSERT
AC_CHECK_FUNCS([memmove strstr]) # cmocka.c
AC_FUNC_STRTOD # argv.c


# Checks for libraries.
PKG_CHECK_MODULES([BLAS],[blas],,true)
AX_BLAS(,AC_MSG_ERROR([missing blas]))
PKG_CHECK_MODULES([CBLAS],[cblas],,true)
AX_CBLAS(,AC_MSG_ERROR([missing cblas]))
#AC_CHECK_HEADERS([cblas.h],[cblas_h_ok=yes;break;])

PKG_CHECK_MODULES([LAPACK],[lapack],,true)
AX_LAPACK(,AC_MSG_ERROR([missing lapack]))
PKG_CHECK_MODULES([LAPACKE],[lapacke],,true)
AX_LAPACKE(,AC_MSG_ERROR([missing lapacke]))

PKG_CHECK_MODULES([UMFPACK],[umfpack],,true) # no pkg-config for umfpack on Ubuntu?
AC_CHECK_HEADERS([suitesparse/umfpack.h],[UMFPACK_CFLAGS="$UMFPACK_CFLAGS -I/usr/include/suitesparse";umfpack_h_ok=yes;break;])
AC_CHECK_HEADERS([umfpack.h],[umfpack_h_ok=yes;break;])
AS_IF([test "x$umfpack_h_ok" == "xyes"],,AC_MSG_ERROR([missing umfpack.h]))
save_LIBS="$LIBS"; LIBS="$UMFPACK_LIBS $LIBS"
AC_SEARCH_LIBS([umfpack_di_symbolic], [umfpack], , AC_MSG_ERROR([missing umfpack library]))
UMFPACK_LIBS="$LIBS"
LIBS="$save_LIBS"

PKG_CHECK_MODULES([CHOLMOD],[cholmod],,true) # no pkg-config for cholmod on Ubuntu?
AC_CHECK_HEADERS([suitesparse/cholmod.h],[CHOLMOD_CFLAGS="$CHOLMOD_CFLAGS -I/usr/include/suitesparse";cholmod_h_ok=yes;break;])
AC_CHECK_HEADERS([cholmod.h],[cholmod_h_ok=yes;break;])
AS_IF([test "x$cholmod_h_ok" == "xyes"],,AC_MSG_ERROR([missing cholmod.h]))
save_LIBS="$LIBS"; LIBS="$CHOLMOD_LIBS $LIBS"
AC_SEARCH_LIBS([cholmod_analyze], [cholmod], , AC_MSG_ERROR([missing cholmod library]))
CHOLMOD_LIBS="$LIBS"
LIBS="$save_LIBS"

AX_CHECK_ZLIB

dnl  dnl Definitions at top of <config.h>.
dnl  AH_TOP([
dnl  /* BEGIN_C_DECLS should be used at the beginning of your declarations,
dnl     so that C++ compilers don't mangle their names.  Use END_C_DECLS at
dnl     the end of C declarations. */
dnl  #undef BEGIN_C_DECLS
dnl  #undef END_C_DECLS
dnl  #ifdef __cplusplus
dnl  # define BEGIN_C_DECLS extern "C" {
dnl  # define END_C_DECLS }
dnl  #else
dnl  # define BEGIN_C_DECLS /* empty */
dnl  # define END_C_DECLS /* empty */
dnl  #endif
dnl
dnl  /* PARAMS is a macro used to wrap function prototypes, so that
dnl     compilers that don't understand ANSI C prototypes still work,
dnl     and ANSI C compilers can issue warnings about type mismatches. */
dnl  #undef PARAMS
dnl  #if defined __STDC__ || defined _AIX \
dnl          || (defined __mips && defined _SYSTYPE_SVR4) \
dnl          || defined WIN32 || defined __cplusplus
dnl  # define PARAMS(protos) protos
dnl  #else
dnl  # define PARAMS(protos) ()
dnl  #endif
dnl  ])
dnl  dnl USAGE:
dnl  dnl #ifndef FOO_H
dnl  dnl #define FOO_H 1
dnl  dnl /* The above macro definitions. */
dnl  dnl #include "..."
dnl  dnl BEGIN_C_DECLS
dnl  dnl int foo PARAMS((void));
dnl  dnl int hello PARAMS((void));
dnl  dnl END_C_DECLS
dnl  dnl #endif /* !FOO_H */

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
AC_MSG_RESULT([       LIBS: $LIBS])
dnl AC_MSG_RESULT([     Shared Libraries: $enable_shared])
dnl AC_MSG_RESULT([     Static Libraries: $enable_static])
AC_MSG_RESULT([Solvers  -----------------------------])
AS_IF([test "x$ax_blas_ok" = xyes -a "x$BLAS_LIBS" != ""],
      AC_MSG_RESULT([       BLAS: $BLAS_LIBS]),
      AC_MSG_RESULT([       BLAS: $ax_blas_ok]))
AS_IF([test "x$ax_cblas_ok" = xyes -a "x$CBLAS_LIBS" != ""],
      AC_MSG_RESULT([      CBLAS: $CBLAS_LIBS]),
      AC_MSG_RESULT([      CBLAS: $ax_cblas_ok]))
AS_IF([test "x$ax_lapack_ok" = xyes -a "x$LAPACK_LIBS" != "x"],
      AC_MSG_RESULT([     LAPACK: $LAPACK_LIBS]),
      AC_MSG_RESULT([     LAPACK: $ax_lapack_ok]))
AS_IF([test "x$ax_lapacke_ok" = xyes -a "x$LAPACKE_LIBS" != "x"],
      AC_MSG_RESULT([    LAPACKE: $LAPACKE_LIBS]),
      AC_MSG_RESULT([    LAPACKE: $ax_lapacke_ok]))
AS_IF([test "x$umfpack_h_ok" = xyes -a "x$UMFPACK_LIBS" != "x"],
      AC_MSG_RESULT([    UMFPACK: $UMFPACK_LIBS ; $UMFPACK_CFLAGS]),
      AC_MSG_RESULT([    UMFPACK: $umfpack_h_ok]))
AS_IF([test "x$cholmod_h_ok" = xyes -a "x$CHOLMOD_LIBS" != "x"],
      AC_MSG_RESULT([    CHOLMOD: $CHOLMOD_LIBS ; $CHOLMOD_CFLAGS]),
      AC_MSG_RESULT([    CHOLMOD: $cholmod_h_ok]))

AC_MSG_RESULT([])
