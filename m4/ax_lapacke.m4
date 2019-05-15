# ===========================================================================
#         http://www.gnu.org/software/autoconf-archive/ax_lapacke.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a library that implements the LAPACK linear-algebra
#   interface (see http://www.netlib.org/lapacke/). On success, it sets the
#   LAPACKE_LIBS output variable to hold the requisite library linkages.
#
#   To link with LAPACK, you should link with:
#
#     $LAPACKE_LIBS $BLAS_LIBS $LIBS $FLIBS
#
#   in that order. BLAS_LIBS is the output variable of the AX_BLAS macro,
#   called automatically. FLIBS is the output variable of the
#   AC_F77_LIBRARY_LDFLAGS macro (called if necessary by AX_BLAS), and is
#   sometimes necessary in order to link with F77 libraries. Users will also
#   need to use AC_F77_DUMMY_MAIN (see the autoconf manual), for the same
#   reason.
#
#   The user may also use --with-lapacke=<lib> in order to use some specific
#   LAPACK library <lib>. In order to link successfully, however, be aware
#   that you will probably need to use the same Fortran compiler (which can
#   be set via the F77 env. var.) as was used to compile the LAPACK and BLAS
#   libraries.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a LAPACK library
#   is found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it
#   is not found. If ACTION-IF-FOUND is not specified, the default action
#   will define HAVE_LAPACK.
#
# LICENSE
#
#   Copyright (c) 2009 Steven G. Johnson <stevenj@alum.mit.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 7

AU_ALIAS([ACX_LAPACKE], [AX_LAPACKE])
AC_DEFUN([AX_LAPACKE], [
AC_REQUIRE([AX_LAPACK])
AC_REQUIRE([AX_CBLAS])
AC_REQUIRE([AX_BLAS])
ax_lapacke_ok=no

AC_ARG_WITH(lapacke,
        [AS_HELP_STRING([--with-lapacke=<lib>], [use LAPACKE library <lib>])])
case $with_lapacke in
        yes | "") ;;
        no) ax_lapacke_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACKE_LIBS="$with_lapacke" ;;
        *) LAPACKE_LIBS="-l$with_lapacke" ;;
esac

LAPACKE_dgels="LAPACKE_dgels"

# We cannot use LAPACKE if LAPACK is not found
if test "x$ax_lapack_ok" != xyes; then
        ax_lapacke_ok=noblas
        LAPACKE_LIBS=""
fi

# First, check LAPACKE_LIBS environment variable
if test "x$LAPACKE_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACKE_LIBS $LAPACK_LIBS $CBLAS_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $LAPACKE_dgels in $LAPACKE_LIBS])
        AC_TRY_LINK_FUNC($LAPACKE_dgels, [ax_lapacke_ok=yes], [LAPACKE_LIBS=""])
        AC_MSG_RESULT($ax_lapacke_ok)
        LIBS="$save_LIBS"
        if test $ax_lapacke_ok = no; then
                LAPACKE_LIBS=""
        fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $ax_lapacke_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $LAPACK_LIBS $CBLAS_LIBS $BLAS_LIBS $FLIBS"
        AC_CHECK_FUNC($LAPACKE_dgels, [ax_lapacke_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapacke in lapacke lapacke_rs6k; do
        if test $ax_lapacke_ok = no; then
                save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $CBLAS_LIBS $BLAS_LIBS $LIBS"
                AC_CHECK_LIB($lapacke, $LAPACKE_dgels,
                    [ax_lapacke_ok=yes; LAPACKE_LIBS="-l$lapacke"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACKE_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_lapacke_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACKE,1,[Define if you have LAPACKE library.]),[$1])
        :
else
        ax_lapacke_ok=no
        $2
fi
])dnl AX_LAPACK
