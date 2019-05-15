# ===========================================================================
#         https://www.gnu.org/software/autoconf-archive/ax_blas.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CBLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a library that implements the CBLAS linear-algebra
#   interface (see http://www.netlib.org/blas/). On success, it sets the
#   CBLAS_LIBS output variable to hold the requisite library linkages.
#
#   To link with CBLAS, you should link with:
#
#     $CBLAS_LIBS $LIBS $FLIBS
#
#   in that order. FLIBS is the output variable of the
#   AC_F77_LIBRARY_LDFLAGS macro (called if necessary by AX_CBLAS), and is
#   sometimes necessary in order to link with F77 libraries. Users will also
#   need to use AC_F77_DUMMY_MAIN (see the autoconf manual), for the same
#   reason.
#
#   Many libraries are searched for, from ATLAS to CXML to ESSL. The user
#   may also use --with-blas=<lib> in order to use some specific CBLAS
#   library <lib>. In order to link successfully, however, be aware that you
#   will probably need to use the same Fortran compiler (which can be set
#   via the F77 env. var.) as was used to compile the CBLAS library.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a CBLAS library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_CBLAS.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
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
#   with this program. If not, see <https://www.gnu.org/licenses/>.
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

#serial 16

AU_ALIAS([ACX_CBLAS], [AX_CBLAS])
AC_DEFUN([AX_CBLAS], [
AC_PREREQ([2.55])
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
AC_REQUIRE([AC_CANONICAL_HOST])
ax_cblas_ok=no

AC_ARG_WITH(blas,
	[AS_HELP_STRING([--with-blas=<lib>], [use CBLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) ax_cblas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) CBLAS_LIBS="$with_blas" ;;
	*) CBLAS_LIBS="-l$with_blas" ;;
esac

sgemm=cblas_sgemm
dgemm=cblas_dgemm

ax_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check CBLAS_LIBS environment variable
if test $ax_cblas_ok = no; then
if test "x$CBLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$CBLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgemm in $CBLAS_LIBS])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$sgemm])], [ax_cblas_ok=yes], [CBLAS_LIBS=""])
	AC_MSG_RESULT($ax_cblas_ok)
	LIBS="$save_LIBS"
fi
fi

# CBLAS linked to by default?  (happens on some supercomputers)
if test $ax_cblas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_MSG_CHECKING([if $sgemm is being linked in already])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$sgemm])], [ax_cblas_ok=yes])
	AC_MSG_RESULT($ax_cblas_ok)
	LIBS="$save_LIBS"
fi

# CBLAS in OpenCBLAS library? (http://xianyi.github.com/OpenCBLAS/)
if test $ax_cblas_ok = no; then
	AC_CHECK_LIB(openblas, $sgemm, [ax_cblas_ok=yes
			                CBLAS_LIBS="-lopenblas"])
fi

# CBLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $ax_cblas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[ax_cblas_ok=yes
			 CBLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# CBLAS in PhiPACK libraries? (requires generic CBLAS lib, too)
if test $ax_cblas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[ax_cblas_ok=yes; CBLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# CBLAS in Intel MKL library?
if test $ax_cblas_ok = no; then
	# MKL for gfortran
	if test x"$ac_cv_fc_compiler_gnu" = xyes; then
		# 64 bit
		if test $host_cpu = x86_64; then
			AC_CHECK_LIB(mkl_gf_lp64, $sgemm,
			[ax_cblas_ok=yes;CBLAS_LIBS="-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread"],,
			[-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread])
		# 32 bit
		elif test $host_cpu = i686; then
			AC_CHECK_LIB(mkl_gf, $sgemm,
				[ax_cblas_ok=yes;CBLAS_LIBS="-lmkl_gf -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_gf -lmkl_sequential -lmkl_core -lpthread])
		fi
	# MKL for other compilers (Intel, PGI, ...?)
	else
		# 64-bit
		if test $host_cpu = x86_64; then
			AC_CHECK_LIB(mkl_intel_lp64, $sgemm,
				[ax_cblas_ok=yes;CBLAS_LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread])
		# 32-bit
		elif test $host_cpu = i686; then
			AC_CHECK_LIB(mkl_intel, $sgemm,
				[ax_cblas_ok=yes;CBLAS_LIBS="-lmkl_intel -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_intel -lmkl_sequential -lmkl_core -lpthread])
		fi
	fi
fi
# Old versions of MKL
if test $ax_cblas_ok = no; then
	AC_CHECK_LIB(mkl, $sgemm, [ax_cblas_ok=yes;CBLAS_LIBS="-lmkl -lguide -lpthread"],,[-lguide -lpthread])
fi

# CBLAS in Apple vecLib library?
if test $ax_cblas_ok = no; then
	save_LIBS="$LIBS"; LIBS="-framework vecLib $LIBS"
	AC_MSG_CHECKING([for $sgemm in -framework vecLib])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$sgemm])], [ax_cblas_ok=yes;CBLAS_LIBS="-framework vecLib"])
	AC_MSG_RESULT($ax_cblas_ok)
	LIBS="$save_LIBS"
fi

# CBLAS in Alpha CXML library?
if test $ax_cblas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [ax_cblas_ok=yes;CBLAS_LIBS="-lcxml"])
fi

# CBLAS in Alpha DXML library? (now called CXML, see above)
if test $ax_cblas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [ax_cblas_ok=yes;CBLAS_LIBS="-ldxml"])
fi

# CBLAS in Sun Performance library?
if test $ax_cblas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
				[CBLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 ax_cblas_ok=yes],[],[-lsunmath])])
	fi
fi

# CBLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $ax_cblas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [ax_cblas_ok=yes; CBLAS_LIBS="-lscs"])
fi

# CBLAS in SGIMATH library?
if test $ax_cblas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [ax_cblas_ok=yes; CBLAS_LIBS="-lcomplib.sgimath"])
fi

# CBLAS in IBM ESSL library? (requires generic CBLAS lib, too)
if test $ax_cblas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[ax_cblas_ok=yes; CBLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic CBLAS library?
if test $ax_cblas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm, [ax_cblas_ok=yes; CBLAS_LIBS="-lcblas"])
fi

AC_SUBST(CBLAS_LIBS)

LIBS="$ax_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_cblas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_CBLAS,1,[Define if you have a CBLAS library.]),[$1])
        :
else
        ax_cblas_ok=no
        $2
fi
])dnl AX_CBLAS
