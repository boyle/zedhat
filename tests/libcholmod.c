/* from SuiteSparse cholmod
 * https://github.com/PetterS/SuiteSparse/blob/master/CHOLMOD/Demo/cholmod_simple.c
 * modified to contain the matrix in-line
 * */

/* ========================================================================== */
/* === Demo/cholmod_simple ================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Demo Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/Demo Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Read in a real symmetric or complex Hermitian matrix from stdin in
 * MatrixMarket format, solve Ax=b where b=[1 1 ... 1]', and print the residual.
 * Usage: cholmod_simple < matrixfile
 */

#include <cholmod.h>

#include <stdarg.h> /* cmocka.h */
#include <stddef.h> /* cmocka.h */
#include <setjmp.h> /* cmocka.h */
#include <stdio.h> /* cmocka.h */
#include <string.h> /* cmocka.h */
#include "cmocka.h"

#include <float.h> /* DBL_EPSILON */
#include <math.h> /* fabs, pow */

int main (void)
{
    cholmod_dense * x, *b, *r ;
    cholmod_factor * L ;
    double one [2] = {1, 0}, m1 [2] = {-1, 0} ;     /* basic scalars */
    cholmod_common c ;
    cholmod_start (&c) ;                /* start CHOLMOD */
    /* A = cholmod_read_sparse (stdin, &c) ; */     /* read in a matrix */
    /* MODIFIED: have a built in matrix as a COO object, then convert to CSR */
    /* MATLAB: A = [ 1 2 3 0; 2 0 0 6; 3 0 0 7; 0 6 7 9 ]; b = ones(4,1); A\b
     * ans = [ -0.2500, -3.1250, 2.5000, 0.2500 ]' */
    double expect[4] = {-0.2500, -3.1250, 2.5000, 0.2500 };
    double Adx [4][4] = { { 1, 2, 3, 0 },
        { 2, 0, 0, 6 },
        { 3, 0, 0, 7 },
        { 0, 6, 7, 9 }
    };
    cholmod_dense Ad = { 4, 4, 16, 4, &Adx, NULL, CHOLMOD_REAL, CHOLMOD_DOUBLE };
    cholmod_sparse * Atmp = cholmod_dense_to_sparse(&Ad, 1, &c);
    if (Atmp == NULL) {
        cholmod_free_sparse (&Atmp, &c) ;
        cholmod_finish (&c) ;
        printf("FAIL\n");
        return (1) ;
    }
    cholmod_sparse * A = cholmod_copy(Atmp, 1, 1, &c);
    cholmod_free_sparse (&Atmp, &c) ;
    cholmod_check_sparse (A, &c);
    /* END MODIFIED: now continue with the usual operations */
    cholmod_print_sparse (A, "A", &c) ;         /* print the matrix */
    if (A == NULL || A->stype == 0) {       /* A must be symmetric */
        cholmod_free_sparse (&A, &c) ;
        cholmod_finish (&c) ;
        printf("FAIL\n");
        return (1) ;
    }
    b = cholmod_ones (A->nrow, 1, A->xtype, &c) ;   /* b = ones(n,1) */
    L = cholmod_analyze (A, &c) ;           /* analyze */
    cholmod_factorize (A, L, &c) ;          /* factorize */
    x = cholmod_solve (CHOLMOD_A, L, b, &c) ;       /* solve Ax=b */
    r = cholmod_copy_dense (b, &c) ;            /* r = b */
    cholmod_sdmult (A, 0, m1, one, x, r, &c) ;      /* r = r-Ax */
    printf ("norm(b-Ax) %8.1e\n",
            cholmod_norm_dense (r, 0, &c)) ;        /* print norm(r) */
    double * tmp = x->x;
    int i;
    for (i = 0; i < 4; i++) {
        printf("expect[%d] = % 10.4f, x[%d] = % 10.4f, delta = %+g\n", i, expect[i], i, tmp[i], expect[i] - tmp[i]);
        if(fabs(expect[i] - tmp[i]) > 100 * DBL_EPSILON) {
            printf("FAIL\n");
            return 1;
        }
    }
    cholmod_free_factor (&L, &c) ;          /* free matrices */
    cholmod_free_sparse (&A, &c) ;
    cholmod_free_dense (&r, &c) ;
    cholmod_free_dense (&x, &c) ;
    cholmod_free_dense (&b, &c) ;
    cholmod_finish (&c) ;               /* finish CHOLMOD */
    printf("PASS\n");
    return (0) ;
}
