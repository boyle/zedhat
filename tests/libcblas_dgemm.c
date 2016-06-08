/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include <stdio.h>
#include <stdlib.h> /* for exit */
#include <math.h>  /* for fabs */
#include <cblas.h> /* for cblas_dgemm */

int
main (void)
{
    double const a[] = { 0.11, 0.12, 0.13,
                         0.21, 0.22, 0.23
                       };
    double const b[] = { 1011, 1012,
                         1021, 1022,
                         1031, 1032
                       };
    double c[] = { 0.00, 0.00,
                   0.00, 0.00
                 };
    /* Compute C = a A B + c C */
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
                 3, 2, 3, /* m, n, k, */
                 1.0, a, 3, b, 2, /* a, A, m, B, n, */
                 0.0, c, 2); /* c, C, k */
    /* A is m x k, B is k x n, C is m x n */
    printf ("[ %g, %g\n", c[0], c[1]);
    printf ("  %g, %g ]\n", c[2], c[3]);
    double const expect[] = { 367.76, 368.12,
                              674.06, 674.72
                            };
    int i;
    for(i = 0; i < 4; i++) {
        if ( fabs(expect[i] - c[i]) > 1e-12 ) {
            printf("FAIL\n");
            return 1;
        }
    }
    printf("PASS\n");
    exit(EXIT_SUCCESS);
}
