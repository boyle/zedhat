/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
/* Based on illustrative sample from UMFPACK's QuickStart.pdf
 * http://userpages.umbc.edu/~rostamia/2003-09-math625/umfpack.html
 * modified 2016, A. Boyle */

#include <stdio.h>
#include <math.h> /* for fabs */
#include <umfpack.h>

void check_errcode(int status)
{
    switch (status) {
    case UMFPACK_ERROR_out_of_memory:
        fprintf(stderr, "out of memory!\n");
        break;
    case UMFPACK_WARNING_singular_matrix:
        fprintf(stderr, "matrix is singular!\n");
        break;
    case UMFPACK_OK:
        break;
    default:
        fprintf(stderr, "UMFPACK error code %d\n", status);
    }
    if(status != UMFPACK_OK) {
        exit(EXIT_FAILURE);
    }
}

int main(void)
{
    /* A = 2   3   0   0   0
           3   0   4   0   6
           0  -1  -3   2   0
           0   0   1   0   0
           0   4   2   0   1 */
    /* as row, column, data, scanned column-wise */
    int Ti[] =    { 0, 1, 0,  2, 4, 1,  2, 3, 4, 2, 1, 4 }; /* row index */
    int Tj[] =    { 0, 0, 1,  1, 1, 2,  2, 2, 2, 3, 4, 4 }; /* col index */
    double Tx[] = { 2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1 }; /* data */
    int Cp[6];
    int Ci[12];
    double Cx[12];
    /* as Compressed Sparse Column (CSC) */
    int Ap[] = { 0, 2, 5, 9, 10, 12 }; /* column non-zero count, skip first idx */
    int Ai[] =    { 0, 1, 0,  2, 4, 1,  2, 3, 4, 2, 1, 4 }; /* row indices */
    double Ax[] = { 2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1 }; /* data */
    double b[] = { 8, 45, -3, 3, 19 }; /* ^T */
    void * Symbolic, *Numeric;
    double x[5];
    int status;
    int i;
    status = umfpack_di_triplet_to_col(5, 5, 12, Ti, Tj, Tx, Cp, Ci, Cx, NULL); /* convert T->A */
    check_errcode(status);
    for(i = 0; i < 6; i++) {
        if(Cp[i] != Ap[i]) {
            printf("FAIL Ap\n");
            exit(EXIT_FAILURE);
        }
    }
    for(i = 0; i < 12; i++) {
        if(Ci[i] != Ai[i]) {
            printf("FAIL Ai\n");
            exit(EXIT_FAILURE);
        }
        if(fabs(Cx[i] - Ax[i]) > 2e-15) {
            printf("FAIL Ax\n");
            exit(EXIT_FAILURE);
        }
    }
    status = umfpack_di_symbolic(5, 5, Ap, Ai, Ax, &Symbolic, NULL, NULL); /* symbolic analysis */
    check_errcode(status);
    status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL); /* LU factorization */
    umfpack_di_free_symbolic(&Symbolic);
    check_errcode(status);
    status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL); /* solve system */
    umfpack_di_free_numeric(&Numeric);
    check_errcode(status);
    /* matlab: expect = A \ b */
    double const expect[] = {1, 2, 3, 4, 5}; /* ^T */
    for (i = 0; i < 5; i++) {
        printf(" %g\n", x[i]);
    }
    for(i = 0; i < 5; i++) {
        if(fabs(expect[i] - x[i]) > 2e-15) {
            printf("FAIL\n");
            exit(EXIT_FAILURE);
        }
    }
    printf("PASS\n");
    exit(EXIT_SUCCESS);
}
