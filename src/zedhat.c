/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include <stdio.h>
#include <stdlib.h> /* malloc, free */
#include <string.h> /* memcpy */
#include <cblas.h>
#include <lapacke.h>
#include "config.h"
#include "argv.h"
#include "matrix.h"

matrix * get_fwd_matrix(const args_t args, const char * name, const char * symbol, const char * units)
{
    matrix * A = malloc_matrix();
    malloc_matrix_name(A, name, symbol, units);
    if (A == NULL) {
        printf("fwd %s: ret = malloc failed (%s)\n", name, args.file[0]);
        return NULL;
    }
    int ret = 1; // FIXME matrix_load(args.file[0], A);
    if(ret) {
        printf("fwd %s: ret = %d (%s)\n", name, ret, args.file[0]);
        free_matrix(A);
        return NULL;
    }
    else {
        return A;
    }
}

#ifdef DEBUG
void matrix_printf(matrix * M, double * D)
{
    int i, j;
    printf("--- %10s (%2s) ---\n", M->name, M->symbol);
    for(i = 0; i < M->m; i++) {
        for(j = 0; j < M->n; j++) {
            printf("\t%g", D[i + (j * M->m)]);
        }
        printf("\n");
    }
    printf("----%10s%6s----\n", "----------", "------");
}
#endif

int main(int argc, char ** argv)
{
    args_t args = {0};
    int ret = parse_argv(argc, argv, &args);
    if(ret != 0) {
        return 1;
    }
    if(args.mode == FORWARD_SOLVER) {
        printf("fwd: || A\\B - X ||_2 == || A X - B ||_2\n");
        double * AA = NULL;
        double * BB = NULL;
        /* load: A X = B */
        matrix * A = get_fwd_matrix(args, "A", NULL, NULL);
        matrix * B = get_fwd_matrix(args, "B", NULL, NULL);
        matrix * X = get_fwd_matrix(args, "X", NULL, NULL);
        if (!A || !B || !X) {
            ret = 1;
            goto fwd_quit;
        }
        const size_t Amn = A->m * A->n;
        const size_t Bmn = B->m * B->n;
        const size_t Xmn = X->m * X->n;
        const size_t BBn = (Bmn > Xmn ? Bmn : Xmn);
        AA = (double *) malloc(Amn * sizeof(double));
        BB = (double *) malloc(BBn * sizeof(double));
        if (!AA || !BB) {
            printf("fail: out of memory\n");
            ret = 1;
            goto fwd_quit;
        }
        /* Find: || A X - B ||_2 */
        /* Compute  C = a A B + c C */
        memcpy(BB, B->dense, Bmn * sizeof(double));
        cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                     A->m, X->n, A->n, /* m, n, k, */
                     +1.0, A->dense, A->m, X->dense, X->m, /* a, A, m, B, n, */
                     -1.0, BB, B->m); /* c, C, k */
        /* Compute  || B ||_2 */
        const double err_mul = cblas_dnrm2( Bmn, BB, 1); /* n, X, incX */
        if (err_mul > args.tol) {
            printf("fail: || A X - B ||_2 = %g\n", err_mul);
            ret = 1;
            goto fwd_quit;
        }
        /* Find BB = A\B; || BB - X ||_2 */
        memcpy(AA, A->dense, Amn * sizeof(double));
        memcpy(BB, B->dense, Bmn * sizeof(double));
        // memset(&BB[Bmn], 0, (Xmn-Bmn)*sizeof(double));
        lapack_int err = LAPACKE_dgels( LAPACK_COL_MAJOR, 'N', /* row/col major, trans, */
                                        A->m, A->n, B->n, /* m, n, n_rhs, */
                                        AA, A->m, BB, B->m ); /* A, m, B, n */
        if (err) {
            printf("fail: DGELS (X=A\\B) info=%d\n", err);
            ret = 1;
            goto fwd_quit;
        }
        /* Compute  || BB - X ||_2 */
        int i;
        for(i = 0; i < Xmn; i++) {
            BB[i] -= X->dense[i];
        }
        const double err_div = cblas_dnrm2( Xmn, BB, 1); /* n, X, incX */
        if (err_div > args.tol) {
            printf("fail: || A\\B - X ||_2 = %g\n", err_div);
            ret = 1;
            goto fwd_quit;
        }
        printf("pass\n");
fwd_quit:
        free_matrix(A);
        free_matrix(B);
        free_matrix(X);
        free(AA);
        free(BB);
    }
    return ret;
}
