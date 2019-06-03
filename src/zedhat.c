/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdio.h>
#include <stdlib.h> /* malloc, free */
#include <string.h> /* memcpy */
#include <assert.h> /* assert */
#include <cblas.h>
#include <lapacke.h>

#include "argv.h"
#include "file.h"
#include "model.h"
#include "matrix.h"

int build_system_matrix(const model * m, matrix * A)
{
    if((m == NULL) || (A == NULL)) {
        return 0;
    }
    const mesh * fwd = &(m->fwd);
    const int dim = fwd->dim;
    const int ne = fwd->n_elems;
    const int se_n = calc_Se_n(dim); /* sparse matrix entries per mesh element */
    const int nn = fwd->n_nodes;
    size_t nnz = se_n * ne;
    if(!malloc_matrix_data(A, COO, nn, nn, nnz)) {
        return 0;
    }
    /* build shape matrices */
    int ret = calc_Se(fwd, A->x.sparse.ia, A->x.sparse.ja, A->x.sparse.a);
    if (ret != 0) {
        return 0;
    }
    assert(nnz == A->x.sparse.na);
    assert(nnz == A->x.sparse.nia);
    assert(nnz == A->x.sparse.nja);
    /* apply conductivity */
    for(int e = 0; e < ne; e++) {
        const double cond = m->params[e]; /* conductivity */
        for(int i = 0; i < se_n; i++) {
            A->x.sparse.a[e * se_n + i] *= cond;
        }
    }
    /* remove ground node */
    int gnd = calc_gnd(1, &nnz, A->x.sparse.ia, A->x.sparse.ja, A->x.sparse.a);
    if(gnd == 0) {
        return 0;
    }
    /* TODO fix redundancy */
    A->x.sparse.na = nnz;
    A->x.sparse.nia = nnz;
    A->x.sparse.nja = nnz;
    return 1;
}

int build_dense_matrix(const model * m, matrix * X)
{
    if((m == NULL) || (X == NULL)) {
        return 0;
    }
    const mesh * fwd = &(m->fwd);
    const int rows = fwd->n_nodes;
    if(!malloc_matrix_data(X, DENSE, rows, 1, rows)) {
        return 0;
    }
    for(int i = 0; i < rows; i++) {
        X->x.dense[i] = 0.0;
    }
    return 1;
}

int build_voltage_matrix(const model * m, matrix * X)
{
    return build_dense_matrix(m, X);
}

int build_current_matrix(const model * m, matrix * B)
{
    return build_dense_matrix(m, B);
}

int main(int argc, char ** argv)
{
    args_t args = {0};
    int ret = parse_argv(argc, argv, &args);
    if(ret != 0) {
        return 1;
    }
    if(args.mode == FORWARD_SOLVER) {
        ret = 1;
        printf("fwd: || A\\B - X ||_2 == || A X - B ||_2\n");
        model * M = malloc_model();
        double * AA = NULL;
        double * BB = NULL;
        matrix * A = malloc_matrix();
        matrix * B = malloc_matrix();
        matrix * X = malloc_matrix();
        int nret = readfile(args.file[0], M);
        if(!nret) {
            printf("error: read %s failed\n", args.file[0]);
            goto fwd_quit;
        }
        if(!malloc_matrix_name(A, "system", "A", "V/A")) {
            printf("error: matrix %s: out of memory\n", "A");
            goto fwd_quit;
        }
        if(!malloc_matrix_name(B, "currents", "B", "A")) {
            printf("error: matrix %s: out of memory\n", "B");
            goto fwd_quit;
        }
        if(!malloc_matrix_name(X, "voltages", "X", "V")) {
            printf("error: matrix %s: out of memory\n", "X");
            goto fwd_quit;
        }
        /* fill in the matrices and vectors */
        if (!build_system_matrix(M, A)) {
            printf("error: failed to build system matrix A\n");
            goto fwd_quit;
        }
        if (!build_current_matrix(M, B)) {
            printf("error: failed to build current matrix B\n");
            goto fwd_quit;
        }
        if (!build_voltage_matrix(M, X)) {
            printf("error: failed to build voltage matrix X\n");
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
            goto fwd_quit;
        }
        /* Find: || A X - B ||_2 */
        /* Compute  C = a A B + c C */
        memcpy(BB, B->x.dense, Bmn * sizeof(double));
        printf_matrix(A);
        printf_matrix(X);
        printf_matrix(B);
        cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
                     A->m, X->n, A->n, /* m, n, k, */
                     +1.0, A->x.dense, A->m, X->x.dense, X->m, /* a, A, m, B, n, */
                     -1.0, BB, B->m); /* c, C, k */
        /* Compute  || B ||_2 */
        const double err_mul = cblas_dnrm2( Bmn, BB, 1); /* n, X, incX */
        if (err_mul > args.tol) {
            printf("fail: || A X - B ||_2 = %g\n", err_mul);
            ret = 1;
            goto fwd_quit;
        }
        /* Find BB = A\B; || BB - X ||_2 */
        memcpy(AA, A->x.dense, Amn * sizeof(double));
        memcpy(BB, B->x.dense, Bmn * sizeof(double));
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
        for(int i = 0; i < Xmn; i++) {
            BB[i] -= X->x.dense[i];
        }
        const double err_div = cblas_dnrm2( Xmn, BB, 1); /* n, X, incX */
        if (err_div > args.tol) {
            printf("fail: || A\\B - X ||_2 = %g\n", err_div);
            ret = 1;
            goto fwd_quit;
        }
        printf("pass\n");
        ret = 0;
fwd_quit:
        free_model(M);
        free_matrix(A);
        free_matrix(B);
        free_matrix(X);
        free(AA);
        free(BB);
    }
    return ret;
}
