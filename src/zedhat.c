/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdio.h>
#include <stdlib.h> /* malloc, free */
#include <string.h> /* memcpy */
#include <assert.h> /* assert */
#include <math.h> /* DBL_EPSILON */
#include <cholmod.h>
//#include <cblas.h>
//#include <lapacke.h>

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
    if(!malloc_matrix_data(A, COO_SYMMETRIC, nn, nn, nnz)) {
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
        printf_matrix(A);
        printf_matrix(X);
        printf_matrix(B);
        printf("CHOLMOD START\n");
        /* CHOLMOD sparse symmetric decomposition */
        cholmod_dense * x = NULL, *b = NULL, *r = NULL;
        cholmod_factor * L = NULL;
        cholmod_common c;
        cholmod_start (&c);                /* start CHOLMOD */
        cholmod_sparse * As = NULL;
        { /* COO to CSC */
            cholmod_triplet * T = cholmod_allocate_triplet(A->m, A->n, A->x.sparse.na, +1, CHOLMOD_REAL, &c);
            assert(T != NULL);
            { /* copy contents of A */
                T->nnz = A->x.sparse.na;
                double * xx = T->x;
                int * ii = T->i;
                int * jj = T->j;
                for(int i = 0; i < T->nnz; i++) {
                    ii[i] = A->x.sparse.ia[i];
                    jj[i] = A->x.sparse.ja[i];
                    xx[i] = A->x.sparse.a[i];
                }
            }
            A = free_matrix(A);
            cholmod_print_triplet (T, "T", &c);         /* print the matrix */
            As = cholmod_triplet_to_sparse(T, T->nnz, &c);
            cholmod_print_sparse (As, "A", &c);         /* print the matrix */
            cholmod_free_triplet(&T, &c);
        }
        assert(As != NULL);
        cholmod_check_sparse (As, &c);
        /* END MODIFIED: now continue with the usual operations */
        assert(As->stype != 0);
        /* create stimulus */
        b = cholmod_zeros(As->nrow, 1, As->xtype, &c);   /* b = zeros(n,1) */
        double * bx = b->x;
        bx[0] = +1;
        bx[50] = -1;
        /* end stimulus */
        L = cholmod_analyze (As, &c);           /* analyze */
        cholmod_factorize (As, L, &c);          /* factorize */
        x = cholmod_solve (CHOLMOD_A, L, b, &c);       /* solve Ax=b */
        r = cholmod_copy_dense (b, &c);            /* r = b */
        double one [2] = {1, 0}, m1 [2] = {-1, 0} ;     /* basic (real-valued) scalars */
        cholmod_sdmult (As, 0, m1, one, x, r, &c);      /* r = r-Ax */
        printf ("norm(b-Ax) %8.1e\n", cholmod_norm_dense (r, 0, &c));        /* print norm(r) */
        cholmod_free_factor (&L, &c);          /* free matrices */
        cholmod_free_sparse (&As, &c);
        cholmod_free_dense (&r, &c);
        cholmod_free_dense (&x, &c);
        cholmod_free_dense (&b, &c);
        cholmod_finish (&c);               /* finish CHOLMOD */
        printf("CHOLMOD DONE\n");
        // /* Find: || A X - B ||_2 */
        // /* Compute  C = a A B + c C */
        // cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
        //              A->m, X->n, A->n, /* m, n, k, */
        //              +1.0, A->x.dense, A->m, X->x.dense, X->m, /* a, A, m, B, n, */
        //              -1.0, BB, B->m); /* c, C, k */
        // /* Compute  || B ||_2 */
        // const double err_mul = cblas_dnrm2( Bmn, BB, 1); /* n, X, incX */
        // if (err_mul > args.tol) {
        //     printf("fail: || A X - B ||_2 = %g\n", err_mul);
        //     ret = 1;
        //     goto fwd_quit;
        // }
        // /* Find BB = A\B; || BB - X ||_2 */
        // memcpy(AA, A->x.dense, Amn * sizeof(double));
        // memcpy(BB, B->x.dense, Bmn * sizeof(double));
        // // memset(&BB[Bmn], 0, (Xmn-Bmn)*sizeof(double));
        // lapack_int err = LAPACKE_dgels( LAPACK_COL_MAJOR, 'N', /* row/col major, trans, */
        //                                 A->m, A->n, B->n, /* m, n, n_rhs, */
        //                                 AA, A->m, BB, B->m ); /* A, m, B, n */
        // if (err) {
        //     printf("fail: DGELS (X=A\\B) info=%d\n", err);
        //     ret = 1;
        //     goto fwd_quit;
        // }
        // /* Compute  || BB - X ||_2 */
        // for(int i = 0; i < Xmn; i++) {
        //     BB[i] -= X->x.dense[i];
        // }
        // const double err_div = cblas_dnrm2( Xmn, BB, 1); /* n, X, incX */
        // if (err_div > args.tol) {
        //     printf("fail: || A\\B - X ||_2 = %g\n", err_div);
        //     ret = 1;
        //     goto fwd_quit;
        // }
        // printf("pass\n");
        ret = 0;
fwd_quit:
        free_model(M);
        free_matrix(A);
        free_matrix(B);
        free_matrix(X);
    }
    return ret;
}
