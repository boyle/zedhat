/* Copyright 2019, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <assert.h> /* assert */
#include <cholmod.h>
#include "matrix.h" /* coo_to_csc */

#include "fwd.h"

#define SUCCESS 1
#define FAILURE 0

static int malloc_system_matrix(const model * mdl, matrix ** A)
{
    *A = malloc_matrix();
    if(!malloc_matrix_name(*A, "system", "A", "V/A")) {
        printf("error: matrix %s: out of memory\n", "A");
        return FAILURE;
    }
    return SUCCESS;
}

static int build_system_matrix(const model * m, int gnd, int frame_idx, matrix * A, matrix * PMAP)
{
    assert(m != NULL);
    assert(A != NULL);
    assert(frame_idx >= 0);
    assert(frame_idx < m->n_params[0]);
    assert(gnd > 0);
    assert(gnd < m->fwd.n_nodes);
    /* build shape matrices */
    /* TODO we could reuse a copy of the expensive to compute shape elements at
     * each iteration, intead of recomputing from scratch here */
    const int rows = calc_sys_size(m); /* = cols */
    size_t nnz = calc_sys_nnz(m);
    if(!malloc_matrix_data(A, COO_SYMMETRIC, rows, rows, nnz)) {
        return FAILURE;
    }
    int ret = calc_sys_elem(&(m->fwd), A->x.sparse.ia, A->x.sparse.ja, A->x.sparse.a);
    if (ret != 0) {
        return FAILURE;
    }
    /* apply conductivity */
    const int ne = m->fwd.n_elems;
    const int dim = m->fwd.dim;
    const int se_n = calc_sys_elem_n(dim); /* sparse matrix entries per mesh element */
    if(PMAP == NULL) {
        assert(ne == m->n_params[0]);
        for(int e = 0; e < ne; e++) {
            const double cond = m->params[e * m->n_params[1] + frame_idx]; /* conductivity */
            for(int i = 0; i < se_n; i++) {
                A->x.sparse.a[e * se_n + i] *= cond;
            }
        }
    }
    else { /* map parameters to elements */
        assert(PMAP->n == m->n_params[0]);
        assert(PMAP->m == ne);
        assert(PMAP->type == CSR);
        for(int e = 0; e < ne; e++) {
            double cond_frac = 0.0;
            for(int j = PMAP->x.sparse.ia[e]; j < PMAP->x.sparse.ia[e + 1]; j++) {
                const int param = PMAP->x.sparse.ja[j];
                const double element_fraction = PMAP->x.sparse.a[j];
                const double conductivity = m->params[frame_idx + param * m->n_params[1]];
                cond_frac += conductivity * element_fraction;
            }
            // printf("[e=%d] σ=%g\n", e, cond_frac);
            for(int i = 0; i < se_n; i++) {
                A->x.sparse.a[e * se_n + i] *= cond_frac;
            }
        }
    }
    /* build CEM */
    int Se_idx = (se_n * ne); /* length of regular FEM entries in COO matrix */
    const int nn = m->fwd.n_nodes;
    for(int i = 0; i < m->n_elec; i++) {
        const int row = m->elec_to_sys[i];
        if(row < nn) {
            continue; /* is PEM */
        }
        const double zc = m->zc[i]; /* Ω/m² */
        const int bc = i + 1; /* electrode# = bc# = i+1 */
        Se_idx += calc_sys_cem(&(m->fwd), bc, zc, row, Se_idx, A->x.sparse.ia, A->x.sparse.ja, A->x.sparse.a);
    }
    assert(nnz == Se_idx);
    /* remove ground node */
    calc_sys_gnd(gnd, nnz, A->x.sparse.ia, A->x.sparse.ja, A->x.sparse.a);
    return SUCCESS;
}

static cholmod_sparse * copy_csc_to_cholmod_sparse(matrix * A, cholmod_common * c)
{
    assert(A != NULL);
    assert(A->type == CSC);
    const int sorted = 1;
    const int packed = 1;
    const int stype = +1;
    const int xtype = CHOLMOD_REAL;
    const int nnz = A->x.sparse.na;
    const int rows = A->m;
    const int cols = A->n;
    cholmod_sparse * As = cholmod_allocate_sparse(rows, cols, nnz, sorted, packed, stype, xtype, c);
    assert(As != NULL);
    { /* copy contents from A */
        double * xx = As->x;
        int * ii = As->i;
        int * jj = As->p;
        for(int i = 0; i < nnz; i++) {
            ii[i] = A->x.sparse.ia[i];
            xx[i] = A->x.sparse.a[i];
        }
        for(int i = 0; i < cols + 1; i++) {
            jj[i] = A->x.sparse.ja[i];
        }
    }
    return As;
}

int build_pmap(model * mdl, matrix ** PMAP)
{
    assert(*PMAP == NULL);
    if(mdl->fwd.n_pmap == 0) {
        return SUCCESS;
    }
    *PMAP = malloc_matrix();
    if(!malloc_matrix_name(*PMAP, "params_to_elems", "PMAP", "")) {
        return FAILURE;
    }
    const int rows = mdl->fwd.n_elems;
    const int cols = mdl->n_params[0];
    const int nnz = mdl->fwd.n_pmap;
    if(!malloc_matrix_data(*PMAP, COO, rows, cols, nnz)) {
        return FAILURE;
    }
    for(int i = 0; i < nnz; i++) {
        (*PMAP)->x.sparse.a[i] = mdl->fwd.pmap_frac[i]; /* TODO faster w/ memcpy? */
        (*PMAP)->x.sparse.ia[i] = mdl->fwd.pmap_elem[i] - 1; /* TODO faster w/ memcpy? */
        (*PMAP)->x.sparse.ja[i] = mdl->fwd.pmap_param[i] - 1; /* TODO faster w/ memcpy? */
    }
    if(!coo_to_csr(*PMAP)) {
        return FAILURE;
    }
    return SUCCESS;
}

int fwd_solve(model * mdl, double * meas)
{
    assert(mdl != NULL);
    assert(meas != NULL);
    if(!check_model(mdl)) {
        return FAILURE;
    }
    const int rows = calc_sys_size(mdl);
    int ret = FAILURE;
    /* CHOLMOD sparse symmetric decomposition */
    cholmod_dense * x = NULL, *b = NULL;
    cholmod_factor * L = NULL;
    cholmod_dense * E = NULL, *Y = NULL;
    cholmod_sparse * xset = NULL;
    cholmod_common c;
    cholmod_sparse * As = NULL;
    matrix * A = NULL;
    matrix * PMAP = NULL;
    cholmod_start (&c);                /* start CHOLMOD */
    if(!malloc_system_matrix(mdl, &A)) {
        printf("error: matrix %s: out of memory\n", "A");
        goto return_result;
    }
    if(!build_pmap(mdl, &PMAP)) {
        printf("error: matrix %s: out of memory\n", "PMAP");
        goto return_result;
    }
    if(!calc_elec_to_sys_map(mdl)) {
        printf("error: %s: out of memory\n", "elec_to_sys");
        goto return_result;
    }
    if(!calc_elec_zc_map(mdl)) {
        printf("error: %s: out of memory\n", "elec_zc");
        goto return_result;
    }
    for(int frame_idx = 0; frame_idx < mdl->n_params[1]; frame_idx++) {
        /* fill in the matrices and vectors */
        const int gnd_node = 1;
        if (!build_system_matrix(mdl, gnd_node, frame_idx, A, PMAP)) {
            printf("error: failed to build system matrix A\n");
            goto return_result;
        }
        // printf_matrix(A);
        if(!coo_to_csc(A)) {
            printf("error: failed to compress system matrix A\n");
            goto return_result;
        }
        cholmod_free_sparse (&As, &c);
        As = copy_csc_to_cholmod_sparse(A, &c);
        cholmod_check_sparse (As, &c);
        /* END MODIFIED: now continue with the usual operations */
        assert(As != NULL);
        assert(As->stype != 0);
        if(frame_idx == 0) {
            /* CHOLMOD analyze; Non-zeros do not change at each iteration
             * (its the same mesh), so we do this expensive step once. */
            L = cholmod_analyze (As, &c);
        }
        cholmod_factorize (As, L, &c);          /* factorize */
        /* create stimulus */
        assert(mdl->n_stimmeas != 0);
        for (int i = 0; i < mdl->n_stimmeas; i++) { /* for each unique stimulus */
#if 0 /* default CHOLMOD solver */
            const int xtype = CHOLMOD_REAL;
            cholmod_free_dense (&b, &c);
            b = cholmod_zeros(rows, 1, xtype, &c);   /* b = zeros(n,1) */
            calc_stim(mdl, i, b->x);
            calc_stim_gnd(mdl, gnd_node, b->x);
            double * bx = b->x;
            for(int j = 0; j < rows; j++) {
                printf("  b[%d]=%g\n", j, bx[j]);
            }
            x = cholmod_solve (CHOLMOD_A, L, b, &c);       /* solve Ax=b */
            /* find solution norm */
            double one [2] = {1, 0}, m1 [2] = {-1, 0} ;     /* scalars [real, imag]*/
            cholmod_sdmult (As, 0, m1, one, x, b, &c);      /* b = b-Ax */
            double norm = cholmod_norm_dense (b, 0, &c);
            if(norm > 1e-12) {
                printf ("error: bad forward solution at frame#%d, measurement#%d: ||b-Ax|| = %8.1e\n", frame_idx, i, norm);
                goto return_result;
            }
            meas[frame_idx * (mdl->n_stimmeas) + i] = calc_meas(mdl, i, x->x) * 10;
            double * xx = x->x;
            printf("potentials\n");
            for(int i = 0; i < rows; i++) {
                printf(" %18g\n", xx[i]);
            }
            cholmod_free_dense (&x, &c);
#else /* alternate CHOLMOD solver for sparse b and sparse x */
            const int n_elec = mdl->n_elec;
            double bs_x[rows];
            int bs_p[2] = {0, n_elec};
            int bs_i[n_elec];
            for(int i = 0; i < n_elec; i++) {
                const int frame_idx = mdl->elec_to_sys[i];
                bs_i[i] = frame_idx;
                bs_x[frame_idx] = 0;
            }
            calc_stim(mdl, i, &(bs_x[0]));
            calc_stim_gnd(mdl, gnd_node, &(bs_x[0]));
            cholmod_dense bs = {
                .nrow = rows, .ncol = 1, .nzmax = rows, .d = rows,
                .x = &(bs_x[0]), .z = NULL, .xtype = CHOLMOD_REAL, .dtype = CHOLMOD_DOUBLE,
            };
            cholmod_sparse bset = {
                .nrow = rows, .ncol = 1, .nzmax = n_elec,
                .p = &(bs_p[0]), .i = &(bs_i[0]),
                .nz = NULL, .x = NULL, .z = NULL,
                .stype = 0, .itype = CHOLMOD_INT, .xtype = CHOLMOD_PATTERN, .dtype = CHOLMOD_DOUBLE,
                .sorted = 0, .packed = 1
            };
            assert(cholmod_check_sparse (&bset, &c));
            assert(cholmod_check_dense (&bs, &c));
            int nret = cholmod_solve2 (CHOLMOD_A, L, &bs, &bset, &x, &xset, &Y, &E, &c);       /* solve Ax=b */
            assert(nret);
            meas[frame_idx * (mdl->n_stimmeas) + i] = calc_meas(mdl, i, x->x);
#endif
        }
    }
    ret = SUCCESS;
return_result:
    cholmod_free_factor (&L, &c);          /* free matrices */
    cholmod_free_sparse (&As, &c);
    cholmod_free_dense (&b, &c);
    cholmod_free_dense (&x, &c);
    cholmod_free_sparse(&xset, &c);
    cholmod_free_dense (&E, &c);
    cholmod_free_dense (&Y, &c);
    cholmod_finish (&c);               /* finish CHOLMOD */
    free_matrix(A);
    free_matrix(PMAP);
    return ret;
}

