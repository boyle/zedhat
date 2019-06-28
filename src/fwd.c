/* Copyright 2019, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <assert.h> /* assert */
#include <string.h> /* memcpy, memset */
#include <cholmod.h>
#include "matrix.h" /* coo_to_csc */

#include "fwd.h"

#define SUCCESS 1
#define FAILURE 0

static void calc_stim_gnd(model const * const m, int gnd, double * b)
{
    assert(gnd > 0);
    assert(gnd <= m->fwd.n_nodes);
    b[gnd - 1] = 0.0;
}

static void calc_stim(model const * const m, int idx, double * b)
{
    assert(b != NULL);
    assert(idx >= 0);
    assert(idx < m->n_stimmeas);
    const int A = m->stimmeas[idx * 4 + 0]; /* current+ */
    const int B = m->stimmeas[idx * 4 + 1]; /* current- */
    const int An = m->elec_to_sys[A - 1];
    const int Bn = m->elec_to_sys[B - 1];
    b[An] = +1;
    b[Bn] = -1;
}

static int cmp_stim(model const * const mdl, int stim1, int stim2)
{
    const int A1 = mdl->stimmeas[stim1 * 4 + 0]; /* current+ */
    const int B1 = mdl->stimmeas[stim1 * 4 + 1]; /* current- */
    const int A2 = mdl->stimmeas[stim2 * 4 + 0]; /* current+ */
    const int B2 = mdl->stimmeas[stim2 * 4 + 1]; /* current- */
    return (A1 == A2) ? (B1 - B2) : (A1 - A2);
}

static double calc_meas(model const * const m, int idx, double * x)
{
    assert(x != NULL);
    assert(idx >= 0);
    assert(idx < m->n_stimmeas);
    const int M = m->stimmeas[idx * 4 + 2]; /* voltage+ */
    const int N = m->stimmeas[idx * 4 + 3]; /* voltage- */
    const int Mn = m->elec_to_sys[M - 1];
    const int Nn = m->elec_to_sys[N - 1];
    const double gain = m->measgain[idx];
    return (x[Mn] - x[Nn]) * gain;
}

static int build_system_matrix(const model * m, int gnd, int frame_idx, matrix * A, matrix ** Acopy, matrix * PMAP)
{
    assert(m != NULL);
    assert(A != NULL);
    assert(frame_idx >= 0);
    assert(gnd > 0);
    assert(gnd < m->fwd.n_nodes);
    /* build shape matrices */
    /* TODO we could reuse a copy of the expensive to compute shape elements at
     * each iteration, intead of recomputing from scratch here */
    const int rows = calc_sys_size(m); /* = cols */
    size_t nnz = calc_sys_nnz(m);
    const int ne = m->fwd.n_elems;
    const int dim = m->fwd.dim;
    const int se_n = calc_sys_elem_n(dim); /* sparse matrix entries per mesh element */
    if(!malloc_matrix_data(A, COO_SYMMETRIC, rows, rows, nnz)) {
        return FAILURE;
    }
    int ret = calc_sys_elem(&(m->fwd), A->x.sparse.ia, A->x.sparse.ja, A->x.sparse.a);
    if (ret != 0) {
        printf("error: bad forward model element#%d\n", ret);
        return FAILURE;
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
    /* store shape matrices for later */
    if((Acopy != NULL) && (*Acopy == NULL)) {
        if(!copy_matrix(A, Acopy, "Asys")) {
            return FAILURE;
        }
    }
    /* apply conductivity */
    if(m->n_params[1] == 0) { /* assume conductivity = 1 S/m */
        return SUCCESS;
    }
    assert(frame_idx < m->n_params[1]);
    if(PMAP == NULL) {
        assert(ne == m->n_params[0]);
        const int rows = m->n_params[0];
        for(int e = 0; e < ne; e++) {
            const double cond = m->params[e + rows * frame_idx]; /* conductivity */
            for(int i = 0; i < se_n; i++) {
                A->x.sparse.a[e * se_n + i] *= cond;
            }
        }
    }
    else { /* map parameters to elements */
        assert(PMAP->n == m->n_params[0]);
        assert(PMAP->m == ne);
        assert(PMAP->type == CSR);
        const int rows = m->n_params[0];
        for(int e = 0; e < ne; e++) {
            double cond_frac = 0.0;
            for(int j = PMAP->x.sparse.ia[e]; j < PMAP->x.sparse.ia[e + 1]; j++) {
                const int param = PMAP->x.sparse.ja[j];
                const double element_fraction = PMAP->x.sparse.a[j];
                const double conductivity = m->params[param + rows * frame_idx];
                cond_frac += conductivity * element_fraction;
            }
            // printf("[e=%d] σ=%g\n", e, cond_frac);
            for(int i = 0; i < se_n; i++) {
                A->x.sparse.a[e * se_n + i] *= cond_frac;
            }
        }
    }
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
    const int nnz = mdl->fwd.n_pmap;
    const int cols = mdl->n_params[0];
    if(!malloc_matrix_data(*PMAP, COO, rows, cols, nnz)) {
        return FAILURE;
    }
    for(int i = 0; i < nnz; i++) {
        assert(mdl->fwd.pmap_elem[i] <= rows);
        assert(mdl->fwd.pmap_param[i] <= cols);
        assert(mdl->fwd.pmap_elem[i] > 0);
        assert(mdl->fwd.pmap_param[i] > 0);
        assert(mdl->fwd.pmap_frac[i] <= 1.0);
        assert(mdl->fwd.pmap_frac[i] >= 0.0);
        (*PMAP)->x.sparse.a[i] = mdl->fwd.pmap_frac[i]; /* TODO faster w/ memcpy? */
        (*PMAP)->x.sparse.ia[i] = mdl->fwd.pmap_elem[i] - 1;
        (*PMAP)->x.sparse.ja[i] = mdl->fwd.pmap_param[i] - 1;
    }
    if(!coo_to_csr(*PMAP)) {
        return FAILURE;
    }
    return SUCCESS;
}

typedef struct {
    int is_initialized;
    int last_frame;
    int rows;
    cholmod_common c;
    cholmod_dense * x;
    cholmod_dense * b;
    cholmod_factor * L;
    cholmod_dense * E;
    cholmod_dense * Y;
    cholmod_sparse * xset;
    cholmod_sparse * As;
    matrix * A;
    matrix * Asys;
    matrix * PMAP;
} fwdsolve_state;

static void free_state(fwdsolve_state * state)
{
    cholmod_free_factor (&state->L, &state->c);          /* free matrices */
    cholmod_free_sparse (&state->As, &state->c);
    cholmod_free_dense (&state->b, &state->c);
    cholmod_free_dense (&state->x, &state->c);
    cholmod_free_sparse(&state->xset, &state->c);
    cholmod_free_dense (&state->E, &state->c);
    cholmod_free_dense (&state->Y, &state->c);
    cholmod_finish (&state->c);               /* finish CHOLMOD */
    state->A = free_matrix(state->A);
    state->Asys = free_matrix(state->Asys);
    state->PMAP = free_matrix(state->PMAP);
    state->is_initialized = 0;
}

static int init_state(model * mdl, fwdsolve_state * state)
{
    assert(mdl != NULL);
    assert(!state->is_initialized);
    if(!check_model(mdl)) {
        return FAILURE;
    }
    /* CHOLMOD sparse symmetric decomposition */
    cholmod_start (&state->c);                /* start CHOLMOD */
    state->A = malloc_matrix();
    if(!malloc_matrix_name(state->A, "system", "A", "V/A")) {
        printf("error: matrix %s: out of memory\n", "A");
        return FAILURE;
    }
    if(!build_pmap(mdl, &state->PMAP)) {
        printf("error: matrix %s: out of memory\n", "PMAP");
        return FAILURE;
    }
    if(!calc_elec_to_sys_map(mdl)) {
        printf("error: %s: out of memory\n", "elec_to_sys");
        return FAILURE;
    }
    if(!calc_elec_zc_map(mdl)) {
        printf("error: %s: out of memory\n", "elec_zc");
        return FAILURE;
    }
    state->rows = calc_sys_size(mdl);
    state->last_frame = -1;
    state->is_initialized = 1;
    return SUCCESS;
}

int fwd_solve(model * mdl, double * meas)
{
    assert(mdl != NULL);
    assert(meas != NULL);
    int ret = FAILURE;
    fwdsolve_state state = {0};
    if(!init_state(mdl, &state)) {
        goto return_result;
    }
    for(int frame_idx = 0; frame_idx < mdl->n_params[1]; frame_idx++) {
        /* fill in the matrices and vectors */
        const int gnd_node = 1;
        if (!build_system_matrix(mdl, gnd_node, frame_idx, state.A, NULL, state.PMAP)) {
            printf("error: failed to build system matrix A\n");
            goto return_result;
        }
        // printf_matrix(state.A);
        if(!coo_to_csc(state.A)) {
            printf("error: failed to compress system matrix A\n");
            goto return_result;
        }
        cholmod_free_sparse (&state.As, &state.c);
        state.As = copy_csc_to_cholmod_sparse(state.A, &state.c);
        cholmod_check_sparse (state.As, &state.c);
        /* END MODIFIED: now continue with the usual operations */
        assert(state.As != NULL);
        assert(state.As->stype != 0);
        if(frame_idx == 0) {
            /* CHOLMOD analyze; Non-zeros do not change at each iteration
             * (its the same mesh), so we do this expensive step once. */
            state.L = cholmod_analyze (state.As, &state.c);
        }
        cholmod_factorize (state.As, state.L, &state.c);          /* factorize */
        /* create stimulus */
        assert(mdl->n_stimmeas != 0);
        const int rows = calc_sys_size(mdl);
        for (int i = 0; i < mdl->n_stimmeas; i++) { /* for each unique stimulus */
            if((i > 0) && (cmp_stim(mdl, i, i - 1) == 0)) { /* same stimulus as last stimmeas row */
                meas[frame_idx * (mdl->n_stimmeas) + i] = calc_meas(mdl, i, state.x->x);
                continue;
            }
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
            assert(cholmod_check_sparse (&bset, &state.c));
            assert(cholmod_check_dense (&bs, &state.c));
            int nret = cholmod_solve2 (CHOLMOD_A, state.L, &bs, &bset, &state.x, &state.xset, &state.Y, &state.E, &state.c);       /* solve Ax=b */
            assert(nret);
            meas[frame_idx * (mdl->n_stimmeas) + i] = calc_meas(mdl, i, state.x->x);
#endif
        }
    }
    ret = SUCCESS;
    // printf_model(mdl);
return_result:
    free_state(&state);
    return ret;
}

static int fwd_solve_node_voltages(model * mdl, int frame_idx, int stim_idx, double * nodal, fwdsolve_state * state)
{
    assert(mdl != NULL);
    assert(frame_idx >= 0);
    assert(stim_idx >= 0);
    assert(stim_idx < mdl->n_stimmeas);
    assert(nodal != NULL);
    assert(state != NULL);
    assert(state->is_initialized);
    if(state->last_frame != frame_idx) {
        /* fill in the matrices and vectors */
        const int gnd_node = 1;
        if (!build_system_matrix(mdl, gnd_node, frame_idx, state->A, &(state->Asys), state->PMAP)) {
            printf("error: failed to build system matrix A\n");
            return FAILURE;
        }
        // printf_matrix(state->Asys);
        // printf_matrix(state->A);
        if(!coo_to_csc(state->A)) {
            printf("error: failed to compress system matrix A\n");
            return FAILURE;
        }
        // printf_matrix(state->A);
        cholmod_free_sparse (&(state->As), &state->c);
        state->As = copy_csc_to_cholmod_sparse(state->A, &state->c);
        cholmod_check_sparse (state->As, &state->c);
        /* END MODIFIED: now continue with the usual operations */
        assert(state->As != NULL);
        assert(state->As->stype != 0);
        if(state->L == NULL) {
            /* CHOLMOD analyze; Non-zeros do not change at each iteration
             * (its the same mesh), so we do this expensive step once. */
            state->L = cholmod_analyze (state->As, &state->c);
        }
        cholmod_factorize (state->As, state->L, &state->c);          /* factorize */
        /* create stimulus */
        assert(mdl->n_stimmeas != 0);
        state->last_frame = frame_idx;
    }
    const int gnd_node = 1;
    const int rows = state->rows;
#if 1 /* default CHOLMOD solver */
    const int xtype = CHOLMOD_REAL;
    cholmod_free_dense (&state->b, &state->c);
    cholmod_free_dense (&state->x, &state->c);
    state->b = cholmod_zeros(rows, 1, xtype, &state->c);   /* b = zeros(n,1) */
    calc_stim(mdl, stim_idx, state->b->x);
    calc_stim_gnd(mdl, gnd_node, state->b->x);
    state->x = cholmod_solve (CHOLMOD_A, state->L, state->b, &state->c);       /* solve Ax=b */
    /* find solution norm */
    // double one [2] = {1, 0}, m1 [2] = {-1, 0} ;     /* scalars [real, imag]*/
    // cholmod_sdmult (state->As, 0, m1, one, state->x, state->b, &state->c);      /* b = b-Ax */
    // double norm = cholmod_norm_dense (state->b, 0, &state->c);
    // if(norm > 1e-12) {
    //     printf ("error: bad forward solution at frame#%d, measurement#%d: ||b-Ax|| = %8.1e\n", frame_idx, stim_idx, norm);
    //     return FAILURE;
    // }
    double * xx = state->x->x;
    for(int i = 0; i < rows; i++) {
        nodal[i] = xx[i];
    }
    cholmod_free_dense (&state->x, &state->c);
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
    calc_stim(mdl, stim_idx, &(bs_x[0]));
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
    assert(cholmod_check_sparse (&bset, &state->c));
    assert(cholmod_check_dense (&bs, &state->c));
    int nret = cholmod_solve2 (CHOLMOD_A, state->L, &bs, &bset, &(state->x), &(state->xset), &(state->Y), &(state->E), &state->c);       /* solve Ax=b */
    assert(nret);
    double * xx = state->x->x;
    for(int i = 0; i < rows; i++) {
        nodal[i] = xx[i];
    }
#endif
    return SUCCESS;
}

static void sdmult(const model * mdl, const int elem_idx, const matrix * A, double a, double * x, double * b) /* b += a(Ax) */
{
    assert(mdl != NULL);
    assert(A != NULL);
    assert(A->type == COO);
    assert(A->symmetric);
    const int n = calc_sys_elem_n(mdl->fwd.dim); /* matrix entries per element */
    const int idx_start = elem_idx * n;
    const int idx_end = (elem_idx + 1) * n;
    for(int idx = idx_start; idx < idx_end; idx++) {
        const int col =  A->x.sparse.ja[idx];
        const int row =  A->x.sparse.ia[idx];
        const double aAx = a * A->x.sparse.a[idx];
        b[row] += aAx * x[col];
        if(col != row) {
            b[col] +=  aAx * x[row];
        }
    }
}

static double fwd_solve_node_stim(const model * mdl, double * b, const int meas_idx, fwdsolve_state * state)
{
    assert(state != NULL);
    assert(state->is_initialized);
    cholmod_free_dense (&state->b, &state->c);
    cholmod_free_dense (&state->x, &state->c);
    const int xtype = CHOLMOD_REAL;
    state->b = cholmod_allocate_dense(state->rows, 1, state->rows, xtype, &state->c);
    memcpy(state->b->x, b, sizeof(double) * state->rows);
    state->x = cholmod_solve (CHOLMOD_A, state->L, state->b, &state->c);       /* solve Ax=b */
    return calc_meas(mdl, meas_idx, state->x->x);
}

#define colmaj(i,j) ((i) + (j)*(rows)) /* column major indexing */
int calc_jacobian(model * mdl, double * J)
{
    /* Jacobian Equation:
     * eqn C.36 A. Boyle 2016 Phd Thesis
     *   J = - T A^-1 dA/dc A^-1 b
     * for measurement selection T, system matrix A, stimulus b, and
     * partial derivative of system matrix with conductivity dA/dc.
     * Note that A^-1 b = x, for nodal voltages x.
     * Note that dA/dc is all zeros except
     * for the element that is changing so we can efficiently leverage the COO
     * system matrix, before conductivity is applied, to grab an element at a
     * time for a unit conductivity change.
     */
    assert(mdl != NULL);
    assert(J != NULL);
    int ret = FAILURE;
    double * x = NULL;
    fwdsolve_state state = {0};
    if(!init_state(mdl, &state)) {
        goto return_result;
    }
    const matrix * PMAP = (state.PMAP);
    const int len = state.rows;
    const int cols = mdl->n_params[0];
    const int rows = mdl->n_stimmeas;
    const int is_pmap = (mdl->fwd.n_pmap > 0);
    if(is_pmap) {
        assert(PMAP->type == CSR);
        assert(PMAP->m == mdl->fwd.n_elems);
        assert(PMAP->n == mdl->n_params[0]);
    }
    memset(J, 0, sizeof(double) * rows * cols);
    x = malloc(sizeof(double) * len);
    assert(x != NULL);
    for(int i = 0; i < rows; i++) { /* i: Jacobian row# = meas#*/
        // if((i == 0) || (cmp_stim(mdl, i, i - 1) != 0)) { /* different stimulus from last stimmeas row */
        memset(x, 0, sizeof(double)*len);
        if(!fwd_solve_node_voltages(mdl, 0, i, x, &state)) {
            goto return_result;
        }
        // }
        for(int e = 0; e < mdl->fwd.n_elems; e++) { /* e: element# */
            double b[len];
            memset(b, 0, sizeof(double)*len);
            sdmult(mdl, e, state.Asys, -1.0, x, b);
            const double meas = fwd_solve_node_stim(mdl, b, i, &state);
            if(is_pmap) {
                for(int k = PMAP->x.sparse.ia[e]; k < PMAP->x.sparse.ia[e + 1]; k++) {
                    const int j = PMAP->x.sparse.ja[k]; /* j: Jacobian col#  = param# */
                    const double element_fraction = PMAP->x.sparse.a[k];
                    const double param = 1.0; // mdl->n_params[0] > 0 ? mdl->params[j] : 1.0;
                    J[colmaj(i, j)] += meas * param * element_fraction;
                }
            }
            else {
                const int j = e; /* e: Jacobian col# = elem# */
                const double param = 1.0; // mdl->n_params[0] > 0 ? mdl->params[j] : 1.0;
                J[colmaj(i, j)] = meas * param;
            }
        }
    }
    ret = SUCCESS;
return_result:
    free(x);
    free_state(&state);
    return ret;
}
