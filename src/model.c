/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <string.h> /* bzero */
#include <stdlib.h> /* free, qsort */
#include <lapacke.h> /* inv: dgetrf, dgetri */
#include <math.h> /* fabs */
#include <limits.h> /* MAX_INT */
#include <assert.h>

#include "matrix.h"
#include "model.h"

/* references:
 * [1] A. Boyle, PhD Thesis, 2016, Geophysical Applications of Electrical Impedance Tomography, Carleton University
 */

#ifdef UNIT_TESTING
extern void * _test_malloc(const size_t size, const char * file, const int line);
extern void _test_free(void * const ptr, const char * file, const int line);
#define malloc(size) _test_malloc(size, __FILE__, __LINE__)
#define free(ptr) _test_free(ptr, __FILE__, __LINE__)
#endif

#define SUCCESS 1
#define FAILURE 0

static void init_mesh(mesh * m)
{
    bzero(m, sizeof(mesh));
}

static void init_model(model * m)
{
    bzero(m, sizeof(model));
}

model * malloc_model()
{
    model * m = malloc(sizeof(model));
    if(m) {
        init_model(m);
    }
    return m;
}

static void free_mesh(mesh * m)
{
    free(m->nodes);
    free(m->elems);
    free(m->matidx);
    free(m->surfaceelems);
    free(m->bc);
}

void * free_model(model * m)
{
    if(m == NULL) {
        return NULL;
    }
    free(m->data);
    free(m->params);
    free(m->stimmeas);
    free_mesh(&(m->fwd));
    init_model(m);
    free(m);
    return NULL;
}

int set_model_data(model * m, int rows, int cols)
{
    assert(rows > 0);
    assert(cols > 0);
    assert(m != NULL);
    m->n_data[0] = rows;
    m->n_data[1] = cols;
    free(m->data);
    m->data = malloc(sizeof(double) * rows * cols);
    if(m->data == NULL) {
        return FAILURE;
    }
    return SUCCESS;
}

int set_model_params(model * m, int rows, int cols)
{
    assert(rows > 0);
    assert(cols > 0);
    assert(m != NULL);
    m->n_params[0] = rows;
    m->n_params[1] = cols;
    free(m->params);
    m->params = malloc(sizeof(double) * rows * cols);
    if(m->params == NULL) {
        return FAILURE;
    }
    return SUCCESS;
}

int set_model_stimmeas(model * m, int rows)
{
    assert(rows > 0);
    assert(m != NULL);
    m->n_stimmeas = rows;
    free(m->stimmeas);
    m->stimmeas = malloc(sizeof(int) * rows * 4);
    if(m->stimmeas == NULL) {
        return FAILURE;
    }
    return SUCCESS;
}

void set_model_hp(model * m, double hp)
{
    assert(hp >= 0.0);
    assert(m != NULL);
    m->hp = hp;
}

void set_mesh_dim(mesh * m, int dim)
{
    assert(dim == 2 || dim == 3);
    assert(m != NULL);
    free_mesh(m);
    init_mesh(m);
    m->dim = dim;
}

int set_mesh_nodes(mesh * m, int n_nodes)
{
    assert(n_nodes >= 0);
    assert(m != NULL);
    assert(m->dim == 2 || m->dim == 3);
    m->n_nodes = n_nodes;
    free(m->nodes);
    m->nodes = malloc(sizeof(double) * n_nodes * m->dim);
    if(m->nodes == NULL) {
        return FAILURE;
    }
    return SUCCESS;
}

int set_mesh_elems(mesh * m, int n_elems)
{
    assert(n_elems >= 0);
    assert(m != NULL);
    assert(m->dim == 2 || m->dim == 3);
    m->n_elems = n_elems;
    free(m->elems);
    free(m->matidx);
    m->elems = malloc(sizeof(int) * n_elems * (m->dim + 1));
    m->matidx = malloc(sizeof(int) * n_elems);
    if((m->elems == NULL) || (m->matidx == NULL)) {
        return FAILURE;
    }
    return SUCCESS;
}

int set_mesh_surfaceelems(mesh * m, int n_se)
{
    assert(n_se >= 0);
    assert(m != NULL);
    assert(m->dim == 2 || m->dim == 3);
    m->n_se = n_se;
    free(m->surfaceelems);
    free(m->bc);
    m->surfaceelems = malloc(sizeof(int) * n_se * m->dim);
    m->bc = malloc(sizeof(int) * n_se);
    if((m->surfaceelems == NULL) || (m->bc == NULL)) {
        return FAILURE;
    }
    return SUCCESS;
}

double * inv(int n, double A[n][n])
{
    double * Ap = &(A[0][0]);
    int ipiv[n + 1];
    bzero(ipiv, (n + 1)*sizeof(int));
    if(LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, Ap, n, ipiv) ||
       LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, Ap, n, ipiv)) {
        return NULL;
    }
    return Ap;
}

double det(int n, double A[n][n])
{
    if(n == 2) {
        return A[0][0] * A[1][1] - A[0][1] * A[1][0];
    }
    double d = 0;
    double sign = 1.0;
    for(int i = 0; i < n; i++, sign *= -1.0) {
        /* construct sub-matrix B */
        double B[n - 1][n - 1];
        for(int j = 1; j < n; j++) {
            for(int k = 0, s = 0; k < n; k++) {
                if(k != i) {
                    B[j - 1][s] = A[j][k];
                    s++;
                }
            }
        }
        d += sign * A[0][i] * det(n - 1, B); /* (-1)^n * a1i * |B| */
    }
    return d;
}

/* shape function
 * Se from [1] eq (C.31) -- first-order tetrahedral/triangular elements
 * Se is symmetric, only the lower triangular portion is calculated (E1^T E1 is symmetric)
 * inputs: nd = # of dimensions (2=2D, 3=3D),
 *         N = ptrs to node coordinates (nd+1 nodes with nd double precision values)
 *             ex: {&n[120], &n[332], &n[556], &n[733]} for 3D, where n[120] = {1.1, 0.1, 3.2}, etc.
 *             NOTE: nodes should be in increasing order so that calc_Se_ij
 *                   gives an upper triangular matrix in the global node numbering
 * output: Se = replaces Se[0..t-1] for t = n(n+1)/2 entries with n=nd+1 elements as
 *         upper triangular symmetric matrix, row-major order
 * returns: NULL if matrix inverse fails, Se otherwise
 */
double * calc_Se_v(const int nd, double const * const * const nodes, double * Se)
{
    const int n = nd + 1;
    double E[n][n]; /* from [1] eq (C.33) */
    for(int i = 0; i < n; i++) {
        E[i][0] = 1;
        for(int j = 1; j < n; j++) {
            E[i][j] = nodes[i][j - 1];
        }
    }
    if(inv(n, E) == NULL) {
        return NULL;
    }
    const double detE = fabs(det(n, E));
    int ndf = 1; /* = !nd = factorial(number of dimensions) */
    for(int i = n - 1; i > 0; i-- ) {
        ndf *= i;
    }
    const double c = 1 / (ndf * detE);
    /* Se = c * E1^T E1, where E1 is E without the top row */
    int m = 0;
    for(int i = 0; i < n; i++ ) {
        for(int j = i; j < n; j++ ) {
            Se[m] = 0;
            for(int k = 1; k < n; k++) {
                Se[m] += E[k][i] * E[k][j];
            }
            /* TODO pull this out and multiply by conductivity D, to save some multiply operations */
            Se[m] *= c;
            m++;
        }
    }
    return Se;
}

void calc_Se_ij(const int nd, int const * const elem, int * ii, int * jj)
{
    const int n = nd + 1;
    int idx = 0;
    for(int i = 0; i < n; i++ ) {
        for(int j = i; j < n; j++ ) {
            ii[idx] = elem[i] - 1;
            jj[idx] = elem[j] - 1;
            idx++;
        }
    }
}

int calc_Se_n(const int nd)
{
    const int n = nd + 1;
    return n * (n + 1) / 2;
}

/* returns 0 on success, > 0 on failure as the singular element number e+1 */
int calc_Se(mesh const * const m, int * ii, int * jj, double * Se)
{
    const int nd = m->dim;
    const int n_elems = m->n_elems;
    const int * elems = m->elems;
    const double * nodes = m->nodes;
    const int n = calc_Se_n(nd);
    for(int i = 0; i < n_elems; i++) {
        calc_Se_ij(nd, elems, ii, jj);
        double const * node_list [4];
        for(int j = 0; j < nd + 1; j ++) {
            int idx = elems[j] - 1;
            if ( idx < 0 ) {
                idx = 0;
            }
            node_list[j] = &(nodes[idx * nd]);
        }
        if(calc_Se_v(nd, node_list, Se) == NULL) {
            return i + 1;
        }
        elems += (nd + 1);
        ii += n;
        jj += n;
        Se += n;
    }
    return 0;
}

/* Removes ground node 'gnd' from COO entries by searching ii and jj.
 * Returns number of entries modified or deleted */
int cmp_int( const void * a, const void * b )
{
    return *(int *)a - *(int *)b;
}
int calc_gnd(const int gnd, int * nnz, int * ii, int * jj, double * Se)
{
    const int gndidx = gnd - 1;
    int ret = 0;
    int idx [*nnz];
    for (int i = 0; i < *nnz; i++ ) {
        if ( ii[i] == gndidx || jj[i] == gndidx) {
            if (ret == 0) { /* replace with (gndidx,gndidx)=1.0 */
                idx[i] = i;
                ii[i] = gndidx;
                jj[i] = gndidx;
                Se[i] = 1.0;
            }
            else {
                idx[i] = INT_MAX; /* remove */
            }
            ret++;
        }
        else {
            idx[i] = i;
        }
    }
    if( ret <= 1 ) {
        return ret;
    }
    qsort(idx, *nnz, sizeof(int), &cmp_int);
    *nnz -= ret - 1;
    /* we only shift entries down in the indices, so
     * its safe to do straight copy without swap */
    for(int i = 0 ; i < *nnz; i++) {
        const int x = idx[i];
        ii[i] = ii[x];
        jj[i] = jj[x];
        Se[i] = Se[x];
    }
    return ret;
}

/* Construct a dense column vector for Neumann (current) stimulus 'amp' on
 * boundary m->bc = bc.
 * Returns number of nodes perturbed: 0 = failure. */
int calc_stim_neumann(mesh const * const m, double amp, int bc, int gnd, double * b)
{
    const int dim = m->dim;
    int cnt_bc_local_nodes = 0;
    for(int i = 0; i < m->n_se; i ++ ) {
        const int bci = m->bc[i];
        if (bci == bc) {
            cnt_bc_local_nodes += dim;
        }
    }
    const double amp_nodal = +amp / (double) cnt_bc_local_nodes;
    const int gndidx = gnd - 1; /* convert node number to index in 'b' */
    for(int i = 0; i < m->n_se; i ++ ) {
        const int bci = m->bc[i];
        if (bci != bc) {
            continue;
        }
        for(int j = 0; j < dim; j++) {
            int idx = m->surfaceelems[i * dim + j] - 1;
            if ( idx == gndidx ) {
                /* ground node: do not apply current */
                cnt_bc_local_nodes--;
            }
            else {
                b[ idx ] += amp_nodal;
            }
        }
    }
    return cnt_bc_local_nodes;
}
