/* Copyright 2018,2019, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <string.h> /* bzero */
#include <stdlib.h> /* free */
#include <lapacke.h> /* inv: dgetrf, dgetri */
#include <math.h> /* fabs, sqrt */
#include <limits.h> /* MAX_INT */
#include <assert.h> /* assert */
#include <stdio.h> /* printf */

#include "matrix.h"
#include "model.h"

/* references:
 * [1] A. Boyle, PhD Thesis, 2016, Geophysical Applications of Electrical Impedance Tomography, Carleton University
 */

#ifdef UNIT_TESTING
extern void * _mock_test_malloc(const size_t size, const char * file, const int line);
extern void _test_free(void * const ptr, const char * file, const int line);
#define malloc(size) _mock_test_malloc(size, __FILE__, __LINE__)
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

void reset_model(model * m)
{
    assert(m != NULL);
    free(m->data);
    free(m->params);
    free(m->stimmeas);
    free(m->elec_to_sys);
    free_mesh(&(m->fwd));
    init_model(m);
}

void * free_model(model * m)
{
    if(m == NULL) {
        return NULL;
    }
    reset_model(m);
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
    m->n_elec = 0;
    free(m->stimmeas);
    free(m->elec_to_sys);
    m->elec_to_sys = NULL;
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
 *             NOTE: nodes should be in increasing order so that calc_sys_elem_ij
 *                   gives an upper triangular matrix in the global node numbering
 * output: Se = replaces Se[0..t-1] for t = n(n+1)/2 entries with n=nd+1 elements as
 *         upper triangular symmetric matrix, row-major order
 * returns: NULL if matrix inverse fails, Se otherwise
 */
static double calc_elem_vol(const int n, double const * const * const nodes, double E[n][n])
{
    for(int i = 0; i < n; i++) {
        E[i][0] = 1;
        for(int j = 1; j < n; j++) {
            E[i][j] = nodes[i][j - 1];
        }
    }
    if(inv(n, E) == NULL) {
        return -1.0;
    }
    return fabs(det(n, E));
}

static double * calc_sys_elem_v(const int nd, double const * const * const nodes, double * Se)
{
    const int n = nd + 1;
    double E[n][n]; /* from [1] eq (C.33) */
    const double detE = calc_elem_vol(n, nodes, E);
    if(detE < 0) {
        return NULL;
    }
    int ndf = 1; /* = !nd = factorial(number of dimensions) */
    for(int i = 2; i <= nd; i++ ) {
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

static void calc_sys_elem_ij(const int nd, int const * const elem, int * ii, int * jj)
{
    const int n = nd + 1;
    int idx = 0;
    for(int i = 0; i < n; i++ ) {
        for(int j = i; j < n; j++ ) {
            int ei = elem[i] - 1;
            int ej = elem[j] - 1;
            /* swap indices (symmetric matrix) if this entry would be in the
             * upper triangular portion of the matrix */
            if(ei < ej) {
                int tmp = ei;
                ei = ej;
                ej = tmp;
            }
            ii[idx] = ei;
            jj[idx] = ej;
            idx++;
        }
    }
}

int calc_sys_elem_n(const int nd)
{
    /* nd=2: 6, nd=3: 10 ... for symmetric (isotropic conductivity)*/
    const int n = nd + 1;
    return n * (n + 1) / 2;
}

/* returns 0 on success, > 0 on failure as the singular element number e+1 */
int calc_sys_elem(mesh const * const m, int * ii, int * jj, double * Se)
{
    const int nd = m->dim;
    const int n_elems = m->n_elems;
    const int * elems = m->elems;
    const double * nodes = m->nodes;
    const int n = calc_sys_elem_n(nd);
    for(int i = 0; i < n_elems; i++) {
        calc_sys_elem_ij(nd, elems, ii, jj);
        double const * node_list [4];
        for(int j = 0; j < nd + 1; j ++) {
            int idx = elems[j] - 1;
            if ( idx < 0 ) {
                idx = 0;
            }
            node_list[j] = &(nodes[idx * nd]);
        }
        if(calc_sys_elem_v(nd, node_list, Se) == NULL) {
            return i + 1;
        }
        elems += (nd + 1);
        ii += n;
        jj += n;
        Se += n;
    }
    return 0;
}

void calc_sys_gnd(const int gnd, size_t nnz, int * ii, int * jj, double * Se)
{
    const int gndidx = gnd - 1;
    int ret = 0;
    for (int i = 0; i < nnz; i++ ) {
        const int is_row = (ii[i] == gndidx);
        const int is_col = (jj[i] == gndidx);
        if ( is_row || is_col) {
            if(ret == 0) {
                /* replace with (gndidx,gndidx)=1.0 */
                ii[i] = gndidx;
                jj[i] = gndidx;
                Se[i] = 1.0;
            }
            else {
                /* drop this entry */
                Se[i] = 0.0;
            }
            ret++;
        }
    }
    assert(ret != 0);
}

#define dot2(a,b) (((a[0])*(b[0])) + ((a[1])*(b[1])))
#define dot3(a,b) (((a[0])*(b[0])) + ((a[1])*(b[1])) + ((a[2])*(b[2])))
#define sum_sq2(a) dot2(a,a)
#define sum_sq3(a) dot3(a,a)
#define cross3i(a,b) (((a[1])*(b[2])) - ((a[2])*(b[1])))
#define cross3j(a,b) (((a[2])*(b[0])) - ((a[0])*(b[2])))
#define cross3k(a,b) (((a[0])*(b[1])) - ((a[1])*(b[0])))

// static void printf_vec(const int nd, double x[nd])
// {
//     printf(" (");
//     int i;
//     for(i = 0; i < nd - 1; i++) {
//         printf("%g, ", x[i]);
//     }
//     printf("%g", x[i]);
//     printf(")\n");
// }

static double calc_elem_area(const int nd, double nodes[nd][nd])
{
    /* len = diff(nodes) .. point to point distances */
    // for(int i = 0; i < nd; i++) {
    //     printf("[%d] ", i); printf_vec(nd, nodes[i]);
    // }
    for(int i = 0; i < nd - 1; i++) {
        for(int j = 0; j < nd; j++) {
            nodes[i][j] -= nodes[i + 1][j]; /* subtract the row below */
        }
        // printf("[%d]-[%d] ", i, i + 1); printf_vec(nd, nodes[i]);
    }
    if(nd == 3) {
        /* cross product of lengths c = a x b*/
        // printf("%gi\n",nodes[0][1]*nodes[1][2] - nodes[0][2]*nodes[1][1]);
        // printf("%gj\n",nodes[0][2]*nodes[1][0] - nodes[0][0]*nodes[1][2]);
        // printf("%gk\n",nodes[0][0]*nodes[1][1] - nodes[0][1]*nodes[1][0]);
        const double tmp[3] = {
            cross3i(nodes[0], nodes[1]),
            cross3j(nodes[0], nodes[1]),
            cross3k(nodes[0], nodes[1]),
        };
        // printf("[0] x [1]: "); printf_vec(nd, tmp);
        assert(dot3(tmp, nodes[0]) == 0); /* confirm c . a = 0 */
        assert(dot3(tmp, nodes[1]) == 0); /* confirm c . b = 0 */
        nodes[0][0] = sum_sq3(tmp);
        // printf("∑[0]²: "); printf_vec(1, nodes[0]);
    }
    else {
        nodes[0][0] = sum_sq2(nodes[0]);
        // printf("∑[0]²: "); printf_vec(1, nodes[0]);
    }
    const double area = sqrt(nodes[0][0]) / (double)(nd - 1);
    return area;
}

static void copy_nodes_from_surfaceelems(mesh const * const m, const int se_idx, const int nd, double nodes[nd][nd])
{
    const int dim = m->dim;
    for(int i = 0; i < dim; i++) {
        const int nidx = m->surfaceelems[se_idx * dim + i] - 1; /* node index */
        for(int j = 0; j < dim; j++) {
            nodes[i][j] = m->nodes[nidx * dim + j];
        }
    }
}

static int max_bc(mesh const * const m)
{
    int max = -1;
    for(int i = 0; i < m->n_se; i++ ) {
        const int bci = m->bc[i];
        if (bci > max) {
            max = bci;
        }
    }
    return max;
}

int calc_elec_to_sys_map(model * m)
{
    assert(m != NULL);
    free(m->elec_to_sys);
    m->n_elec = 0;
    m->elec_to_sys = NULL;
    const int dim = m->fwd.dim;
    const int n_bc = max_bc(&(m->fwd)) + 1;
    if(n_bc < 1) {
        return SUCCESS;
    }
    int cnt[n_bc];
    int nodes[n_bc];
    int pem_nodes[n_bc];
    for(int i = 0; i < n_bc; i++) {
        cnt[i] = 0;
        nodes[i] = 0;
        pem_nodes[i] = 0;
    }
    for(int i = 0; i < m->fwd.n_se; i ++ ) {
        const int bci = m->fwd.bc[i];
        cnt[bci]++;
        for(int j = 0; j < dim; j++) {
            const int nidx = m->fwd.surfaceelems[i * dim + j];
            if(nidx > 0) {
                nodes[bci]++;
                pem_nodes[bci] = nidx - 1;
            }
        }
    }
    /* we assume the BC=0 is the generic boundary with no electrodes, so start at i=1 */
    int n_elec = 0;
    for(int i = 1; i < n_bc; i++) {
        if(cnt[i] > 0) {
            n_elec++;
        }
    }
    m->elec_to_sys = malloc(sizeof(int) * n_elec);
    if(m->elec_to_sys == NULL) {
        return FAILURE;
    }
    m->n_elec = n_elec;
    int n_cem = 0;
    int electrode = 0;
    for(int i = 1; i < n_bc; i++) {
        if((cnt[i] == 1) && (nodes[i] == 1)) { /* PEM */
            m->elec_to_sys[electrode++] = pem_nodes[i];
        }
        else if(cnt[i] > 0) { /* CEM */
            m->elec_to_sys[electrode++] = n_cem + m->fwd.n_nodes;
            n_cem++;
        }
    }
    assert(electrode == n_elec);
    return SUCCESS;
}

static int count_se_for_this_bc(mesh const * const m, int bc)
{
    int n = 0;
    for(int i = 0; i < m->n_se; i ++ ) {
        const int bci = m->bc[i];
        if (bci == bc) {
            n++;
        }
    }
    return n;
}

/* Construct a dense column vector for Neumann (current) stimulus 'amp' on
 * boundary m->bc = bc.
 * Returns number of nodes perturbed: 0 = failure. */
int calc_stim_neumann(mesh const * const m, double amp, int bc, double * b)
{
    assert(m != NULL);
    assert(b != NULL);
    const int dim = m->dim;
    const int n = count_se_for_this_bc(m, bc);
    int cnt_bc_local_nodes = n * dim;
    int list[n];
    double area[n];
    double total_area = 0;
    for(int i = 0, j = 0; i < m->n_se; i ++ ) {
        const int bci = m->bc[i];
        if (bci != bc) {
            continue;
        }
        double tmp[dim][dim];
        copy_nodes_from_surfaceelems(m, i, dim, tmp);
        area[j] = calc_elem_area(dim, tmp);
        total_area += area[j];
        // printf("bc#%d: area[%d]=%g\n", bc, j, area[j]);
        list[j] = i;
        j++;
    }
    /* set boundary condition */
    for(int j = 0; j < n; j++ ) {
        const int i = list[j];
        for(int k = 0; k < dim; k++) {
            const int nidx = m->surfaceelems[i * dim + k] - 1;
            // printf("#%d nidx=%d %g*%g/%g %d\n",j,nidx,amp, area[j],total_area,dim);
            b[ nidx ] += amp * (area[j] / total_area) / (double) dim;
        }
    }
    return cnt_bc_local_nodes;
}

void calc_stim_gnd(model const * const m, int gnd, double * b)
{
    assert(gnd > 0);
    assert(gnd <= m->fwd.n_nodes);
    b[gnd - 1] = 0.0;
}

void calc_stim(model const * const m, int idx, double * b)
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

double calc_meas(model const * const m, int idx, double * x)
{
    assert(x != NULL);
    assert(idx >= 0);
    assert(idx < m->n_stimmeas);
    const int M = m->stimmeas[idx * 4 + 2]; /* voltage+ */
    const int N = m->stimmeas[idx * 4 + 3]; /* voltage- */
    const int Mn = m->elec_to_sys[M - 1];
    const int Nn = m->elec_to_sys[N - 1];
    return x[Mn] - x[Nn];
}

int calc_sys_cem_n(const int nd)
{
    /* nn = 2D: 3x3 sym, 3D: 4x4 sym
     * n_{bc=i} = ∑(n_se_{bc=i})x(nn) for n_se_{bc=i} = number of surface elements with BC=i
     * total CEM system matrix elements is then
     *   n = ∑_{i=1..ne} n_{bc=i} */
    return (nd == 2) ? 6 : 10;
}

int calc_sys_cem(mesh const * const m, int bc, double zc, int cem_node, int Se_idx, int * ii, int * jj, double * Se)
{
    /* For each segment of the CEM add
     * 2D:    cem   --   n1    n2
     * A_2d = 6     0    -3    -3   * 1/6 * length / zc
     *        0     0     0     0
     *       -3     0     2     1
     *       -3     0     1     2
     * 3D:    cem   --   n1    n2    n3
     * A_3d = 12    0    -4    -4    -4 * 1/12 * area / zc
     *         0    0     0     0     0
     *        -4    0     2     1     1
     *        -4    0     1     2     1
     *        -4    0     1     1     2
     */
    assert(m != NULL);
    const int Se_idx_start = Se_idx;
    const int dim = m->dim;
    const int n = count_se_for_this_bc(m, bc);
    int list[n];
    double area[n];
    double total_area = 0;
    for(int i = 0, j = 0; i < m->n_se; i ++ ) {
        const int bci = m->bc[i];
        if (bci != bc) {
            continue;
        }
        double tmp[dim][dim];
        copy_nodes_from_surfaceelems(m, i, dim, tmp);
        area[j] = calc_elem_area(dim, tmp);
        total_area += area[j];
        // printf("bc#%d: area[%d]=%g\n", bc, j, area[j]);
        list[j] = i;
        j++;
    }
    for(int j = 0; j < n; j++ ) {
        const int i = list[j];
        const double scale = (dim == 2) ? 6.0 : 12.0;
        const double area_zc = area[j] / zc;
        Se[Se_idx] = area_zc;
        ii[Se_idx] = cem_node;
        jj[Se_idx] = cem_node;
        Se_idx++;
        for(int k = 0; k < dim; k++) {
            const int nidx = m->surfaceelems[i * dim + k] - 1;
            Se[Se_idx] = ((dim == 2) ? -3.0 : -4.0) / scale * area_zc;
            ii[Se_idx] = cem_node;
            jj[Se_idx] = nidx;
            Se_idx++;
            Se[Se_idx] = +2.0 / scale * area_zc;
            ii[Se_idx] = nidx;
            jj[Se_idx] = nidx;
            Se_idx++;
            for(int p = 0; p < k; p++) {
                const int pidx = m->surfaceelems[i * dim + p] - 1;
                Se[Se_idx] = +1.0 / scale * area_zc;
                ii[Se_idx] = pidx;
                jj[Se_idx] = nidx;
                Se_idx++;
            }
        }
    }
    const int nnz_added = calc_sys_cem_n(dim) * n;
    assert(Se_idx - Se_idx_start == nnz_added);
    return nnz_added;
}

enum se_type {UNKNOWN, NONE, PEM, CEM};
static int check_bc(mesh const * const msh)
{
    assert(msh != NULL);
    const int dim = msh->dim;
    const int n_bc = max_bc(msh) + 1;
    enum se_type type[n_bc];
    for(int i = 0; i < n_bc; i++) {
        type[i] = UNKNOWN;
    }
    for(int i = 0; i < msh->n_se; i ++ ) {
        const int bci = msh->bc[i];
        /* a PEM is recorded in surfaceelems with a single node
         * as the first entry and the remaining entries as zero;
         * so we test the second entry to see if its zero */
        int cnt = 0;
        for(int j = 0; j < msh->dim; j++) {
            cnt += (msh->surfaceelems[i * dim + j] > 0) ? 1 : 0;
        }
        enum se_type newtype;
        if(cnt == msh->dim) {
            newtype = CEM;
        }
        else if(cnt == 1) {
            newtype = PEM;
        }
        else {
            printf("error: bad number of nodes on boundary %d\n", bci);
            return FAILURE;    /* bad node#s */
        }
        if(type[bci] == UNKNOWN) {
            type[bci] = newtype;
        }
        else if(type[bci] != newtype) {
            printf("error: mixed PEM and CEM on boundary %d\n", bci);
            return FAILURE;    /* mixed PEM and CEM syntax in same  BC */
        }
        else if(newtype == PEM) {
            printf("error: multiple PEM surface elements have the same boundary#%d, see surface element %d\n", bci, i + 1);
            return FAILURE;    /* duplicate PEM for same BC */
        }
        if((newtype == PEM) && (msh->surfaceelems[i * dim + 0] <= 0)) {
            printf("error: for PEM on boundary %d, must have node# on first entry\n", bci);
            return FAILURE;    /* must have se[0] be the node for a PEM */
        }
    }
    /* expect all BC to be used */
    for(int bc = 0; bc < n_bc; bc++) {
        if(type[bc] == UNKNOWN) {
            printf("error: unused boundary %d\n", bc);
            return FAILURE;
        }
    }
    return SUCCESS;
}

int check_model(model const * const mdl)
{
    return check_bc(&(mdl->fwd));
}

int calc_sys_size(model const * const m)
{
    const int dim = m->fwd.dim;
    const int n_bc = max_bc(&(m->fwd)) + 1;
    enum se_type type[n_bc];
    for(int i = 0; i < n_bc; i++) {
        type[i] = UNKNOWN;
    }
    type[0] = NONE; /* skip bc=0, assuming its the exterior */
    for(int i = 0; i < m->fwd.n_se; i ++ ) {
        const int bci = m->fwd.bc[i];
        if (type[bci] != UNKNOWN) {
            continue;
        }
        type[bci] = (m->fwd.surfaceelems[i * dim + 1] == 0) ? PEM : CEM;
    }
    int n_cem = 0;
    for(int bc = 0; bc < n_bc; bc++) {
        n_cem += (type[bc] == CEM) ? 1 : 0;
    }
    return (m->fwd.n_nodes + n_cem);
}

size_t calc_sys_nnz(model const * const m)
{
    const int dim = m->fwd.dim;
    const int n_bc = max_bc(&(m->fwd)) + 1;
    enum se_type type[n_bc];
    for(int i = 0; i < n_bc; i++) {
        type[i] = UNKNOWN;
    }
    int cnt[n_bc];
    for(int i = 0; i < n_bc; i++) {
        cnt[i] = 0;
    }
    for(int i = 0; i < m->fwd.n_se; i ++ ) {
        const int bci = m->fwd.bc[i];
        type[bci] = (m->fwd.surfaceelems[i * dim + 1] == 0) ? PEM : CEM;
        cnt[bci]++;
    }
    type[0] = NONE; /* force skipping bc=0, assuming its the exterior */
    /* count surfaceelems that are CEM */
    int se_cem = 0;
    for(int bc = 0; bc < n_bc; bc++) {
        se_cem += (type[bc] == CEM) ? cnt[bc] : 0;
    }
    const int ne = m->fwd.n_elems;
    const int se_n = calc_sys_elem_n(dim); /* sparse matrix entries per mesh element */
    const int nnz_per_se = calc_sys_cem_n(dim);
    return (se_n * ne) + (se_cem * nnz_per_se);
}
