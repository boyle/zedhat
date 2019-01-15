/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#include <string.h> /* bzero */
#include <stdlib.h> /* free, qsort */
#include <lapacke.h> /* inv: dgetrf, dgetri */
#include <math.h> /* fabs */
#include <limits.h> /* MAX_INT */

#include "matrix.h"
#include "model.h"

/* references:
 * [1] A. Boyle, PhD Thesis, 2016, Geophysical Applications of Electrical Impedance Tomography, Carleton University
 */

void mesh_init(mesh * m)
{
    bzero(m, sizeof(mesh));
}

void mesh_free(mesh * m)
{
    free(m->nodes);
    free(m->elems);
    free(m->matidx);
    free(m->surfaceelems);
    free(m->bc);
    mesh_init(m);
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
    int i, j, k, s;
    for(i = 0; i < n; i++, sign *= -1.0) {
        /* construct sub-matrix B */
        double B[n - 1][n - 1];
        for(j = 1; j < n; j++) {
            for(k = s = 0; k < n; k++) {
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
    int i, j, k, m;
    double E[n][n]; /* from [1] eq (C.33) */
    for(i = 0; i < n; i++) {
        E[i][0] = 1;
        for(j = 1; j < n; j++) {
            E[i][j] = nodes[i][j - 1];
        }
    }
    if(inv(n, E) == NULL) {
        return NULL;
    }
    const double detE = fabs(det(n, E));
    int ndf = 1; /* = !nd = factorial(number of dimensions) */
    for( i = n - 1; i > 0; i-- ) {
        ndf *= i;
    }
    const double c = 1 / (ndf * detE);
    /* Se = c * E1^T E1, where E1 is E without the top row */
    m = 0;
    for( i = 0; i < n; i++ ) {
        for( j = i; j < n; j++ ) {
            Se[m] = 0;
            for( k = 1; k < n; k++) {
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
    int i, j;
    for( i = 0; i < n; i++ ) {
        for( j = i; j < n; j++ ) {
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
    int i, j;
    const int n = calc_Se_n(nd);
    for( i = 0; i < n_elems; i++) {
        calc_Se_ij(nd, elems, ii, jj);
        double const * node_list [4];
        for( j = 0; j < nd + 1; j ++) {
            node_list[j] = &(nodes[(elems[j] - 1) * nd]);
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
 * Side effect: sorts ii, jj, Se
 * Returns number of entries removed */
int cmp_int( const void * a, const void * b )
{
    return *(int *)a - *(int *)b;
}
int calc_gnd(const int gnd, int * nnz, int * ii, int * jj, double * Se)
{
    int ret = 0;
    int i;
    int idx [*nnz];
    for ( i = 0; i < *nnz; i++ ) {
        if ( ii[i] == gnd || jj[i] == gnd) {
            idx[i] = INT_MAX;
            ret++;
        }
        else {
            idx[i] = i;
            if ( ii[i] > gnd ) {
                ii[i]--;
            }
            if ( jj[i] > gnd ) {
                jj[i]--;
            }
        }
    }
    qsort(idx, *nnz, sizeof(int), &cmp_int);
    for(i = 0 ; i < *nnz - ret; i++) {
        const int x = idx[i];
        ii[i] = ii[x];
        jj[i] = jj[x];
        Se[i] = Se[x];
    }
    nnz -= ret;
    return ret;
}
