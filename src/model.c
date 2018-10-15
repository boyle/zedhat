/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#include <string.h> /* bzero */
#include <stdlib.h> /* free */
#include <lapacke.h> /* inv: dgetrf, dgetri */

#include "model.h"

/* references:
 * [1] A. Boyle, PhD Thesis, 2016, Geophysical Applications of Electrical Impedance Tomography, Carleton University
 */

void mesh_init(struct mesh * m)
{
    bzero(m, sizeof(struct mesh));
}

void mesh_free(struct mesh * m)
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
