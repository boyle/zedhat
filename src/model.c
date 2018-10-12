/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#include <string.h> /* memset */
#include <stdlib.h> /* free */
#include <lapacke.h> /* inv: dgetrf, dgetri */

#include "model.h"

void mesh_init(struct mesh * m)
{
    memset(m, 0, sizeof(struct mesh));
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

/*
double* A inv(double * A, int n)
{
    int ipiv[n + 1] = {0};
    lapack_int ret;
    if(LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv)) {
        return NULL;
    }
    if(LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, A, n, ipiv)) {
        return NULL;
    }
    return A;
}
*/

#if 0
double det(int n, double A[3][3])
{
    int i;
    double d = 0;
    for (i = 0; i < 3; i++) {
        d += A[0][i] * (A[1][(i + 1) % 3] * A[2][(i + 2) % 3] - A[2][(i + 1) % 3] * A[1][(i + 2) % 3]);
    }
    return d;
}
#else
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
#endif
