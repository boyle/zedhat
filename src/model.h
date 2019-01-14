/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#ifndef __MODEL_H__
#define __MODEL_H__

typedef struct mesh {
    int dim;
    double * nodes; /* n_nodes x dim */
    int n_nodes;
    int * elems; /* n_elems x (dim+1) */
    int * matidx; /* n_elems x 1 */
    int n_elems;
    int * surfaceelems; /* n_se x dim */
    int * bc; /* n_se x 1 */
    int n_se;
} mesh_t;

void mesh_init(mesh_t * m);
void mesh_free(mesh_t * m);

double det(int n, double A[n][n]);
double * inv(int n, double A[n][n]);
double * calc_Se_v(const int nd, double const * const * const nodes, double * Se);
void calc_Se_ij(const int nd, int const * const elem, int * ii, int * jj);
int calc_Se_n(const int nd);
int calc_Se(const int nd, const int n_elems, double const * const nodes, int const * elems, int * ii, int * jj, double * Se);
#endif
