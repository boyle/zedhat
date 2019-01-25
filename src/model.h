/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#ifndef __MODEL_H__
#define __MODEL_H__

typedef struct mesh_type {
    int dim;
    double * nodes; /* n_nodes x dim */
    int n_nodes;
    int * elems; /* n_elems x (dim+1) */
    int * matidx; /* n_elems x 1 */
    int n_elems;
    int * surfaceelems; /* n_se x dim */
    int * bc; /* n_se x 1 */
    int n_se;
} mesh;

void mesh_init(mesh * m);
void mesh_free(mesh * m);

double det(int n, double A[n][n]);
double * inv(int n, double A[n][n]);
double * calc_Se_v(const int nd, double const * const * const nodes, double * Se);
int calc_Se_n(const int nd);
int calc_Se(mesh const * const mdl, int * ii, int * jj, double * Se);
int calc_gnd(const int gnd, int * nnz, int * ii, int * jj, double * Se);
int calc_stim_neumann(mesh const * const m, double amp, int bc, int gnd, double * b);
#endif
