/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#ifndef __MODEL_H__
#define __MODEL_H__

enum geomtype {DEFAULT};
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
    enum geomtype type;
} mesh;

typedef enum enum_map_method_type {MAP_FWD, MAP_FWD_MATIDX, MAP_REC, MAP_REC_MATIDX} enum_map_method;

typedef struct model_type {
    mesh fwd; /* forward model */
    mesh rec; /* reconstruction model (optional) */
    enum_map_method map_method;
    int p; /* number of parameters ~ fwd.n_elems or mapping from rec/matidx to fwd */
    int m; /* number of data measurements */
    int n; /* number of frames */
    double * data; /* m x n */
    double * param; /* p x n */
    double hp; /* hyperparameter */
} model;

void model_init(model * m);
void model_free(model * m);

void mesh_init(mesh * m);
void mesh_free(mesh * m);

double det(int n, double A[n][n]);
double * inv(int n, double A[n][n]);
int calc_Se_n(const int nd);
int calc_Se(mesh const * const mdl, int * ii, int * jj, double * Se);
int calc_gnd(const int gnd, int * nnz, int * ii, int * jj, double * Se);
int calc_stim_neumann(mesh const * const m, double amp, int bc, int gnd, double * b);
#endif
