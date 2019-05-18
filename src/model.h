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


typedef struct model_type {
    mesh fwd; /* forward model */
    int n_data [2]; /* number of data measurements */
    int n_params [2]; /* number of image frames */
    int n_stimmeas; /* number of stimmeas measurements*/
    double * data; /* m x n */
    double * params; /* p x k */
    int * stimmeas; /* n_stimmeas x 4: sp sn mp mn */
    double hp; /* hyperparameter */
    int format; /* zedhat file format */
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
