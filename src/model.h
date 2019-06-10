/* Copyright 2018,2019, Alistair Boyle, 3-clause BSD License */
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

#define MAX_MODEL_FORMAT 1
typedef struct model_type {
    mesh fwd; /* forward model */
    int n_data [2]; /* number of data measurements */
    int n_params [2]; /* number of image frames */
    int n_stimmeas; /* number of stimmeas measurements*/
    double * data; /* m x n */
    double * params; /* p x k */
    int * stimmeas; /* n_stimmeas x 4: sp sn mp mn */
    double hp; /* hyperparameter */
} model;

model * malloc_model();
void * free_model(model * m);
/* these handle free-ing of any old data and then malloc of new mesh variables */
int set_model_data(model * m, int rows, int cols);
int set_model_params(model * m, int rows, int cols);
int set_model_stimmeas(model * m, int rows);
void set_model_hp(model * m, double hp);
void set_mesh_dim(mesh * m, int dim);
int set_mesh_nodes(mesh * m, int n_nodes);
int set_mesh_elems(mesh * m, int n_elems);
int set_mesh_surfaceelems(mesh * m, int n_se);

double det(int n, double A[n][n]);
double * inv(int n, double A[n][n]);
int calc_Se_n(const int nd);
int calc_Se(mesh const * const mdl, int * ii, int * jj, double * Se);
int calc_gnd(const int gnd, size_t * nnz, int * ii, int * jj, double * Se);
int calc_stim_neumann(mesh const * const m, double amp, int bc, int gnd, double * b);
#endif
