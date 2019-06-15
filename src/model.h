/* Copyright 2018,2019, Alistair Boyle, 3-clause BSD License */
#ifndef __MODEL_H__
#define __MODEL_H__

#include <stddef.h> /* size_t */

typedef struct mesh_type {
    int dim; /* dimensionality: 2 = 2D, 3 = 3D */
    /* nodes */
    double * nodes; /* size: n_nodes x dim */
    int n_nodes;
    /* volume elements */
    int * elems; /* size: n_elems x (dim+1) */
    int * matidx; /* size: n_elems x 1 */
    int n_elems;
    /* surface elements */
    int * surfaceelems; /* size: n_se x dim */
    int * bc; /* size: n_se x 1 */
    int n_se;
    /* parameter map
     * optional COO sparse matrix describing mapping from parameters to elements */
    int * pmap_param; /* size: n_pmap x 1 */
    int * pmap_elem; /* size: n_pmap x 1 */
    double * pmap_frac; /* size: n_pmap x 1 */
    int n_pmap;
} mesh;

typedef struct model_type {
    mesh fwd; /* forward model */
    mesh rec; /* reconstruction model */
    double hp; /* hyperparameter */
    /* electrode# to system matrix row/col# */
    int * elec_to_sys;
    int n_elec; /* number of electrodes */
    /* measurement data, n frames */
    double * data; /* m x n */
    int n_data [2]; /* number of data measurements m */
    /* parameters (e.g. conductivity) */
    double * params; /* p x k */
    int n_params [2]; /* number of image frames */
    /* stimulus and measurement sequence */
    int * stimmeas; /* n_stimmeas x 4: sp sn mp mn */
    double * measgain; /* n_stimmeas x 1: measurement gain assuming unit (1 Amp) stimulus */
    int n_stimmeas; /* number of stimmeas measurements*/
    /* contact impedances */
    int * zcbc;
    double * zc;
    int n_zc;
} model;

model * malloc_model();
void * free_model(model * m);
void reset_model(model * m);
/* these handle free-ing of any old data and then malloc of new mesh variables */
int set_model_data(model * m, int rows, int cols);
int set_model_params(model * m, int rows, int cols);
int set_model_stimmeas(model * m, int rows);
void set_mesh_dim(mesh * m, int dim);
int set_mesh_nodes(mesh * m, int n_nodes);
int set_mesh_elems(mesh * m, int n_elems);
int set_mesh_surfaceelems(mesh * m, int n_se);
int set_mesh_pmap(mesh * m, int n_pmap);
int set_model_zc(model * m, int n_zc);

int calc_elec_zc_map(model * m);
int calc_elec_to_sys_map(model * m);

double det(int n, double A[n][n]);
double * inv(int n, double A[n][n]);
int calc_sys_size(model const * const m);
size_t calc_sys_nnz(model const * const m);
int calc_sys_elem_n(const int nd);
int calc_sys_cem_n(const int nd);
int calc_sys_elem(mesh const * const mdl, int * ii, int * jj, double * Se);
int calc_sys_cem(mesh const * const m, int bc, double zc, int cem_node, int Se_idx, int * ii, int * jj, double * Se);
void calc_sys_gnd(const int gnd, size_t nnz, int * ii, int * jj, double * Se);
int calc_stim_neumann(mesh const * const m, double amp, int bc, double * b);
void calc_stim_gnd(model const * const m, int gnd, double * b);
void calc_stim(model const * const m, int idx, double * b);
double calc_meas(model const * const m, int idx, double * x);

int check_model(model const * const m);
#endif
