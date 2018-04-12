/* Copyright 2017--2018, Alistair Boyle, 3-clause BSD License */

#include <stdlib.h> /* malloc, free */
#include <string.h> /* strncmp, memset, strdup */
#include <matio.h>

#include "config.h"
#include "matrix.h"

int matrix_load(const char * file, matrix_t * matrix)
{
    int ret = 0;
    mat_t * in = NULL;
    matvar_t * t;
    if (matrix == NULL) {
        return 1;    /* bad ptr */
    }
    in = Mat_Open(file, MAT_ACC_RDONLY);
    if (in == NULL) {
        return 2;    /* bad input filename */
    }
    do {
        t = Mat_VarReadNextInfo(in);
    }
    while ( (t != NULL) &&
            strncmp(t->name, matrix->name, strlen(matrix->name)) );
    if (t == NULL) {
        return 3;    /* end of list, matrix->name not found */
    }
    t = Mat_VarRead(in, t->name);
    /* Mat_VarPrint(t, 1); */
    /* now check that we like what we found */
    matrix->scale = 1.0;
    if (t->isComplex) {
        ret = 4; /* expect 'real' data */
        goto _matrix_load_quit;
    }
    if (t->isLogical) {
        ret = 5; /* expect 'real' data */
        goto _matrix_load_quit;
    }
    if (t->data_type != MAT_T_DOUBLE) {
        ret = 6; /* expect 'double' data */
        goto _matrix_load_quit;
    }
    matrix->m = 0;
    matrix->n = 1;
    for (int i = 0; i < t->rank; i++) {
        switch(i) {
        case 0: matrix->m = t->dims[i]; break;
        case 1: matrix->n = t->dims[i]; break;
        default: ret = 6; goto _matrix_load_quit; /* too many dimensions */
        }
    }
    switch (t->class_type) {
    case MAT_C_DOUBLE: /* steal 'dense' matrix data */
        matrix->type = DENSE;
        matrix->dense = (double *) t->data;
        t->data = NULL;
        break;
    case MAT_C_SPARSE:
        /* TODO proper loading of sparse matrices */
        matrix->type = DENSE;
        matrix->dense = calloc(matrix->m * matrix->n, sizeof(double));
        if(matrix->dense == NULL) {
            ret = 7;
            goto _matrix_load_quit;
        }
        mat_sparse_t * ts = (mat_sparse_t *) t->data;
        /* Compressed Sparse Column: col = j+1; row = k+1 */
        for (int j = 0; j < ts->njc - 1; j++) { /* matlab likes CSC */
            for (int i = ts->jc[j]; i < ts->jc[j + 1] && i < ts->ndata; i++) {
                const int k = ts->ir[i];
                const double dd = ((double *) ts->data)[i];
                matrix->dense[ (j * matrix->m ) + k ] = dd;
            }
        }
        break;
    default:
        ret = 7; /* expect dense or sparse matrices */
        goto _matrix_load_quit;
    }
_matrix_load_quit: /* and clean up */
    Mat_VarFree(t);
    Mat_Close(in);
    return ret;
}

/*
int matrix_save(const char * file, const matrix_t matrix)
{
}
*/

void matrix_sparse_free(matrix_sparse_t * sparse)
{
    free(sparse->a);
    free(sparse->ia);
    free(sparse->ja);
    sparse->a = NULL;
    sparse->ia = NULL;
    sparse->ja = NULL;
}

matrix_t * matrix_malloc(const char * name)
{
    matrix_t * matrix = malloc(sizeof(matrix_t));
    if ( matrix == NULL ) {
        return NULL;
    }
    memset(matrix, 0, sizeof(matrix_t));
    matrix->name = strdup(name);
    matrix->units = strdup("□");
    matrix->symbol = strdup("□");
    return matrix;
}

void matrix_free(matrix_t * matrix)
{
    if(matrix == NULL) {
        return;
    }
    switch (matrix->type) {
    case DENSE:
        free(matrix->dense);
        matrix->dense = NULL;
        break;
    case CSR:
    case CSC:
    case COO:
        matrix_sparse_free(matrix->sparse);
        matrix->sparse = NULL;
        break;
    case IDENTITY:
        break;
    }
    free(matrix->symbol);
    free(matrix->name);
    free(matrix->units);
    matrix->symbol = NULL;
    matrix->name = NULL;
    matrix->units = NULL;
    free(matrix);
}
