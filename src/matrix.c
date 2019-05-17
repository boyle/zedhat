/* Copyright 2017--2018, Alistair Boyle, 3-clause BSD License */

#include <stdlib.h> /* malloc, free */
#include <string.h> /* strncmp, memset, strdup */
#include <matio.h>

#include "config.h"
#include "matrix.h"

int matrix_load(const char * file, matrix * M)
{
    int ret = 0;
    mat_t * in = NULL;
    matvar_t * t = NULL;
    if (M == NULL) {
        return 1;    /* bad ptr */
    }
    in = Mat_Open(file, MAT_ACC_RDONLY);
    if (in == NULL) {
        return 2;    /* bad input filename */
    }
    do {
        Mat_VarFree(t);
        t = Mat_VarReadNextInfo(in);
    } while ( (t != NULL) &&
              strncmp(t->name, M->name, strlen(M->name)) );
    if (t == NULL) {
        ret = 3;    /* end of list, M->name not found */
        goto _matrix_load_quit;
    }
    Mat_VarFree(t);
    Mat_Rewind(in);
    t = Mat_VarRead(in, M->name);
    if (t == NULL) {
        ret = 4;
        goto _matrix_load_quit;
    }
    /* Mat_VarPrint(t, 1); */
    /* now check that we like what we found */
    M->scale = 1.0;
    if (t->isComplex) {
        ret = 5; /* expect 'real' data */
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
    M->m = 0;
    M->n = 1;
    int i;
    for (i = 0; i < t->rank; i++) {
        switch(i) {
        case 0: M->m = t->dims[i]; break;
        case 1: M->n = t->dims[i]; break;
        default: ret = 6; goto _matrix_load_quit; /* too many dimensions */
        }
    }
    switch (t->class_type) {
    case MAT_C_DOUBLE: /* steal 'dense' matrix data */
        M->type = DENSE;
        M->dense = (double *) t->data;
        t->data = NULL;
        break;
    case MAT_C_SPARSE:
        /* TODO proper loading of sparse matrices */
        M->type = DENSE;
        M->dense = calloc(M->m * M->n, sizeof(double));
        if(M->dense == NULL) {
            ret = 7;
            goto _matrix_load_quit;
        }
        mat_sparse_t * ts = (mat_sparse_t *) t->data;
        if (ts == NULL) {
            ret = 8;
            goto _matrix_load_quit;
        }
        /* Compressed Sparse Column: col = j+1; row = k+1 */
        int i, j;
        for (j = 0; j < ts->njc - 1; j++) { /* matlab likes CSC */
            for (i = ts->jc[j]; i < ts->jc[j + 1] && i < ts->ndata; i++) {
                const int k = ts->ir[i];
                const double dd = ((double *) ts->data)[i];
                M->dense[ (j * M->m ) + k ] = dd;
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
int matrix_save(const char * file, const matrix * M)
{
}
*/

void matrix_sparse_free(matrix_sparse * sparse)
{
    free(sparse->a);
    free(sparse->ia);
    free(sparse->ja);
    sparse->a = NULL;
    sparse->ia = NULL;
    sparse->ja = NULL;
}

matrix * matrix_malloc(const char * name, const char * symbol, const char * units)
{
    matrix * M = malloc(sizeof(matrix));
    if ( M == NULL ) {
        return NULL;
    }
    memset(M, 0, sizeof(matrix));
    if(name) {
        M->name = strdup(name);
    }
    if(symbol) {
        M->units = strdup(symbol);
    }
    if(units) {
        M->symbol = strdup(units);
    }
    return M;
}

void matrix_free(matrix * M)
{
    if(M == NULL) {
        return;
    }
    switch (M->type) {
    case DENSE:
        free(M->dense);
        M->dense = NULL;
        break;
    case CSR:
    case CSC:
    case COO:
        matrix_sparse_free(M->sparse);
        M->sparse = NULL;
        break;
    case IDENTITY:
        break;
    }
    free(M->symbol);
    free(M->name);
    free(M->units);
    M->symbol = NULL;
    M->name = NULL;
    M->units = NULL;
    free(M);
}
