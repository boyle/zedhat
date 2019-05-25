/* Copyright 2017--2018, Alistair Boyle, 3-clause BSD License */

#include <stdlib.h> /* malloc, free */
#include <string.h> /* strncmp, memset, strdup */

#include "config.h"
#include "matrix.h"

void matrix_init(matrix * m)
{
    bzero(m, sizeof(matrix));
}

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
