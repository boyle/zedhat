/* Copyright 2017--2018, Alistair Boyle, 3-clause BSD License */

#include <stdlib.h> /* malloc, free */
#include <string.h> /* strncpy, strncmp, memset */
#include <matio.h>

#include "config.h"
#include "matrix.h"

int matrix_load(const char * file, matrix_t * matrix)
{
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
    Mat_VarPrint(t, 1);
    Mat_VarFree(t);
    Mat_Close(in);
    return 0;
}

int matrix_save(const char * file, const matrix_t matrix)
{
    //mat_t * in = NULL;
    mat_t * out = Mat_Create(file, "");
    if(out == NULL) {
        return 1;    /* bad output filename */
    }
    //matvar_t * t = Mat_VarReadNextInfo(in);
    Mat_Close(out);
    return 1;
    return 0;
}

void matrix_sparse_free(matrix_sparse_t * sparse)
{
    free(sparse->a);
    free(sparse->ia);
    free(sparse->ja);
    sparse->a = NULL;
    sparse->ia = NULL;
    sparse->ja = NULL;
    sparse->nnz = 0;
    sparse->n = 0;
}

matrix_t * matrix_malloc(const char * name)
{
    matrix_t * matrix = malloc(sizeof(matrix_t));
    if ( matrix == NULL ) {
        return matrix;
    }
    memset(matrix, 0, sizeof(matrix_t));
    matrix->name = name;
    matrix->symbol =  "□";
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
        matrix_sparse_free(&(matrix->sparse));
        break;
    case IDENTITY:
        break;
    }
    matrix->m = 0;
    matrix->n = 0;
    matrix->symbol =  "□";
    free(matrix);
}
