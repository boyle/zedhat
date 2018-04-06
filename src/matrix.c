/* Copyright 2017, Alistair Boyle, 3-clause BSD License */

#include <stdlib.h> /* free */
#include <string.h> /* strncpy */
#include <matio.h>

#include "config.h"
#include "matrix.h"

int matrix_load(const char * file, matrix_t * matrix)
{
    mat_t * in = NULL;
    matvar_t * t;
    matrix->scale = 1.0;
    in = Mat_Open(file, MAT_ACC_RDONLY);
    if (in == NULL) {
        return 1; /* bad input filename */
    }
    t = Mat_VarReadNextInfo(in);
    if (t == NULL) {
        return 2;    /* end of list */
    }
    t = Mat_VarRead(in, t->name);
    Mat_VarPrint(t, 1);
    Mat_VarFree(t);
    Mat_Close(in);
    return 1;
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

void matrix_free(matrix_t * matrix)
{
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
    case IDENTITY: break;
    }
    matrix->type = IDENTITY;
    matrix->scale = 1.0;
    matrix->m = 0;
    matrix->n = 0;
    matrix->transposed = false;
    matrix->symmetric = false;
    matrix->units = "SI units?";
    matrix->name = "name?";
    strncpy(matrix->symbol, "□", 4);
    matrix->symbol[4] = '\0';
}
