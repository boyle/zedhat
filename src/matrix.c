/* Copyright 2017, Alistair Boyle, 3-clause BSD License */

#include <stdio.h>
#include <unistd.h> /* for access() */
#include <string.h> /* for strcmp */
#include <matio.h>
#include <stdint.h> /* uint32_t */
#include <stdbool.h> /* bool */
#include "config.h"
/* testa.mat <- a
 * testb.mat <- b
 * testab.mat <- {a,b}
 * a =    1     2     3
 *        4     5     6
 * b = (1,1)        1
 *     (1,2)        2
 *     (2,2)        5 */

typedef enum matrix_type {IDENTITY, DENSE, CSR, CSC, COO} matrix_type_t;

typdef struct matrix_sparse {
    double * a;   /* data */
    uint32_t * ia; /* row index */
    uint32_t * ja; /* column index */
    uint32_t nnz; /* non-zero entries */
    uint32_t n; /* allocated space */
} matrix_sparse_t;

typedef struct matrix {
    double scale;
    matrix_type_t type;
    uint32_t m; /* rows */
    uint32_t n; /* cols */
    bool transposed;
    bool symmetric;
    char symbol[5]; // potentially a UTF-8 4-byte string, null-terminated
    const char * name;
    const char * units;
    union {
        double * dense;    /* dense array, vector or matrix */
        matrix_sparse_t sparse; /* sparse matrix (CSR, CSC, COO) */
    }
} matrix_t;

int matrix_load(const char * file, matrix_t * matrix)
{
    matrix->scale = 1.0;
    mat_t * in = NULL;
    in = Mat_Open(argv[1], MAT_ACC_RDONLY);
    if(in == NULL) {
        fprintf(stderr, "error: %s: failed to open\n", file);
        return 1; /* bad input filename */
    }
    matvar_t * t = Mat_VarReadNextInfo(in);
    if(t != NULL) {
        return 2;    /* end of list */
    }
    t = Mat_VarRead(in, t->name);
    Mat_VarPrint(t, 1);
    Mat_VarFree(t);
    Mat_Close(in);
}

int matrix_save(const char * file, const matrix_t matrix)
{
    mat_t * out = Mat_Create(file, "");
    if(out == NULL) {
        return 1;    /* bad output filename */
    }
    matvar_t * t = Mat_VarReadNextInfo(in);
    Mat_Close(out);
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
    }
    free(matrix->name);
    matrix->units = "SI units?";
    matrix->name = "?";
    matrix->symbol = "□";
    matrix->type = IDENTITY;
    matrix->scale = 1.0;
    matrix->m = 0;
    matrix->n = 0;
    matrix->transposed = FALSE;
    matrix->symmetric = FALSE;
}
