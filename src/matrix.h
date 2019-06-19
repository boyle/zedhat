/* Copyright 2017, Alistair Boyle, 3-clause BSD License */
#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stdbool.h> /* bool */
#include "config.h"

enum matrix_type {IDENTITY, /* DIAGONAL, */ DENSE, /* DENSE_SYMMETRIC, */ CSR, CSR_SYMMETRIC, CSC, CSC_SYMMETRIC, COO, COO_SYMMETRIC, MAX_MATRIX_TYPE};
/* NOTE: *_SYMMETRIC values are immediately converted to their base type and set/check the bool symmetric flag */

typedef struct matrix_sparse_type {
    double * a;   /* data */
    int * ia; /* row index */
    int * ja; /* column index */
    size_t nia; /* length of ia */
    size_t nja; /* length of ja */
    size_t na; /* length of a */
} matrix_sparse;

typedef struct matrix_t {
    double scale;
    enum matrix_type type;
    size_t m; /* rows */
    size_t n; /* cols */
    // bool transposed;
    bool symmetric;
    char * symbol; // potentially a UTF-8 4-byte string, null-terminated
    char * name;
    char * units;
    union {
        double * dense;    /* dense array, vector or matrix; column-major order (x[1] is row 2, column 1) */
        matrix_sparse sparse; /* sparse matrix (CSR, CSC, COO) */
    } x;
} matrix;

matrix * malloc_matrix(void);
int malloc_matrix_name(matrix * M, const char * name, const char * symbol, const char * units);
int malloc_matrix_data(matrix * M, enum matrix_type type, const size_t rows, const size_t cols, const size_t nnz);
void free_matrix_data(matrix * M);
matrix * free_matrix(matrix * M);
//void matrix_transpose(matrix * M);
void printf_matrix(matrix const * const A);
int coo_to_csc(matrix * A);
int coo_to_csr(matrix * A);
int copy_matrix(matrix const * const src, matrix ** dest, const char * symbol);

#endif
