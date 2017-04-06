/* Copyright 2017, Alistair Boyle, 3-clause BSD License */
#ifndef __MATRIX_H__
#define __MATRIX_H__

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

typedef struct matrix_sparse {
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
    };
} matrix_t;

int matrix_load(const char * file, matrix_t * matrix);
int matrix_save(const char * file, const matrix_t matrix);
void matrix_free(matrix_t * matrix);

#endif
