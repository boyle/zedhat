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
    unsigned int * ia; /* row index */
    unsigned int * ja; /* column index */
    size_t nia; /* length of ia */
    size_t nja; /* length of ja */
    size_t na; /* length of a */
} matrix_sparse_t;

typedef struct matrix {
    double scale;
    matrix_type_t type;
    size_t m; /* rows */
    size_t n; /* cols */
    bool transposed;
    bool symmetric;
    char * symbol; // potentially a UTF-8 4-byte string, null-terminated
    char * name;
    char * units;
    union {
        double * dense;    /* dense array, vector or matrix */
        matrix_sparse_t * sparse; /* sparse matrix (CSR, CSC, COO) */
    };
} matrix_t;

int matrix_load(const char * file, matrix_t * matrix);
/* int matrix_save(const char * file, const matrix_t matrix); */

matrix_t * matrix_malloc(const char * name, const char * symbol, const char * units);
void matrix_free(matrix_t * matrix);

#endif
