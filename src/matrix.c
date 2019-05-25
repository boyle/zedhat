/* Copyright 2017--2018, Alistair Boyle, 3-clause BSD License */

#include <stdlib.h> /* malloc, free */
#include <string.h> /* strncmp, strdup, bzero */
#include <assert.h> /* assert */

#include "config.h"
#include "matrix.h"

void matrix_init(matrix * m)
{
    bzero(m, sizeof(matrix));
}

int m_strdup(char ** ptr, const char * str)
{
    free(*ptr);
    *ptr = NULL;
    if(str == NULL) {
        return 1;
    }
    *ptr = strdup(str);
    if(*ptr == NULL) {
        return 0;    /* fail */
    }
    return 1; /* success */
}

int matrix_name(matrix * M, const char * name, const char * symbol, const char * units)
{
    if ( M == NULL ) {
        return 0; /* abort */
    }
    int ret1 = m_strdup(&(M->name), name);
    int ret2 = m_strdup(&(M->symbol), symbol);
    int ret3 = m_strdup(&(M->units), units);
    return ret1 && ret2 && ret3;
}

matrix * matrix_malloc(enum matrix_type type, const size_t rows, const size_t cols, const size_t nnz)
{
    matrix * M = malloc(sizeof(matrix));
    if ( M == NULL ) {
        return NULL;
    }
    matrix_init(M);
    M->scale = 1.0;
    M->type = type;
    M->m = rows;
    M->n = cols;
    /* malloc data */
    switch(type) {
    case IDENTITY:
        assert(nnz == 0);
        assert(rows == cols);
        break;
    case DIAGONAL:
        assert(nnz == (rows < cols ? rows : cols));
        M->dense = malloc(sizeof(double) * nnz);
        break;
    case DENSE_SYMMETRIC:
        assert(cols == rows);
        assert(nnz == ((rows * (rows - 1)) / 2 + rows));
        M->symmetric = true;
        M->type -= 4; /* strip _SYMMETRIC */
        M->dense = malloc(sizeof(double) * nnz);
        break;
    case DENSE:
        assert(nnz == rows * cols);
        M->dense = malloc(sizeof(double) * nnz);
        break;
    case CSR_SYMMETRIC:
        M->symmetric = true;
        M->type -= 4; /* strip _SYMMETRIC */
    case CSR:
        M->sparse.a = malloc(sizeof(double) * nnz);
        M->sparse.ia = malloc(sizeof(unsigned int) * (rows + 1));
        M->sparse.ja = malloc(sizeof(unsigned int) * nnz);
        break;
    case CSC_SYMMETRIC:
        M->symmetric = true;
        M->type -= 4; /* strip _SYMMETRIC */
    case CSC:
        M->sparse.a = malloc(sizeof(double) * nnz);
        M->sparse.ia = malloc(sizeof(unsigned int) * nnz);
        M->sparse.ja = malloc(sizeof(unsigned int) * (cols + 1));
        break;
    case COO_SYMMETRIC:
        M->symmetric = true;
        M->type -= 4; /* strip _SYMMETRIC */
    case COO:
        M->sparse.a = malloc(sizeof(double) * nnz);
        M->sparse.ia = malloc(sizeof(unsigned int) * nnz);
        M->sparse.ja = malloc(sizeof(unsigned int) * nnz);
        break;
    case MAX_MATRIX_TYPE:
        assert(false); /* LCOV_EXCL_LINE */
        break;
    }
    /* handle any malloc failures */
    if((M->type == DENSE) || (M->type == DIAGONAL)) {
        if (M->dense == NULL) {
            matrix_free(M);
            return NULL;
        }
    }
    else { /* COO, CSC, CSR */
        if ((M->sparse.a == NULL) || (M->sparse.ia == NULL) || (M->sparse.ja == NULL)) {
            matrix_free(M);
            return NULL;
        }
    }
    return M;
}


void matrix_free(matrix * M)
{
    if(M == NULL) {
        return;
    }
    assert(M->type < MAX_MATRIX_TYPE);
    switch (M->type) {
    case CSR_SYMMETRIC:
    case CSC_SYMMETRIC:
    case COO_SYMMETRIC:
        M->type -= 4;
        break;
    default: /* NOP */
        break;
    }
    switch (M->type) {
    case DIAGONAL:
    case DENSE:
        free(M->dense);
        break;
    case CSR:
    case CSC:
    case COO:
        free(M->sparse.a);
        free(M->sparse.ia);
        free(M->sparse.ja);
        break;
    default:
        assert(false); /* LCOV_EXCL_LINE */
    }
    free(M->symbol);
    free(M->name);
    free(M->units);
    matrix_init(M);
    free(M);
}

void matrix_transpose(matrix * M)
{
    assert(M != NULL);
    M->transposed = ! M->transposed;
}
