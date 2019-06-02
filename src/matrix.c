/* Copyright 2017-2019, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdlib.h> /* malloc, free */
#include <string.h> /* strncmp, strdup, bzero */
#include <assert.h> /* assert */

#include "matrix.h"

#ifdef UNIT_TESTING
extern void * _mock_test_malloc(const size_t size, const char * file, const int line);
extern void _test_free(void * const ptr, const char * file, const int line);
extern char * _test_strdup(const char * ptr, const char * file, const int line);
#define malloc(size) _mock_test_malloc(size, __FILE__, __LINE__)
#define free(ptr) _test_free(ptr, __FILE__, __LINE__)
#undef strdup
#define strdup(ptr) _test_strdup(ptr, __FILE__, __LINE__)
#endif

#define SUCCESS 1
#define FAILURE 0

void matrix_init(matrix * m)
{
    bzero(m, sizeof(matrix));
}

matrix * malloc_matrix()
{
    matrix * M = malloc(sizeof(matrix));
    if ( M == NULL ) {
        return NULL;
    }
    matrix_init(M);
    return M;
}

int m_strdup(char ** ptr, const char * str)
{
    free(*ptr);
    *ptr = NULL;
    if(str == NULL) {
        return SUCCESS;
    }
    *ptr = strdup(str);
    if(*ptr == NULL) {
        return FAILURE;
    }
    return SUCCESS;
}

int malloc_matrix_name(matrix * M, const char * name, const char * symbol, const char * units)
{
    if ( M == NULL ) {
        return FAILURE;
    }
    int ret1 = m_strdup(&(M->name), name);
    int ret2 = m_strdup(&(M->symbol), symbol);
    int ret3 = m_strdup(&(M->units), units);
    return ret1 && ret2 && ret3;
}

int malloc_matrix_data(matrix * M, enum matrix_type type, const size_t rows, const size_t cols, const size_t nnz)
{
    if(M == NULL) {
        return FAILURE;
    }
    M->scale = 1.0;
    M->type = type;
    M->m = rows;
    M->n = cols;
    assert(rows > 0);
    assert(cols > 0);
    assert(nnz <= cols * rows);
    /* malloc data */
    assert(type < MAX_MATRIX_TYPE);
    switch(type) {
    case IDENTITY:
        assert(nnz == 0);
        assert(rows == cols);
        break;
//  case DIAGONAL:
//      assert(nnz == (rows < cols ? rows : cols));
//      M->x.dense = malloc(sizeof(double) * nnz);
//      break;
//  case DENSE_SYMMETRIC:
//      assert(cols == rows);
//      assert(nnz == ((rows * (rows - 1)) / 2 + rows));
//      M->symmetric = true;
//      M->type--; /* strip _SYMMETRIC */
//      M->x.dense = malloc(sizeof(double) * nnz);
//      break;
    case DENSE:
        assert(nnz == rows * cols);
        M->x.dense = malloc(sizeof(double) * nnz);
        break;
//  case CSR_SYMMETRIC:
//      M->symmetric = true;
//      M->type--; /* strip _SYMMETRIC */
//  case CSR:
//      M->x.sparse.a = malloc(sizeof(double) * nnz);
//      M->x.sparse.ia = malloc(sizeof(unsigned int) * (rows + 1));
//      M->x.sparse.ja = malloc(sizeof(unsigned int) * nnz);
//      break;
//  case CSC_SYMMETRIC:
//      M->symmetric = true;
//      M->type--; /* strip _SYMMETRIC */
//  case CSC:
//      M->x.sparse.a = malloc(sizeof(double) * nnz);
//      M->x.sparse.ia = malloc(sizeof(unsigned int) * nnz);
//      M->x.sparse.ja = malloc(sizeof(unsigned int) * (cols + 1));
//      break;
    case COO_SYMMETRIC:
        M->symmetric = true;
        M->type--; /* strip _SYMMETRIC */
    case COO:
        M->x.sparse.a = malloc(sizeof(double) * nnz);
        M->x.sparse.ia = malloc(sizeof(unsigned int) * nnz);
        M->x.sparse.ja = malloc(sizeof(unsigned int) * nnz);
        M->x.sparse.na = nnz;
        M->x.sparse.nia = nnz;
        M->x.sparse.nja = nnz;
        break;
    case MAX_MATRIX_TYPE: /* LCOV_EXCL_LINE */
        break; /* LCOV_EXCL_LINE */
    }
    /* handle any malloc failures */
    if(/*(*/M->type == DENSE/*) || (M->type == DIAGONAL)*/) {
        if (M->x.dense == NULL) {
            return FAILURE;
        }
    }
    else if (M->type != IDENTITY) { /* COO, CSC, CSR */
        if ((M->x.sparse.a == NULL) || (M->x.sparse.ia == NULL) || (M->x.sparse.ja == NULL)) {
            return FAILURE;
        }
    }
    return SUCCESS;
}


matrix * free_matrix(matrix * M)
{
    if(M == NULL) {
        return NULL;
    }
    assert(M->type < MAX_MATRIX_TYPE);
    switch (M->type) {
    // TODO case CSR_SYMMETRIC:
    // TODO case CSC_SYMMETRIC:
    case COO_SYMMETRIC:
        M->type--;
        break;
    default: /* NOP */
        break;
    }
    switch (M->type) {
    // TODO case DIAGONAL:
    case DENSE:
        free(M->x.dense);
        break;
    // TODO case CSR:
    // TODO case CSC:
    case COO:
        free(M->x.sparse.a);
        free(M->x.sparse.ia);
        free(M->x.sparse.ja);
        break;
    default: /* NOP */
        break;
    }
    free(M->symbol);
    free(M->name);
    free(M->units);
    matrix_init(M);
    free(M);
    return NULL;
}

/*
void matrix_transpose(matrix * M)
{
    assert(M != NULL);
    M->transposed = ! M->transposed;
}
*/
