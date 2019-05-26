/* Copyright 2017--2018, Alistair Boyle, 3-clause BSD License */

#include <stdlib.h> /* malloc, free */
#include <string.h> /* strncmp, strdup, bzero */
#include <assert.h> /* assert */

#include "config.h"
#include "matrix.h"

#ifdef UNIT_TESTING
extern void* _test_malloc(const size_t size, const char* file, const int line);
extern void _test_free(void * const ptr, const char * file, const int line);
#define malloc(size) _test_malloc(size, __FILE__, __LINE__)
#define free(ptr) _test_free(ptr, __FILE__, __LINE__)
#endif

#define SUCCESS 1
#define FAILURE 0

void matrix_init(matrix * m)
{
    bzero(m, sizeof(matrix));
}

matrix* malloc_matrix() {
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

int malloc_matrix_data(matrix*M, enum matrix_type type, const size_t rows, const size_t cols, const size_t nnz)
{
   if(M== NULL)
      return FAILURE;
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
//  case DIAGONAL:
//      assert(nnz == (rows < cols ? rows : cols));
//      M->dense = malloc(sizeof(double) * nnz);
//      break;
//  case DENSE_SYMMETRIC:
//      assert(cols == rows);
//      assert(nnz == ((rows * (rows - 1)) / 2 + rows));
//      M->symmetric = true;
//      M->type -= 4; /* strip _SYMMETRIC */
//      M->dense = malloc(sizeof(double) * nnz);
//      break;
    case DENSE:
        assert(nnz == rows * cols);
        M->dense = malloc(sizeof(double) * nnz);
        break;
//  case CSR_SYMMETRIC:
//      M->symmetric = true;
//      M->type -= 4; /* strip _SYMMETRIC */
//  case CSR:
//      M->sparse.a = malloc(sizeof(double) * nnz);
//      M->sparse.ia = malloc(sizeof(unsigned int) * (rows + 1));
//      M->sparse.ja = malloc(sizeof(unsigned int) * nnz);
//      break;
//  case CSC_SYMMETRIC:
//      M->symmetric = true;
//      M->type -= 4; /* strip _SYMMETRIC */
//  case CSC:
//      M->sparse.a = malloc(sizeof(double) * nnz);
//      M->sparse.ia = malloc(sizeof(unsigned int) * nnz);
//      M->sparse.ja = malloc(sizeof(unsigned int) * (cols + 1));
//      break;
    case COO_SYMMETRIC:
        M->symmetric = true;
        M->type -= 4; /* strip _SYMMETRIC */
    case COO:
        M->sparse.a = malloc(sizeof(double) * nnz);
        M->sparse.ia = malloc(sizeof(unsigned int) * nnz);
        M->sparse.ja = malloc(sizeof(unsigned int) * nnz);
        M->sparse.na = nnz;
        M->sparse.nia = nnz;
        M->sparse.nja = nnz;
        break;
    case MAX_MATRIX_TYPE:
        assert(false); /* LCOV_EXCL_LINE */
        break;
    }
    /* handle any malloc failures */
    if((M->type == DENSE) /*|| (M->type == DIAGONAL)*/) {
        if (M->dense == NULL) {
            return FAILURE;
        }
    }
    else { /* COO, CSC, CSR */
        if ((M->sparse.a == NULL) || (M->sparse.ia == NULL) || (M->sparse.ja == NULL)) {
            return FAILURE;
        }
    }
    return SUCCESS;
}


void free_matrix(matrix * M)
{
    if(M == NULL) {
        return;
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
        free(M->dense);
        break;
    // TODO case CSR:
    // TODO case CSC:
    case COO:
        free(M->sparse.a);
        free(M->sparse.ia);
        free(M->sparse.ja);
        break;
    default: /* NOP */
        break;
    }
    free(M->symbol);
    free(M->name);
    free(M->units);
    matrix_init(M);
    free(M);
}

/*
void matrix_transpose(matrix * M)
{
    assert(M != NULL);
    M->transposed = ! M->transposed;
}
*/
