/* Copyright 2017-2019, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdlib.h> /* malloc, free */
#include <string.h> /* strncmp, strdup, bzero */
#include <assert.h> /* assert */
#include <stdio.h> /* printf */
#include <cholmod.h> /* cholmod_* */

#include "matrix.h"

#ifdef UNIT_TESTING
extern void * _mock_test_malloc(const size_t size, const char * file, const int line);
extern void * _mock_test_realloc(void * const ptr, const size_t size, const char * file, const int line);
extern void _test_free(void * const ptr, const char * file, const int line);
extern char * _test_strdup(const char * ptr, const char * file, const int line);
#define malloc(size) _mock_test_malloc(size, __FILE__, __LINE__)
#define realloc(ptr, size) _mock_test_realloc(ptr, size, __FILE__, __LINE__)
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

static void free_matrix_data(matrix * M)
{
    assert(M != NULL);
    assert(M->type < MAX_MATRIX_TYPE);
    switch (M->type) {
    // TODO case CSR_SYMMETRIC:
    case CSC_SYMMETRIC:
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
        M->x.dense = NULL;
        break;
    // TODO case CSR:
    case CSC:
    case COO:
        free(M->x.sparse.a);
        free(M->x.sparse.ia);
        free(M->x.sparse.ja);
        M->x.sparse.a = NULL;
        M->x.sparse.ia = NULL;
        M->x.sparse.ja = NULL;
        break;
    default: /* NOP */
        break;
    }
    M->type = IDENTITY;
}

matrix * free_matrix(matrix * M)
{
    if(M == NULL) {
        return NULL;
    }
    free_matrix_data(M);
    free(M->symbol);
    free(M->name);
    free(M->units);
    matrix_init(M);
    free(M);
    return NULL;
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
    free_matrix_data(M);
    M->scale = 1.0;
    M->type = type;
    M->m = rows;
    M->n = cols;
    assert(rows > 0);
    assert(cols > 0);
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
//      M->x.sparse.ia = malloc(sizeof(int) * (rows + 1));
//      M->x.sparse.ja = malloc(sizeof(int) * nnz);
//      break;
    case CSC_SYMMETRIC:
        M->symmetric = true;
        M->type--; /* strip _SYMMETRIC */
    case CSC:
        M->x.sparse.a = malloc(sizeof(double) * nnz);
        M->x.sparse.ia = malloc(sizeof(int) * nnz);
        M->x.sparse.ja = malloc(sizeof(int) * (cols + 1));
        M->x.sparse.na = nnz;
        M->x.sparse.nia = nnz;
        M->x.sparse.nja = (cols + 1);
        break;
    case COO_SYMMETRIC:
        M->symmetric = true;
        M->type--; /* strip _SYMMETRIC */
    case COO:
        M->x.sparse.a = malloc(sizeof(double) * nnz);
        M->x.sparse.ia = malloc(sizeof(int) * nnz);
        M->x.sparse.ja = malloc(sizeof(int) * nnz);
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

/*
void matrix_transpose(matrix * M)
{
    assert(M != NULL);
    M->transposed = ! M->transposed;
}
*/

void printf_matrix(matrix const * const A)
{
    assert(A != NULL);
    printf("matrix %s: %zux%zu\n", A->symbol, A->m, A->n);
    printf(" ");
    if(A->name != NULL) {
        printf(" %s", A->name);
    }
    if(A->units != NULL) {
        printf(" [%s]", A->units);
    }
    if(A->symmetric) {
        printf(" symmetric");
    }
    switch(A->type) {
    case IDENTITY: printf(" IDENTITY\n"); break;
    case DENSE: printf(" DENSE\n"); break;
    case COO_SYMMETRIC: printf(" COO_SYMMETRIC\n"); break;
    case COO: printf(" COO (nnz=%zu)\n", A->x.sparse.na); break;
    case CSC_SYMMETRIC: printf(" CSC_SYMMETRIC\n"); break;
    case CSC: printf(" CSC (nnz=%zu)\n", A->x.sparse.na); break;
    case MAX_MATRIX_TYPE: printf(" MAX\n"); break; /* LCOV_EXCL_LINE */
    }
    if(A->type == DENSE) {
        const int rows = A->m;
        const int cols = A->n;
        const int maxi = rows < 10 ? rows : 10;
        const int maxj = cols < 10 ? cols : 10;
        for(int i = 0; i < maxi; i++) {
            printf("  ");
            for(int j = 0; j < maxj; j++) {
                printf(" %14.6g", A->x.dense[j * (A->m) + i]);
            }
            printf("\n");
        }
    }
    else if(A->type == COO) {
        const int nnz = A->x.sparse.na;
        const int maxi = nnz; // < 10 ? nnz : 10;
        for(int i = 0; i < maxi; i++) {
            printf("  (%4d,%4d) = %14.6e\n", A->x.sparse.ia[i], A->x.sparse.ja[i], A->x.sparse.a[i]);
        }
    }
    else if(A->type == CSC) {
        const int nnz = A->x.sparse.na;
        const int cols = A->n;
        int maxi = nnz; // < 10 ? nnz : 10;
        printf(" column indices: ");
        for(int i = 0; i < cols + 1; i++) {
            printf(" %d", A->x.sparse.ja[i]);
        }
        printf("\n");
        for(int i = 0; i < maxi; i++) {
            int col;
            for(col = 0; i > A->x.sparse.ja[col]; col++); /* search for col in compressed list */
            printf("  (%4d,%4d) = %14.6e\n", A->x.sparse.ia[i], col, A->x.sparse.a[i]);
        }
    }
}

// /* This is an *in-place* COO to CSC convertor.
//  * We need a qsort_r/s() to handle three parallel arrays. The interfaces and
//  * capabilities vary accross platforms (linux/bsd/macos/win); we start with our
//  * own recursive implementation. This implementation follows the Lomuto
//  * partition scheme (Hoare is better, and there are more optimizations); see
//  * https://en.wikipedia.org/wiki/Quicksort */
// #define cmp_coo(a,b,ii,jj) ( (a[ii] != a[jj]) ? (a[ii] - a[jj]) : (b[ii] - b[jj]) )
// #define swap(type, a, i, j) do { const type tmp = a[i]; a[i] = a[j]; a[j] = tmp; } while(0)
// static int my_qsort_partition(matrix * A, int lo, int hi)
// {
//     /* sort by column (ja) then row (ia) */
//     int i = lo;
//     for(int j = lo; j < hi - 1; j++) {
//         if(cmp_coo(A->x.sparse.ja, A->x.sparse.ia, i, j) < 0) {
//             swap(int, A->x.sparse.ja, i, j);
//             swap(int, A->x.sparse.ia, i, j);
//             swap(double, A->x.sparse.a, i, j);
//             i++;
//         }
//     }
//     swap(int, A->x.sparse.ja, i, hi);
//     swap(int, A->x.sparse.ia, i, hi);
//     swap(double, A->x.sparse.a, i, hi);
//     return i;
// }
//
// static void my_qsort_recurse(matrix * A, int lo, int hi)
// {
//     if(lo >= hi) {
//         return;
//     }
//     int p = my_qsort_partition(A, lo, hi);
//     my_qsort_recurse(A, lo, p - 1);
//     my_qsort_recurse(A, p + 1, hi);
// }
//
// static void my_qsort(matrix * A)
// {
//     my_qsort_recurse(A, 0, A->x.sparse.na - 1);
// }
//
//
// int coo_to_csc(matrix * A)
// {
//     my_qsort(A);
//     /* merge duplicate indices, compact the array */
//     int k = 0;
//     int rowk = A->x.sparse.ia[k];
//     int colk = A->x.sparse.ja[k];
//     for(int i = 1; i < A->x.sparse.na; i++) {
//         const int rowi = A->x.sparse.ia[i];
//         const int coli = A->x.sparse.ja[i];
//         if((rowk == rowi) && (colk == coli)) {
//             A->x.sparse.a[k] += A->x.sparse.a[i];
//         }
//         else {
//             const double vali = A->x.sparse.a[i];
//             A->x.sparse.a[k] = vali;
//             A->x.sparse.ia[k] = rowi;
//             A->x.sparse.ja[k] = coli; /* ... destroyed next*/
//             k++;
//             rowk = rowi;
//             colk = coli;
//         }
//     }
//        /* enlarge the column array if necessary */
//     if(A->x.sparse.nja < A->n) {
//         int * ptr = realloc(A->x.sparse.ja, sizeof(int) * A->n + 1);
//         if(ptr == NULL) {
//             return FAILURE;
//         }
//         else {
//             A->x.sparse.ja = ptr;
//         }
//     }
//     /* construct column pointers */
//     /* ... this loop could be merged with above */
//     int colj = 0;
//     for(int i = 0; i < k; i++) {
//         const int coli = A->x.sparse.ja[i];
//         if(coli != colj) {
//             A->x.sparse.ja[++colj] = i;
//         }
//     }
//     A->x.sparse.ja[A->n] = k;
//     A->x.sparse.ja[0] = 0;
//     /* resize the compacted arrays */
//     A->x.sparse.na = k;
//     A->x.sparse.nia = k;
//     A->x.sparse.nja = A->n + 1;
//     double * x1 = realloc(A->x.sparse.a, sizeof(double) * A->x.sparse.na);
//     if(x1 != NULL) {
//         A->x.sparse.a = x1;
//     }
//     int * x2 = realloc(A->x.sparse.ia, sizeof(int) * A->x.sparse.nia);
//     if(x2 != NULL) {
//         A->x.sparse.ia = x2;
//     }
//     int * x3 = realloc(A->x.sparse.ja, sizeof(int) * A->x.sparse.nja);
//     if(x3 != NULL) {
//        A->x.sparse.ja = x3;
//     }
//     return SUCCESS;
// }

int coo_to_csc(matrix * A)
{
    assert(A != NULL);
    assert(A->type == COO);
    const int xtype = CHOLMOD_REAL;
    cholmod_common c;
    cholmod_start (&c);                /* start CHOLMOD */
    cholmod_triplet * T = cholmod_allocate_triplet(A->m, A->n, A->x.sparse.na, +1, xtype, &c);
    assert(T != NULL);
    { /* copy contents from A */
        T->nnz = A->x.sparse.na;
        double * xx = T->x;
        int * ii = T->i;
        int * jj = T->j;
        for(int i = 0; i < T->nnz; i++) {
            ii[i] = A->x.sparse.ia[i];
            jj[i] = A->x.sparse.ja[i];
            xx[i] = A->x.sparse.a[i];
        }
    }
    //cholmod_print_triplet (T, "T", &c);
    assert(cholmod_check_triplet(T, &c));
    cholmod_sparse * As = cholmod_triplet_to_sparse(T, T->nnz, &c);
    //cholmod_print_sparse (As, "A", &c);
    assert(cholmod_check_sparse(As, &c));
    { /* copy contents to A */
        const int cols = As->ncol;
        double * xx = As->x;
        int * ii = As->i;
        int * jj = As->p;
        const int nnz = jj[cols];
        void * ptr;
        {
            ptr = realloc(A->x.sparse.a, sizeof(double) * nnz);
            if(ptr == NULL) {
                return FAILURE;
            }
            A->x.sparse.na = nnz;
            A->x.sparse.a = ptr;
        }
        {
            ptr = realloc(A->x.sparse.ia, sizeof(int) * nnz);
            if(ptr == NULL) {
                return FAILURE;
            }
            A->x.sparse.nia = nnz;
            A->x.sparse.ia = ptr;
        }
        {
            ptr = realloc(A->x.sparse.ja, sizeof(int) * (cols + 1));
            if(ptr == NULL) {
                return FAILURE;
            }
            A->x.sparse.ja = ptr;
            A->x.sparse.nja = cols + 1;
        }
        for(int i = 0; i < A->x.sparse.na; i++) {
            A->x.sparse.a[i] = xx[i];
            A->x.sparse.ia[i] = ii[i];
        }
        for(int i = 0; i < A->x.sparse.nja; i++) {
            A->x.sparse.ja[i] = jj[i];
        }
        A->type = CSC;
    }
    cholmod_free_sparse(&As, &c);
    cholmod_free_triplet(&T, &c);
    cholmod_finish (&c);               /* finish CHOLMOD */
    return SUCCESS;
}
