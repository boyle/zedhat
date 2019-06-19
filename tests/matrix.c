/* Copyright 2019, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdarg.h> /* cmocka.h */
#include <stddef.h> /* cmocka.h */
#include <setjmp.h> /* cmocka.h */
#include <stdint.h> /* cmocka.h */

#include <string.h> /* strlen + strcpy = mock strdup*/
#include <stdio.h> /* printf */

#include "cmocka.h"
#include "matrix.h"

void * _mock_test_malloc(const size_t size, const char * file, const int line)
{
    if(mock()) {
        printf("  malloc x\n");
        return NULL;
    }
    else {
        printf("  malloc ✓\n");
        return _test_malloc(size, file, line);
    }
}

void * _mock_test_realloc(void * const ptr, const size_t size, const char * file, const int line)
{
    if(mock()) {
        printf("  realloc x\n");
        return NULL;
    }
    else {
        printf("  realloc ✓\n");
        return _test_realloc(ptr, size, file, line);
    }
}

char * _test_strdup(const char * str, const char * file, const int line)
{
    size_t size = strlen(str) + 1;
    char * n = _mock_test_malloc(size, file, line);
    if(n != NULL) {
        strcpy(n, str);
    }
    return n;
}

void test_malloc_matrix(void ** state)
{
    matrix * M = NULL;
    will_return(_mock_test_malloc, 0);
    M = malloc_matrix();
    assert_non_null(M); /* success */
    M = free_matrix(M);
    will_return(_mock_test_malloc, 1);
    M = malloc_matrix();
    assert_null(M); /* failure */
    M = free_matrix(M);
}

void test_malloc_matrix_name_null(void ** state)
{
    int ret = malloc_matrix_name(NULL, "A", "B", "C");
    assert_int_equal(ret, 0);
}

int bitcnt(int x)
{
    int sum = 0;
    while(x != 0) {
        sum += x & 1;
        x = x >> 1;
    }
    return sum;
}

void test_malloc_matrix_name_happy(void ** state)
{
    will_return(_mock_test_malloc, 0);
    matrix * M = malloc_matrix();
    assert_non_null(M);
    for(int i = 0; i < 8; i++) {
        if(bitcnt(i) > 0) {
            will_return_count(_mock_test_malloc, 0, bitcnt(i));
        }
        int ret = malloc_matrix_name(M, i & 1 ? "A" : NULL, i & 2 ? "B" : NULL, i & 4 ? "C" : NULL);
        printf("%d='b%d%d%d (%d)", i, (i >> 0) & 1, (i >> 1) & 1, (i >> 2) & 1, bitcnt(i));
        printf_matrix(M);
        assert_int_equal(ret, 1);
    }
    M = free_matrix(M);
}

void test_malloc_matrix_name_sad(void ** state)
{
    will_return(_mock_test_malloc, 0);
    matrix * M = malloc_matrix();
    assert_non_null(M);
    for(int i = 1; i < 8; i++) {
        if(bitcnt(i) > 1) {
            will_return_count(_mock_test_malloc, 0, bitcnt(i) - 1);
        }
        will_return(_mock_test_malloc, 1);
        int ret = malloc_matrix_name(M, i & 1 ? "A" : NULL, i & 2 ? "B" : NULL, i & 4 ? "C" : NULL);
        assert_int_equal(ret, 0);
    }
    M = free_matrix(M);
}

void test_malloc_matrix_data_null(void ** state)
{
    int ret = malloc_matrix_data(NULL, 0, 0, 0, 0);
    assert_int_equal(ret, 0);
}

void test_malloc_matrix_data_identity(void ** state)
{
    will_return(_mock_test_malloc, 0);
    matrix * M = malloc_matrix();
    assert_non_null(M);
    int ret = malloc_matrix_data(M, IDENTITY, 1, 1, 0);
    assert_int_equal(ret, 1);
    matrix * C = NULL;
    will_return(_mock_test_malloc, 0);
    assert_int_equal(copy_matrix(M, &C, NULL), 1);
    C = free_matrix(C);
    M = free_matrix(M);
}

void test_malloc_matrix_data_dense(void ** state)
{
    int ret;
    matrix * M = NULL;
    /* bad data malloc */
    will_return(_mock_test_malloc, 0);
    M = malloc_matrix();
    will_return(_mock_test_malloc, 1);
    ret = malloc_matrix_data(M, DENSE, 2, 2, 4);
    assert_int_equal(ret, 0);
    /* good data malloc */
    will_return(_mock_test_malloc, 0);
    ret = malloc_matrix_data(M, DENSE, 2, 2, 4);
    for(int i = 0; i < M->m * M->n; i++) {
        M->x.dense[i] = i;
    }
    printf_matrix(M);
    assert_int_equal(ret, 1);
    /* good copy */
    printf("copy1 - good malloc\n");
    matrix * C = NULL;
    printf_matrix(C); /* should be okay printing from a NULL ptr */
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 0);
    assert_int_equal(copy_matrix(M, &C, NULL), 1);
    printf_matrix(C);
    /* bad copy */
    printf("copy2 - bad malloc\n");
    will_return(_mock_test_malloc, 0);
    will_return(_mock_test_malloc, 1);
    assert_int_equal(copy_matrix(M, &C, NULL), 0);
    printf_matrix(C);
    C = free_matrix(C);
    M = free_matrix(M);
}

void core_malloc_matrix_data(enum matrix_type type)
{
    for(int i = 0; i < 6; i++) {
        printf("i=%d\n", i);
        will_return(_mock_test_malloc, 0);
        matrix * M = malloc_matrix();
        if(i > 0) {
            will_return_count(_mock_test_malloc, 0, (i < 4) ? i : 3);
        }
        if(i < 3) {
            will_return(_mock_test_malloc, 1);
        }
        int ret = malloc_matrix_data(M, type, 2, 3, 2);
        assert_int_equal(ret, i < 3 ? 0 : 1);
        if(ret) {
            for(int i = 0; i < M->x.sparse.na; i++) {
                M->x.sparse.a[i] = i * 10;
            }
            for(int i = 0; i < M->x.sparse.nia; i++) {
                M->x.sparse.ia[i] = i % 2;
            }
            for(int i = 0; i < M->x.sparse.nja; i++) {
                M->x.sparse.ja[i] = (3 - i) % 3;
            }
            printf_matrix(M);
            assert_int_equal(M->x.sparse.na, 2);
            if(type == COO) {
                assert_int_equal(M->x.sparse.nia, 2);
                assert_int_equal(M->x.sparse.nja, 2);
            }
            else if(type == CSC) {
                assert_int_equal(M->x.sparse.nia, 2);
                assert_int_equal(M->x.sparse.nja, 4);
            }
            else if(type == CSR) {
                assert_int_equal(M->x.sparse.nia, 3);
                assert_int_equal(M->x.sparse.nja, 2);
            }
            else {
                fail();
            }
        }
        if(i > 0) {
            will_return_count(_mock_test_malloc, 0, i);
        }
        if(i < 5) {
            will_return(_mock_test_malloc, 1);
        }
        matrix * C = NULL;
        ret = copy_matrix(M, &C, "COPY");
        assert_int_equal(ret, i < 5 ? 0 : 1);
        if(ret)  {
            printf_matrix(C);
        }
        M = free_matrix(M);
        C = free_matrix(C);
    }
}

void test_malloc_matrix_data_coo(void ** state)
{
    core_malloc_matrix_data(COO);
}

void test_malloc_matrix_data_csc(void ** state)
{
    core_malloc_matrix_data(CSC);
}

void test_malloc_matrix_data_csr(void ** state)
{
    core_malloc_matrix_data(CSR);
}

void core_malloc_matrix_data_symmetric(enum matrix_type type)
{
    will_return(_mock_test_malloc, 0);
    matrix * M = malloc_matrix();
    will_return_count(_mock_test_malloc, 0, 3);
    int ret = malloc_matrix_data(M, type, 3, 3, 20);
    for(int i = 0; i < M->x.sparse.na; i++) {
        M->x.sparse.a[i] = i * 10;
    }
    for(int i = 0; i < M->x.sparse.nia; i++) {
        M->x.sparse.ia[i] = i % 3;
    }
    for(int i = 0; i < M->x.sparse.nja; i++) {
        M->x.sparse.ja[i] = (i + 2) % 3;
    }
    if(type == CSC_SYMMETRIC) {
        for(int i = 0; i < 4; i++) {
            M->x.sparse.ja[i] = i;
        }
        M->x.sparse.ja[3] = 20;
    }
    else if(type == CSR_SYMMETRIC) {
        for(int i = 0; i < 4; i++) {
            M->x.sparse.ia[i] = i;
        }
        M->x.sparse.ia[3] = 20;
    }
    printf_matrix(M);
    assert_int_equal(ret, 1);
    M->type = type;
    printf_matrix(M);
    M = free_matrix(M);
}

void test_malloc_matrix_data_coo_symmetric(void ** state)
{
    core_malloc_matrix_data_symmetric(COO_SYMMETRIC);
}

void test_malloc_matrix_data_csc_symmetric(void ** state)
{
    core_malloc_matrix_data_symmetric(CSC_SYMMETRIC);
}

void test_malloc_matrix_data_csr_symmetric(void ** state)
{
    core_malloc_matrix_data_symmetric(CSR_SYMMETRIC);
}

static matrix * build_coo()
{
    matrix * M = malloc_matrix();
    assert_non_null(M); /* success */
    assert_int_equal(malloc_matrix_data(M, COO, 2, 2, 2), 1);
    M->x.sparse.a[0] = 1.0;
    M->x.sparse.a[1] = 2.0;
    M->x.sparse.ia[0] = 0;
    M->x.sparse.ia[1] = 1;
    M->x.sparse.ja[0] = 1;
    M->x.sparse.ja[1] = 0;
    printf_matrix(M);
    return M;
}

void test_malloc_matrix_data_coo_to_csc(void ** state)
{
    for(int i = 0; i < 4; i++) { /* fail 1st, 2nd, 3rd realloc; then success */
        will_return_count(_mock_test_malloc, 0, 4);
        matrix * M = build_coo();
        /* now test coo_to_csc */
        if(i > 0) {
            will_return_count(_mock_test_realloc, 0, i);
        }
        if(i < 3) {
            will_return(_mock_test_realloc, 1);
        }
        int ret = coo_to_csc(M);
        assert_int_equal(ret, (i < 3) ? 0 : 1);
        printf_matrix(M);
        free_matrix(M);
    }
}

void test_malloc_matrix_data_coo_to_csr(void ** state)
{
    for(int i = 0; i < 4; i++) { /* fail 1st, 2nd, 3rd realloc; then success */
        will_return_count(_mock_test_malloc, 0, 4);
        matrix * M = build_coo();
        /* now test coo_to_csr */
        if(i > 0) {
            will_return_count(_mock_test_realloc, 0, i);
        }
        if(i < 3) {
            will_return(_mock_test_realloc, 1);
        }
        int ret = coo_to_csr(M);
        assert_int_equal(ret, (i < 3) ? 0 : 1);
        printf_matrix(M);
        free_matrix(M);
    }
}

int main(void)
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_malloc_matrix),
        cmocka_unit_test(test_malloc_matrix_name_null),
        cmocka_unit_test(test_malloc_matrix_name_happy),
        cmocka_unit_test(test_malloc_matrix_name_sad),
        cmocka_unit_test(test_malloc_matrix_data_null),
        cmocka_unit_test(test_malloc_matrix_data_identity),
        cmocka_unit_test(test_malloc_matrix_data_coo),
        cmocka_unit_test(test_malloc_matrix_data_coo_symmetric),
        cmocka_unit_test(test_malloc_matrix_data_csc),
        cmocka_unit_test(test_malloc_matrix_data_csc_symmetric),
        cmocka_unit_test(test_malloc_matrix_data_coo_to_csc),
        cmocka_unit_test(test_malloc_matrix_data_csr),
        cmocka_unit_test(test_malloc_matrix_data_csr_symmetric),
        cmocka_unit_test(test_malloc_matrix_data_coo_to_csr),
        cmocka_unit_test(test_malloc_matrix_data_dense),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
