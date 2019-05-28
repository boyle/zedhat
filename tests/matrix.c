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

void * __real__test_malloc(const size_t size, const char * file, const int line);
void * __wrap__test_malloc(size_t size)
{
    if(mock()) {
        return NULL;
    }
    else {
        return __real__test_malloc(size, __FILE__, __LINE__);
    }
}

char * __wrap_strdup(const char * str)
{
    size_t size = strlen(str) + 1;
//        char * n = __real__test_malloc(size, __FILE__, __LINE__);
    char * n = malloc(size);
    if(n != NULL) {
        strcpy(n, str);
    }
    return n;
}

void test_malloc_matrix(void ** state)
{
    matrix * M = NULL;
    will_return(__wrap__test_malloc, 0);
    M = malloc_matrix();
    assert_non_null(M); /* success */
    M = free_matrix(M);
    will_return(__wrap__test_malloc, 1);
    M = malloc_matrix();
    assert_null(M); /* failure */
    M = free_matrix(M);
}

void printf_matrix(const matrix * M)
{
    printf("%s [%3s]: %zux%zu %5s -- %s\n",
           M->symbol ? M->symbol : "?",
           M->units ? M->units : "?",
           M->m, M->n,
           M->type == DENSE ? "DENSE" : M->type == COO ? "COO" : "?",
           M->name ? M->name : "?");
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
    will_return(__wrap__test_malloc, 0);
    matrix * M = malloc_matrix();
    assert_non_null(M);
    int i;
    for(i = 0; i < 8; i++) {
        if(bitcnt(i) > 0) {
            will_return_count(__wrap__test_malloc, 0, bitcnt(i));
        }
        int ret = malloc_matrix_name(M, i & 1 ? "A" : NULL, i & 2 ? "B" : NULL, i & 4 ? "C" : NULL);
        printf("%d='b%d%d%d ", i, (i >> 0) & 1, (i >> 1) & 1, (i >> 2) & 1);
        printf_matrix(M);
        assert_int_equal(ret, 1);
    }
    M = free_matrix(M);
}

void test_malloc_matrix_name_sad(void ** state)
{
    will_return(__wrap__test_malloc, 0);
    matrix * M = malloc_matrix();
    assert_non_null(M);
    int i;
    for(i = 1; i < 8; i++) {
        if(bitcnt(i) > 1) {
            will_return_count(__wrap__test_malloc, 0, bitcnt(i) - 1);
        }
        will_return(__wrap__test_malloc, 1);
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
    will_return(__wrap__test_malloc, 0);
    matrix * M = malloc_matrix();
    assert_non_null(M);
    int ret = malloc_matrix_data(M, IDENTITY, 1, 1, 0);
    assert_int_equal(ret, 1);
    M = free_matrix(M);
}

void test_malloc_matrix_data_dense(void ** state)
{
    int ret;
    will_return(__wrap__test_malloc, 0);
    matrix * M = malloc_matrix();
    will_return(__wrap__test_malloc, 0);
    ret = malloc_matrix_data(M, DENSE, 2, 2, 4);
    assert_int_equal(ret, 1);
    M = free_matrix(M);
    will_return(__wrap__test_malloc, 0);
    M = malloc_matrix();
    will_return(__wrap__test_malloc, 1);
    ret = malloc_matrix_data(M, DENSE, 2, 2, 4);
    assert_int_equal(ret, 0);
    M = free_matrix(M);
}

void test_malloc_matrix_data_coo(void ** state)
{
    int i;
    for(i = 0; i < 8; i++) {
        will_return(__wrap__test_malloc, 0);
        matrix * M = malloc_matrix();
        will_return(__wrap__test_malloc, (i >> 0) & 1);
        will_return(__wrap__test_malloc, (i >> 1) & 1);
        will_return(__wrap__test_malloc, (i >> 2) & 1);
        int ret = malloc_matrix_data(M, COO, 1, 1, 1);
        assert_int_equal(ret, bitcnt(i) > 0 ? 0 : 1);
        M = free_matrix(M);
    }
}

void test_malloc_matrix_data_coo_symmetric(void ** state)
{
    will_return(__wrap__test_malloc, 0);
    matrix * M = malloc_matrix();
    will_return_count(__wrap__test_malloc, 0, 3);
    int ret = malloc_matrix_data(M, COO_SYMMETRIC, 1, 1, 1);
    assert_int_equal(ret, 1);
    M->type = COO_SYMMETRIC;
    M = free_matrix(M);
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
        cmocka_unit_test(test_malloc_matrix_data_dense),
        cmocka_unit_test(test_malloc_matrix_data_coo),
        cmocka_unit_test(test_malloc_matrix_data_coo_symmetric),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
