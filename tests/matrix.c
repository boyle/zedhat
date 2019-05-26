/* Copyright 2019, Alistair Boyle, 3-clause BSD License */
#include <stdarg.h> /* cmocka.h */
#include <stddef.h> /* cmocka.h */
#include <setjmp.h> /* cmocka.h */
#include <stdint.h> /* cmocka.h */
#include "cmocka.h"

#include <string.h> /* strlen + strcpy = mock strdup*/
#include <stdio.h> /* printf */

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
    if(mock()) {
        return NULL;
    }
    else {
        size_t size = strlen(str) + 1;
        char * n = __real__test_malloc(size, __FILE__, __LINE__);
        if(n != NULL) {
            strcpy(n, str);
        }
        return n;
    }
}

typedef struct {
    matrix * M;
} mystate;

static int group_setup( void ** state)
{
    *state = calloc(1, sizeof(mystate));
    return (*state == NULL);
}
static int group_teardown(void ** state)
{
   if(state == NULL)
      return 1;
    free(*state);
    *state = NULL;
    return 0;
}

static int test_teardown(void ** state)
{
   if(state == NULL) return 1;
    mystate * s = *state;
    free_matrix(s->M);
    s->M = NULL;
    return 0;
}

static int test_setup(void ** state)
{
   if(state == NULL) return 1;
    mystate * s = *state;
    will_return(__wrap__test_malloc, 0);
    s->M = malloc_matrix();
    return (s->M == NULL);
}

void test_malloc_matrix(void ** state)
{
    matrix * M = NULL;
    will_return(__wrap__test_malloc, 0);
    M = malloc_matrix();
    assert_non_null(M); /* success */
    free_matrix(M);
    will_return(__wrap__test_malloc, 1);
    M = malloc_matrix();
    assert_null(M); /* failure */
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

void test_malloc_matrix_name(void ** state)
{
    matrix * M = NULL;
    will_return(__wrap__test_malloc, 0);
    M = malloc_matrix();
    assert_non_null(M); /* success */
    int ret;
    ret = malloc_matrix_name(NULL, "A", "B", "C");
    assert_int_equal(ret, 0);
    int i;
    for(i = 0; i < 4; i++) { /* success */
        printf("i=%d ", i);
        if(i > 0) {
            will_return_count(__wrap_strdup, 0, i);
        }
        ret = malloc_matrix_name(M, i > 0 ? "A" : NULL, i > 1 ? "B" : NULL, i > 2 ? "C" : NULL);
        printf_matrix(M);
        assert_int_equal(ret, 1);
        if(i > 1) {
            will_return_count(__wrap_strdup, 0, i - 1);
        }
        if(i > 0) {
            will_return(__wrap_strdup, 1);
        }
        ret = malloc_matrix_name(M, i > 0 ? "A" : NULL, i > 1 ? "B" : NULL, i > 2 ? "C" : NULL);
        assert_int_equal(ret, i < 1);
    }
    assert_int_equal(ret, 0);
    free_matrix(M);
}


void test_malloc_matrix_data_null(void ** state)
{
    int ret = malloc_matrix_data(NULL, 0, 0, 0, 0);
    assert_int_equal(ret, 0);
}

void test_malloc_matrix_data_identity(void ** state)
{
    mystate * s = *state;
    assert_non_null(s->M);
    int ret = malloc_matrix_data(s->M, IDENTITY, 1, 1, 0);
    assert_int_equal(ret, 1);
}

void test_malloc_matrix_data_dense(void ** state)
{
    mystate * s = *state;
    will_return_count(__wrap__test_malloc, 0, 1);
    int ret = malloc_matrix_data(s->M, DENSE, 2, 2, 4);
    assert_int_equal(ret, 1);
    test_teardown(state);
}

void test_malloc_matrix_data_coo(void ** state)
{
    mystate * s = *state;
    will_return_count(__wrap__test_malloc, 0, 3);
    int ret = malloc_matrix_data(s->M, COO, 1, 1, 1);
    assert_int_equal(ret, 1);
    test_teardown(state);
}

void test_malloc_matrix_data_coo_symmetric(void ** state)
{
    mystate * s = *state;
    will_return_count(__wrap__test_malloc, 0, 3);
    int ret = malloc_matrix_data(s->M, COO_SYMMETRIC, 1, 1, 1);
    assert_int_equal(ret, 1);
    test_teardown(state);
}

int main(void)
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_malloc_matrix),
        cmocka_unit_test(test_malloc_matrix_name),
        cmocka_unit_test_setup_teardown(test_malloc_matrix_data_null, test_setup, test_teardown),
        cmocka_unit_test_setup_teardown(test_malloc_matrix_data_identity, test_setup, test_teardown),
        cmocka_unit_test_setup_teardown(test_malloc_matrix_data_dense, test_setup, test_teardown),
        cmocka_unit_test_setup_teardown(test_malloc_matrix_data_coo, test_setup, test_teardown),
        cmocka_unit_test_setup_teardown(test_malloc_matrix_data_coo_symmetric, test_setup, test_teardown),
    };
    return cmocka_run_group_tests(tests, group_setup, group_teardown);
}
