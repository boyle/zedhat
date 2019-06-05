/* Copyright 2018-2019, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdarg.h> /* cmocka.h */
#include <stddef.h> /* cmocka.h */
#include <setjmp.h> /* cmocka.h */
#include <stdint.h> /* cmocka.h */

#include <stdio.h> /* cmocka.h */
#include <string.h> /* cmocka.h, memcpy */
#include <float.h> /* DBL_EPSILON */
#include <math.h> /* fabs, pow */

#include <cholmod.h>

#include "cmocka.h"
#include "fwd.h"

void * _mock_test_malloc(size_t size, const char * file, const int line)
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

void test_2d_resistor_cem (void ** state)
{
    /* 2D square:
     *    first four nodes at z=0
     *    first two elems */
    double nodes[4][2] = {
        {0, 0},
        {1, 0},
        {0, 1},
        {1, 1},
    };
    int elems[2][3] = { /* from netgen cube.geo */
        {1, 2, 4},
        {1, 3, 4},
    };
    double cond[2] = {1, 1};
    int bc[4] = {0, 1, 0, 2};
    int se[4][2] = {
        {1, 2},
        {2, 4},
        {4, 3},
        {3, 1},
    };
    int stimmeas[2][4] = {
        {1, 2, 1, 2},
        {2, 1, 1, 2},
    };
    model m = {{ 0 }};
    m.fwd.dim = 2;
    m.fwd.elems = &(elems[0][0]);
    m.fwd.nodes = &(nodes[0][0]);
    m.fwd.n_elems = 2;
    m.fwd.n_nodes = 4;
    m.fwd.bc = &(bc[0]);
    m.fwd.surfaceelems = &(se[0][0]);
    m.fwd.n_se = 4;
    m.stimmeas = &(stimmeas[0][0]);
    m.n_stimmeas = 2;
    m.params = &(cond[0]);
    m.n_params[0] = 2;
    m.n_params[1] = 1;
    double meas[2] = {0};
    will_return_count(_mock_test_malloc, 0, 5 + 3);
    will_return_count(_mock_test_realloc, 0, 3);
    assert_int_equal(fwd_solve(&m, &meas[0]), 1);
    test_free(m.elec_to_sys);
    printf("meas = %g %g\n", meas[0], meas[1]);
    const double expect = 1.0 + 2.0 * 1e-3;
    for(int i = 0; i < 2; i++) {
        const double expecti = i == 0 ? +expect : -expect;
        printf("Δ = %g - %g = %g\n", meas[i] / 10.0, expecti, meas[i] / 10.0 - expecti);
    }
    const double reltol = 1e3 * DBL_EPSILON;
    printf("relative tolerance = %g\n", reltol);
    assert_double_equal(meas[0] / 10.0, +(1.0 + 2.0 * 1e-3), reltol); /* TODO rm /10 */
    assert_double_equal(meas[1] / 10.0, -(1.0 + 2.0 * 1e-3), reltol); /* TODO rm /10 */
}

void test_3d_resistor_cem (void ** state)
{
    /* 3D cube:
     *    all nodes used as in elems */
    double nodes[8][3] = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 1},
    };
    int elems[6][4] = { /* from netgen cube.geo */
        {2, 5, 1, 4},
        {1, 5, 3, 4},
        {3, 5, 7, 8},
        {8, 4, 5, 6},
        {2, 6, 5, 4},
        {8, 4, 3, 5},
    };
    double cond[6] = {1, 1, 1, 1, 1, 1};
    int bc[12] = {2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    int se[12][3] = {
        {1, 2, 5},
        {2, 5, 6},
        {4, 3, 8},
        {3, 8, 7},
        {6, 8, 2},
        {8, 2, 4},
        {1, 5, 3},
        {5, 3, 7},
        {5, 7, 8},
        {5, 8, 6},
        {1, 3, 2},
        {3, 2, 4},
    };
    int stimmeas[2][4] = {
        {1, 2, 1, 2},
        {2, 1, 1, 2},
    };
    model m = {{ 0 }};
    m.fwd.dim = 3;
    m.fwd.elems = &(elems[0][0]);
    m.fwd.nodes = &(nodes[0][0]);
    m.fwd.n_elems = 6;
    m.fwd.n_nodes = 8;
    m.fwd.bc = &(bc[0]);
    m.fwd.surfaceelems = &(se[0][0]);
    m.fwd.n_se = 12;
    m.stimmeas = &(stimmeas[0][0]);
    m.n_stimmeas = 2;
    m.params = &(cond[0]);
    m.n_params[0] = 6;
    m.n_params[1] = 1;
    double meas[2] = {0};
    will_return_count(_mock_test_malloc, 0, 5 + 3);
    will_return_count(_mock_test_realloc, 0, 3);
    assert_int_equal(fwd_solve(&m, &meas[0]), 1);
    test_free(m.elec_to_sys);
    printf("meas = %g %g\n", meas[0], meas[1]);
    const double expect = 1.0 + 2.0 * 1e-3;
    for(int i = 0; i < 2; i++) {
        const double expecti = i == 0 ? +expect : -expect;
        printf("Δ = %g - %g = %g\n", meas[i] / 10, expecti, meas[i] / 10 - expecti);
    }
    const double reltol = 1e3 * DBL_EPSILON;
    printf("relative tolerance = %g\n", reltol);
    assert_double_equal(meas[0] / 10, +expect, reltol); /* TODO rm /10 */
    assert_double_equal(meas[1] / 10, -expect, reltol); /* TODO rm /10 */
}

/* TODO currently no PEM test cases */

void test_fail (void ** state)
{
    /* 2D square:
     *    first four nodes at z=0
     *    first two elems */
    double nodes[4][2] = {
        {0, 0},
        {1, 0},
        {0, 1},
        {1, 1},
    };
    int elems[2][3] = { /* from netgen cube.geo */
        {1, 2, 4},
        {1, 3, 4},
    };
    double cond[2] = {1, 1};
    int bc[4] = {0, 1, 0, 2};
    int se[4][2] = {
        {1, 2},
        {2, 4},
        {4, 3},
        {3, 0}, /* PEM for line coverage; no check on measurements */
    };
    int stimmeas[2][4] = {
        {1, 2, 1, 2},
        {2, 1, 1, 2},
    };
    model m = {{ 0 }};
    m.fwd.dim = 2;
    m.fwd.elems = &(elems[0][0]);
    m.fwd.nodes = &(nodes[0][0]);
    m.fwd.n_elems = 2;
    m.fwd.n_nodes = 4;
    m.fwd.bc = &(bc[0]);
    m.fwd.surfaceelems = &(se[0][0]);
    m.fwd.n_se = 4;
    m.stimmeas = &(stimmeas[0][0]);
    m.n_stimmeas = 2;
    m.params = &(cond[0]);
    m.n_params[0] = 2;
    m.n_params[1] = 1;
    double meas[2] = {0};
    /* check_model fail */
    se[3][0] = 0;
    assert_int_equal(fwd_solve(&m, &meas[0]), 0);
    se[3][0] = 3;
    /* various malloc failures */
    will_return(_mock_test_malloc, 1);
    assert_int_equal(fwd_solve(&m, &meas[0]), 0);
    for(int i = 1; i < 8; i++) {
        will_return(_mock_test_malloc, 0);
        will_return(_mock_test_malloc, (i >> 0) & 1);
        will_return(_mock_test_malloc, (i >> 1) & 1);
        will_return(_mock_test_malloc, (i >> 2) & 1);
        assert_int_equal(fwd_solve(&m, &meas[0]), 0);
    }
    will_return_count(_mock_test_malloc, 0, 4);
    will_return(_mock_test_malloc, 1);
    assert_int_equal(fwd_solve(&m, &meas[0]), 0);
    will_return_count(_mock_test_malloc, 0, 4 + 3);
    will_return(_mock_test_malloc, 1);
    assert_int_equal(fwd_solve(&m, &meas[0]), 0);
    will_return_count(_mock_test_malloc, 0, 4 + 3);
    will_return_count(_mock_test_realloc, 0, 2);
    will_return(_mock_test_realloc, 1);
    assert_int_equal(fwd_solve(&m, &meas[0]), 0);
    /* build system matrix failure: singular element volume */
    elems[0][1] = 1;
    will_return_count(_mock_test_malloc, 0, 4 + 3);
    //will_return_count(_mock_test_realloc, 0, 3);
    assert_int_equal(fwd_solve(&m, &meas[0]), 0);
    elems[0][1] = 2;
    /* passing example again */
    will_return_count(_mock_test_malloc, 0, 4 + 3);
    will_return_count(_mock_test_realloc, 0, 3);
    assert_int_equal(fwd_solve(&m, &meas[0]), 1);
    test_free(m.elec_to_sys);
}


int main(void)
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_2d_resistor_cem),
        cmocka_unit_test(test_3d_resistor_cem),
        cmocka_unit_test(test_fail),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
