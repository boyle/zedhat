/* Copyright 2018-2019, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdarg.h> /* cmocka.h */
#include <stddef.h> /* cmocka.h */
#include <setjmp.h> /* cmocka.h */
#include <stdint.h> /* cmocka.h */

#include <stdio.h> /* cmocka.h */
#include <string.h> /* cmocka.h, memcpy */
#include <float.h> /* DBL_EPSILON */

#include "cmocka.h"
#include "inv.h"

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
    double cond[2] = {0.5, 0.5};
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
    double measgain[2] = {2.5, 1};
    double zc[2] = {1e-2, 1e-3};
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
    m.measgain = &(measgain[0]);
    m.n_stimmeas = 2;
    m.params = &(cond[0]);
    m.n_params[0] = 2;
    m.n_params[1] = 1;
    m.n_zc = 2;
    m.zc = &(zc[0]);
    m.hp = 1e-5;
    double meas[4] = {
        5.002749999984084, -2.001099999993633, /* b_1 = v_h (@ 0.5 S/m)*/
        4.548204545438010, -1.819281818175204, /* b_2 = v_i (@ 0.55 S/m) */
    };
    m.n_data[0] = 2;
    m.n_data[1] = 2;
    m.data = &(meas[0]);
    double x[2] = {0};
    printf_model(&m, 2);
    will_return_always(_mock_test_malloc, 0);
    will_return_always(_mock_test_realloc, 0);
    assert_int_equal(inv_solve(&m, &x[0]), 1);
    printf("x = %g %g\n", x[0], x[1]);
    const double expect = 0.045454545454872;
    for(int i = 0; i < 2; i++) {
        printf("Δ = %g - %g = %g\n", x[i], expect, x[i] - expect);
    }
    const double reltol = 2e3 * DBL_EPSILON;
    printf("relative tolerance = %g\n", reltol);
    test_free(m.elec_to_sys);
    assert_double_equal(x[0], expect, reltol);
    assert_double_equal(x[1], expect, reltol);
}

// #define NUM_MALLOCS_PMAP (5 + 3 + 3 + 4 + 1)
// #define NUM_REALLOCS_PMAP (3 + 3)
// void test_2d_resistor_pmap (void ** state)
// {
//     /* 2D square:
//      *    first four nodes at z=0
//      *    first two elems */
//     double nodes[4][2] = {
//         {0, 0},
//         {1, 0},
//         {0, 1},
//         {1, 1},
//     };
//     int elems[2][3] = { /* from netgen cube.geo */
//         {1, 2, 4},
//         {1, 3, 4},
//     };
//     double cond[1] = {0.5};
//     int bc[4] = {0, 1, 0, 2};
//     int se[4][2] = {
//         {1, 2},
//         {2, 4},
//         {4, 3},
//         {3, 1},
//     };
//     int stimmeas[2][4] = {
//         {1, 2, 1, 2},
//         {2, 1, 1, 2},
//     };
//     double measgain[2] = {1, 1};
//     int pmap_elem[2] = {1, 2};
//     int pmap_param[2] = {1, 1};
//     double pmap_frac[2] = {1.0, 1.0};
//     model m = {{ 0 }};
//     m.fwd.dim = 2;
//     m.fwd.elems = &(elems[0][0]);
//     m.fwd.nodes = &(nodes[0][0]);
//     m.fwd.n_elems = 2;
//     m.fwd.n_nodes = 4;
//     m.fwd.bc = &(bc[0]);
//     m.fwd.surfaceelems = &(se[0][0]);
//     m.fwd.n_se = 4;
//     m.fwd.n_pmap = 2;
//     m.fwd.pmap_elem = &(pmap_elem[0]);
//     m.fwd.pmap_param = &(pmap_param[0]);
//     m.fwd.pmap_frac = &(pmap_frac[0]);
//     m.stimmeas = &(stimmeas[0][0]);
//     m.measgain = &(measgain[0]);
//     m.n_stimmeas = 2;
//     m.params = &(cond[0]);
//     m.n_params[0] = 1;
//     m.n_params[1] = 1;
//     double meas[2] = {0};
//     will_return_always(_mock_test_malloc, 0);
//     will_return_always(_mock_test_realloc, 0);
//     assert_int_equal(fwd_solve(&m, &meas[0]), 1);
//     printf("meas = %g %g\n", meas[0], meas[1]);
//     const double default_zc = 1e-2;
//     const double expect = 2.0 + 2.0 * default_zc;
//     for(int i = 0; i < 2; i++) {
//         const double expecti = i == 0 ? +expect : -expect;
//         printf("Δ = %g - %g = %g\n", meas[i], expecti, meas[i] - expecti);
//     }
//     const double reltol = 4e3 * DBL_EPSILON;
//     printf("relative tolerance = %g\n", reltol);
//     assert_double_equal(meas[0], +expect, reltol);
//     assert_double_equal(meas[1], -expect, reltol);
//     double J[2][1] = {{0}};
//     assert_int_equal(calc_jacobian(&m, (&J[0][0])), 1);
//     printf("Jacobian:\n");
//     for(int i = 0; i < m.n_stimmeas; i++) {
//         for(int j = 0; j < 1; j++) {
//             printf(" %18g", J[i][j]);
//         }
//         printf("\n");
//     }
//     assert_double_equal(J[0][0], -2.0 * 2, reltol);
//     assert_double_equal(J[1][0], +2.0 * 2, reltol);
//     test_free(m.elec_to_sys);
//     test_free(m.zc);
// }
//
// void test_3d_resistor_cem (void ** state)
// {
//     /* 3D cube:
//      *    all nodes used as in elems */
//     double nodes[8][3] = {
//         {0, 0, 0},
//         {1, 0, 0},
//         {0, 1, 0},
//         {1, 1, 0},
//         {0, 0, 1},
//         {1, 0, 1},
//         {0, 1, 1},
//         {1, 1, 1},
//     };
//     int elems[6][4] = { /* from netgen cube.geo */
//         {2, 5, 1, 4},
//         {1, 5, 3, 4},
//         {3, 5, 7, 8},
//         {8, 4, 5, 6},
//         {2, 6, 5, 4},
//         {8, 4, 3, 5},
//     };
//     double cond[6] = {1, 1, 1, 1, 1, 1};
//     int bc[12] = {2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
//     int se[12][3] = {
//         {1, 2, 5},
//         {2, 5, 6},
//         {4, 3, 8},
//         {3, 8, 7},
//         {6, 8, 2},
//         {8, 2, 4},
//         {1, 5, 3},
//         {5, 3, 7},
//         {5, 7, 8},
//         {5, 8, 6},
//         {1, 3, 2},
//         {3, 2, 4},
//     };
//     int stimmeas[2][4] = {
//         {1, 2, 1, 2},
//         {2, 1, 1, 2},
//     };
//     double measgain[2] = {1, 1};
//     model m = {{ 0 }};
//     m.fwd.dim = 3;
//     m.fwd.elems = &(elems[0][0]);
//     m.fwd.nodes = &(nodes[0][0]);
//     m.fwd.n_elems = 6;
//     m.fwd.n_nodes = 8;
//     m.fwd.bc = &(bc[0]);
//     m.fwd.surfaceelems = &(se[0][0]);
//     m.fwd.n_se = 12;
//     m.stimmeas = &(stimmeas[0][0]);
//     m.measgain = &(measgain[0]);
//     m.n_stimmeas = 2;
//     m.params = &(cond[0]);
//     m.n_params[0] = 6;
//     m.n_params[1] = 1;
//     double meas[2] = {0};
//     will_return_always(_mock_test_malloc, 0);
//     will_return_always(_mock_test_realloc, 0);
//     assert_int_equal(fwd_solve(&m, &meas[0]), 1);
//     printf("meas = %g %g\n", meas[0], meas[1]);
//     const double default_zc = 1e-2;
//     const double expect = 1.0 + 2.0 * default_zc;
//     for(int i = 0; i < 2; i++) {
//         const double expecti = i == 0 ? +expect : -expect;
//         printf("Δ = %g - %g = %g\n", meas[i], expecti, meas[i] - expecti);
//     }
//     const double reltol = 1e3 * DBL_EPSILON;
//     printf("relative tolerance = %g\n", reltol);
//     assert_double_equal(meas[0], +expect, reltol);
//     assert_double_equal(meas[1], -expect, reltol);
//     double J[2][6] = {{0}};
//     assert_int_equal(calc_jacobian(&m, (&J[0][0])), 1);
//     printf("Jacobian:\n");
//     for(int i = 0; i < m.n_stimmeas; i++) {
//         for(int j = 0; j < m.fwd.n_elems; j++) {
//             printf(" %18g", J[i][j]);
//         }
//         printf("\n");
//     }
//     for(int j = 0; j < m.fwd.n_elems; j++) {
//         assert_double_equal(J[0][j], -1.0 / 6.0, reltol);
//     }
//     for(int j = 0; j < m.fwd.n_elems; j++) {
//         assert_double_equal(J[1][1], +1.0 / 6.0, reltol);
//     }
//     test_free(m.elec_to_sys);
//     test_free(m.zc);
// }
//
// /* TODO currently no PEM test cases */

void test_inv_solve_fail (void ** state)
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
    double measgain[2] = {1, 1};
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
    m.measgain = &(measgain[0]);
    m.n_stimmeas = 2;
    m.params = &(cond[0]);
    m.n_params[0] = 2;
    m.n_params[1] = 1;
    m.hp = 1e-5;
    double meas[4] = {0};
    m.n_data[0] = 2;
    m.n_data[1] = 2;
    m.data = &(meas[0]);
    double x[2] = {0};
    printf("\n"); printf_model(&m, 2); printf("\n");
    will_return(_mock_test_malloc, 1);
    assert_int_equal(inv_solve(&m, &x[0]), 0);
    test_free(m.elec_to_sys);
    test_free(m.zc);
    m.zc = NULL;
    m.elec_to_sys = NULL;
}

int main(void)
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_2d_resistor_cem),
//        cmocka_unit_test(test_3d_resistor_cem),
//        cmocka_unit_test(test_2d_resistor_pmap),
        cmocka_unit_test(test_inv_solve_fail),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
