/* Copyright 2018-2019, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdarg.h> /* cmocka.h */
#include <stddef.h> /* cmocka.h */
#include <setjmp.h> /* cmocka.h */
#include <stdint.h> /* cmocka.h */

#include <stdio.h> /* cmocka.h */
#include <string.h> /* cmocka.h, memcpy */
#include <float.h> /* DBL_EPSILON */

#include <cholmod.h>

#include "cmocka.h"
#include "fwd.h"

/* determined from "happy" runs */
#define NUM_FWD_MALLOCS (9)
#define NUM_FWD_REALLOCS (3)
#define NUM_FWD_MALLOCS_PMAP (16)
#define NUM_FWD_REALLOCS_PMAP (6)

#define NUM_JAC_MALLOCS (17)
#define NUM_JAC_REALLOCS (3)
#define NUM_JAC_MALLOCS_PMAP (24)
#define NUM_JAC_REALLOCS_PMAP (6)

static int malloc_count = 0;
static int malloc_fail_at = 0;

static int realloc_count = 0;
static int realloc_fail_at = 0;

void reset_mock_count()
{
    malloc_count = 0;
    realloc_count = 0;
    malloc_fail_at = 0;
    realloc_fail_at = 0;
}

void * _mock_test_malloc(size_t size, const char * file, const int line)
{
    malloc_count++;
    if(malloc_count == malloc_fail_at) {
        printf("  malloc x (#%d) @ %s:%d\n", malloc_count, file, line);
        return NULL;
    }
    else {
        printf("  malloc ✓ (#%d) @ %s:%d\n", malloc_count, file, line);
        return _test_malloc(size, file, line);
    }
}

void * _mock_test_realloc(void * const ptr, const size_t size, const char * file, const int line)
{
    realloc_count++;
    if(realloc_count == realloc_fail_at) {
        printf("  realloc x (#%d) @ %s:%d\n", realloc_count, file, line);
        return NULL;
    }
    else {
        printf("  realloc ✓ (#%d) @ %s:%d\n", realloc_count, file, line);
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

void test_malloc_and_free (void ** state)
{
    reset_mock_count();
    void *p = _mock_test_malloc(100, __FILE__, __LINE__);
    void *q = _mock_test_malloc(200, __FILE__, __LINE__);
    p = _mock_test_realloc(p, 150, __FILE__, __LINE__);
    test_free(p);
    test_free(q);
    assert_int_equal(malloc_count, 2);
    assert_int_equal(realloc_count, 1);
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
    double measgain[2] = {2.5, 1};
    double zc[2] = {1e-4, 1e-3};
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
    double meas[2] = {0};
    m.n_data[0] = 2;
    m.n_data[1] = 1;
    m.data = &(meas[0]);
    reset_mock_count();
    assert_int_equal(fwd_solve(&m, &meas[0]), 1); /* prove an initial success */
    printf("fwd_solve()  mallocs = %d, reallocs = %d\n", malloc_count, realloc_count);
    printf("meas = %g %g\n", meas[0], meas[1]);
    const double default_zc = zc[0] + zc[1];
    const double expect = 1.0 + default_zc;
    for(int i = 0; i < 2; i++) {
        const double expecti = (i == 0 ? +expect * 2.5 : -expect);
        printf("Δ = %g - %g = %g\n", meas[i], expecti, meas[i] - expecti);
    }
    const double reltol = 2e4 * DBL_EPSILON;
    printf("relative tolerance = %g\n", reltol);
    assert_double_equal(meas[0], +expect * 2.5, reltol);
    assert_double_equal(meas[1], -expect, reltol);
    printf_model(&m, 2);
    reset_mock_count();
    double J[4] = {0};
    assert_int_equal(calc_jacobian(&m, (&J[0])), 1);
    printf("calc_jacobian()  mallocs = %d, reallocs = %d\n", malloc_count, realloc_count);
    printf("Jacobian:\n");
    for(int i = 0; i < m.n_stimmeas; i++) {
        for(int j = 0; j < m.fwd.n_elems; j++) {
            printf(" %18g", J[i + j * m.n_stimmeas]);
        }
        printf("\n");
    }
    assert_double_equal(J[0], -1.25, reltol);
    assert_double_equal(J[2], -1.25, reltol);
    assert_double_equal(J[1], +0.50, reltol);
    assert_double_equal(J[3], +0.50, reltol);
    test_free(m.elec_to_sys);
}

void test_2d_resistor_pmap (void ** state)
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
    double cond[1] = {0.5};
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
    double measgain[2] = {1, 1};
    int pmap_elem[2] = {1, 2};
    int pmap_param[2] = {1, 1};
    double pmap_frac[2] = {1.0, 1.0};
    model m = {{ 0 }};
    m.fwd.dim = 2;
    m.fwd.elems = &(elems[0][0]);
    m.fwd.nodes = &(nodes[0][0]);
    m.fwd.n_elems = 2;
    m.fwd.n_nodes = 4;
    m.fwd.bc = &(bc[0]);
    m.fwd.surfaceelems = &(se[0][0]);
    m.fwd.n_se = 4;
    m.fwd.n_pmap = 2;
    m.fwd.pmap_elem = &(pmap_elem[0]);
    m.fwd.pmap_param = &(pmap_param[0]);
    m.fwd.pmap_frac = &(pmap_frac[0]);
    m.stimmeas = &(stimmeas[0][0]);
    m.measgain = &(measgain[0]);
    m.n_stimmeas = 2;
    m.params = &(cond[0]);
    m.n_params[0] = 1;
    m.n_params[1] = 1;
    double meas[2] = {0};
    reset_mock_count();
    assert_int_equal(fwd_solve(&m, &meas[0]), 1); /* prove an initial success */
    printf("fwd_solve() PMAP mallocs = %d, reallocs = %d\n", malloc_count, realloc_count);
    printf("meas = %g %g\n", meas[0], meas[1]);
    const double default_zc = 1e-2;
    const double expect = 2.0 + 2.0 * default_zc;
    for(int i = 0; i < 2; i++) {
        const double expecti = i == 0 ? +expect : -expect;
        printf("Δ = %g - %g = %g\n", meas[i], expecti, meas[i] - expecti);
    }
    const double reltol = 4e3 * DBL_EPSILON;
    printf("relative tolerance = %g\n", reltol);
    assert_double_equal(meas[0], +expect, reltol);
    assert_double_equal(meas[1], -expect, reltol);
    reset_mock_count();
    double J[2] = {0};
    assert_int_equal(calc_jacobian(&m, (&J[0])), 1);
    printf("calc_jacobian() PMAP mallocs = %d, reallocs = %d\n", malloc_count, realloc_count);
    printf("Jacobian:\n");
    for(int i = 0; i < m.n_stimmeas; i++) {
        for(int j = 0; j < 1; j++) {
            printf(" %18g", J[i + j * m.n_stimmeas]);
        }
        printf("\n");
    }
    assert_double_equal(J[0], -2.0 * 2, reltol);
    assert_double_equal(J[1], +2.0 * 2, reltol);
    test_free(m.elec_to_sys);
    test_free(m.zc);
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
    double measgain[2] = {1, 1};
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
    m.measgain = &(measgain[0]);
    m.n_stimmeas = 2;
    m.params = &(cond[0]);
    m.n_params[0] = 6;
    m.n_params[1] = 1;
    double meas[2] = {0};
    assert_int_equal(fwd_solve(&m, &meas[0]), 1);
    printf("meas = %g %g\n", meas[0], meas[1]);
    const double default_zc = 1e-2;
    const double expect = 1.0 + 2.0 * default_zc;
    for(int i = 0; i < 2; i++) {
        const double expecti = i == 0 ? +expect : -expect;
        printf("Δ = %g - %g = %g\n", meas[i], expecti, meas[i] - expecti);
    }
    const double reltol = 1e3 * DBL_EPSILON;
    printf("relative tolerance = %g\n", reltol);
    assert_double_equal(meas[0], +expect, reltol);
    assert_double_equal(meas[1], -expect, reltol);
    double J[12] = {0};
    assert_int_equal(calc_jacobian(&m, (&J[0])), 1);
    printf("Jacobian:\n");
    for(int i = 0; i < m.n_stimmeas; i++) {
        for(int j = 0; j < m.fwd.n_elems; j++) {
            printf(" %18g", J[i + j * m.n_stimmeas]);
        }
        printf("\n");
    }
    for(int j = 0; j < m.fwd.n_elems; j++) {
        assert_double_equal(J[0 + j * m.n_stimmeas], -1.0 / 6.0, reltol);
    }
    for(int j = 0; j < m.fwd.n_elems; j++) {
        assert_double_equal(J[1 + j * m.n_stimmeas], +1.0 / 6.0, reltol);
    }
    test_free(m.elec_to_sys);
    test_free(m.zc);
}

/* TODO currently no PEM test cases */

void test_fwd_solve_fail (void ** state)
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
    double cond[1] = {1};
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
    int pmap_elem[2] = {1, 2};
    int pmap_param[2] = {1, 1};
    double pmap_frac[2] = {1.0, 1.0};
    model m = {{ 0 }};
    m.fwd.dim = 2;
    m.fwd.elems = &(elems[0][0]);
    m.fwd.nodes = &(nodes[0][0]);
    m.fwd.n_elems = 2;
    m.fwd.n_nodes = 4;
    m.fwd.bc = &(bc[0]);
    m.fwd.surfaceelems = &(se[0][0]);
    m.fwd.n_se = 4;
    m.fwd.n_pmap = 2;
    m.fwd.pmap_elem = &(pmap_elem[0]);
    m.fwd.pmap_param = &(pmap_param[0]);
    m.fwd.pmap_frac = &(pmap_frac[0]);
    m.stimmeas = &(stimmeas[0][0]);
    m.measgain = &(measgain[0]);
    m.n_stimmeas = 2;
    m.params = &(cond[0]);
    m.n_params[0] = 1;
    m.n_params[1] = 1;
    double meas[2] = {0};
    printf("\n"); printf_model(&m, 2); printf("\n");
    /* check_model fail: bad boundary condition */
    se[3][0] = 0;
    reset_mock_count();
    assert_int_equal(fwd_solve(&m, &meas[0]), 0);
    assert_int_equal(malloc_count, 0);
    assert_int_equal(realloc_count, 0);
    se[3][0] = 3;
    /* various malloc failures */
    printf("malloc sad\n");
    for(int i = 0, mallocs = 0, reallocs = 0; i < NUM_FWD_MALLOCS_PMAP + NUM_FWD_REALLOCS_PMAP; i++) {
        reset_mock_count();
        switch(i) {
        case 11:
        case 12:
        case 13:
        case 19:
        case 20:
        case 21:
            reallocs++;
            realloc_fail_at = reallocs;
            break;
        default:
            mallocs++;
            malloc_fail_at = mallocs;
        }
        printf("[i=%d] mallocs = %d/%d, reallocs = %d/%d\n", i,
               mallocs, NUM_FWD_MALLOCS_PMAP,
               reallocs, NUM_FWD_REALLOCS_PMAP);
        assert_int_equal(fwd_solve(&m, &meas[0]), 0);
        // expected malloc & realloc
        assert_int_equal(malloc_count, mallocs);
        assert_int_equal(realloc_count, reallocs);
    }
    /* build system matrix failure: singular element volume */
    printf("singular element (zero volume)\n");
    reset_mock_count();
    elems[0][1] = 1;
    assert_int_equal(fwd_solve(&m, &meas[0]), 0);
    assert_int_equal(malloc_count, NUM_FWD_MALLOCS_PMAP);
    assert_int_equal(realloc_count, NUM_FWD_REALLOCS_PMAP-3);
    elems[0][1] = 2;
    /* passing example again */
    printf("happy\n");
    reset_mock_count();
    assert_int_equal(fwd_solve(&m, &meas[0]), 1);
    assert_int_equal(malloc_count, NUM_FWD_MALLOCS_PMAP);
    assert_int_equal(realloc_count, NUM_FWD_REALLOCS_PMAP);
    test_free(m.elec_to_sys);
    test_free(m.zc);
}

void test_calc_jacobian_fail (void ** state)
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
    double cond[1] = {1};
    int bc[4] = {0, 1, 0, 2};
    int se[4][2] = {
        {1, 2},
        {2, 4},
        {4, 3},
        {3, 0},
    };
    int stimmeas[2][4] = {
        {1, 2, 1, 2},
        {2, 1, 1, 2},
    };
    double measgain[2] = {1, 1};
    int pmap_elem[2] = {1, 2};
    int pmap_param[2] = {1, 1};
    double pmap_frac[2] = {1.0, 1.0};
    model m = {{ 0 }};
    m.fwd.dim = 2;
    m.fwd.elems = &(elems[0][0]);
    m.fwd.nodes = &(nodes[0][0]);
    m.fwd.n_elems = 2;
    m.fwd.n_nodes = 4;
    m.fwd.bc = &(bc[0]);
    m.fwd.surfaceelems = &(se[0][0]);
    m.fwd.n_se = 4;
    m.fwd.n_pmap = 2;
    m.fwd.pmap_elem = &(pmap_elem[0]);
    m.fwd.pmap_param = &(pmap_param[0]);
    m.fwd.pmap_frac = &(pmap_frac[0]);
    m.stimmeas = &(stimmeas[0][0]);
    m.measgain = &(measgain[0]);
    m.n_stimmeas = 2;
    m.params = &(cond[0]);
    m.n_params[0] = 1;
    m.n_params[1] = 1;
    double J[2] = {0};
    printf("\n"); printf_model(&m, 2); printf("\n");
    printf("happy\n");
    reset_mock_count();
    assert_int_equal(calc_jacobian(&m, &J[0]), 1);
    assert_int_equal(malloc_count, NUM_JAC_MALLOCS_PMAP);
    assert_int_equal(realloc_count, NUM_JAC_REALLOCS_PMAP);
    /* check_model fail: bad boundary condition */
    printf("bad bc\n");
    reset_mock_count();
    se[3][0] = 0;
    assert_int_equal(calc_jacobian(&m, &J[0]), 0);
    assert_int_equal(malloc_count, 0);
    assert_int_equal(realloc_count, 0);
    se[3][0] = 3;
    /* various malloc failures */
    printf("malloc sad\n");
    for(int i = 0, mallocs = 0, reallocs = 0; i < NUM_JAC_MALLOCS_PMAP + NUM_JAC_REALLOCS_PMAP; i++) {
        reset_mock_count();
        switch(i) {
        case 11:
        case 12:
        case 13:
        case 27:
        case 28:
        case 29:
            reallocs++;
            realloc_fail_at = reallocs;
            break;
        default:
            mallocs++;
            malloc_fail_at = mallocs;
        }
        printf("[i=%d] mallocs = %d/%d, reallocs = %d/%d\n", i,
               mallocs, NUM_JAC_MALLOCS_PMAP,
               reallocs, NUM_JAC_REALLOCS_PMAP);
        assert_int_equal(calc_jacobian(&m, &J[0]), 0);
        // expected malloc & realloc
        assert_int_equal(malloc_count, mallocs);
        assert_int_equal(realloc_count, reallocs);
    }
    /* build system matrix failure: singular element volume */
    printf("singular element (zero volume)\n");
    reset_mock_count();
    elems[0][1] = 1;
    assert_int_equal(calc_jacobian(&m, &J[0]), 0);
    assert_int_equal(malloc_count, NUM_JAC_MALLOCS);
    assert_int_equal(realloc_count, NUM_JAC_REALLOCS);
    elems[0][1] = 2;
    /* passing example again */
    printf("happy\n");
    reset_mock_count();
    assert_int_equal(calc_jacobian(&m, &J[0]), 1);
    assert_int_equal(malloc_count, NUM_JAC_MALLOCS_PMAP);
    assert_int_equal(realloc_count, NUM_JAC_REALLOCS_PMAP);
    test_free(m.elec_to_sys);
    test_free(m.zc);
}


int main(void)
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_malloc_and_free),
        cmocka_unit_test(test_2d_resistor_cem),
        cmocka_unit_test(test_3d_resistor_cem),
        cmocka_unit_test(test_2d_resistor_pmap),
        cmocka_unit_test(test_fwd_solve_fail),
        cmocka_unit_test(test_calc_jacobian_fail),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
