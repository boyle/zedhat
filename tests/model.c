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
#include <stdlib.h> /* qsort */

#include <cholmod.h>

#include "cmocka.h"
#include "model.h"

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

// TODO FIXME these undef shouldn't be required
#undef malloc
#undef free

#define assert_mat_equal(m, n, a, b, delta) do { \
   for (int i=0; i < m*n; i++) { \
         assert_float_equal(a[i],b[i],delta); \
   }\
} while(0)

void test_model(void ** state)
{
    model * m;
    free_model(NULL);
    will_return(_mock_test_malloc, 0);
    m = malloc_model();
    assert_non_null(m);
    m = free_model(m);
    assert_null(m);
    will_return(_mock_test_malloc, 1);
    m = malloc_model();
    assert_null(m);
}

double simple_det3(int n, double A[3][3])
{
    assert_int_equal(n, 3);
    double d = 0;
    for (int i = 0; i < 3; i++) {
        d += A[0][i] * (A[1][(i + 1) % 3] * A[2][(i + 2) % 3] - A[2][(i + 1) % 3] * A[1][(i + 2) % 3]);
    }
    return d;
}

void test_det2(void ** state)
{
    double I[2][2] = {
        {1, 0},
        {0, 1},
    };
    double A[2][2] = {
        {1, 9},
        {4, 4},
    };
    double At[2][2] = {
        {1, 4},
        {9, 4},
    };
    double B[2][2] = {
        {0.1, 5.3},
        {1.1, 2.1},
    };
    double AB[2][2] = {{0}};
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                AB[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    printf("|I| = %g\n", det(2, I));
    printf("|A| = %g\n", det(2, A));
    printf("|At| = %g\n", det(2, At));
    printf("|B| = %g\n", det(2, B));
    printf("|AB| = %g\n", det(2, AB));
    assert_float_equal( det(2, I), 1.0, 0.0);
    assert_float_equal( det(2, A), det(2, At), 0.0);
    assert_float_equal( det(2, AB), det(2, A)*det(2, B), 2e4 * DBL_EPSILON);
    printf("checks |I| = 1; |A| = |At|; |A*B| = |A|*|B|\n");
}

void test_det3(void ** state)
{
    double I[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };
    double A[3][3] = {
        {1, 9, 3},
        {4, 4, 6},
        {7, 8, 9},
    };
    double At[3][3] = {
        {1, 4, 7},
        {9, 4, 8},
        {3, 6, 9},
    };
    double B[3][3] = {
        {0.1, 5.3, 3.3},
        {1.1, 2.1, 7.2},
        {4.1, 0.8, 9.1},
    };
    double AB[3][3] = {{0}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                AB[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    printf("|I| = %g\n", det(3, I));
    printf("|A| = %g\n", det(3, A));
    printf("|At| = %g\n", det(3, At));
    printf("|B| = %g\n", det(3, B));
    printf("|AB| = %g\n", det(3, AB));
    assert_float_equal( det(3, I), 1.0, 0.0);
    assert_float_equal( det(3, A), det(3, At), 0.0);
    assert_float_equal( det(3, AB), det(3, A)*det(3, B), 2e4 * DBL_EPSILON);
    printf("checks |I| = 1; |A| = |At|; |A*B| = |A|*|B|\n");
    /* compare to simple_det3 when n == 3 */
    assert_float_equal( det(3, I), simple_det3(3, I), 0.0);
    assert_float_equal( det(3, A), simple_det3(3, A), 0.0);
    assert_float_equal( det(3, B), simple_det3(3, B), 0.0);
    assert_float_equal( det(3, At), simple_det3(3, At), 0.0);
    assert_float_equal( det(3, AB), simple_det3(3, AB), 0.0);
}

void test_det4(void ** state)
{
    double I[4][4] = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
    };
    double A[4][4] = {
        {1, 2, 3, 4},
        {4, 4, 6, 3},
        {7, 8, 2, 9},
        {1, 4, 7, 3},
    };
    double At[4][4] = {
        {1, 4, 7, 1},
        {2, 4, 8, 4},
        {3, 6, 2, 7},
        {4, 3, 9, 3},
    };
    double B[4][4] = {
        {0.1, 1.3, 3.3, 0.1},
        {1.1, 2.1, 0.3, 0.5},
        {4.1, 0.8, 0.2, 1.2},
        {2.2, 2.5, 0.6, 0.7},
    };
    double AB[4][4] = {{0}};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                AB[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    printf("|I| = %g\n", det(4, I));
    printf("|A| = %0.16g\n", det(4, A));
    printf("|At| = %g\n", det(4, At));
    printf("|B| = %0.16g\n", det(4, B));
    printf("|AB| = %0.16g\n", det(4, AB));
    assert_float_equal( det(4, I), 1.0, 0.0);
    assert_float_equal( det(4, A), det(4, At), 0.0);
    assert_float_equal( det(4, AB), det(4, A)*det(4, B), 2e5 * DBL_EPSILON);
    printf("checks |I| = 1; |A| = |At|; |A*B| = |A|*|B|\n");
}

void printf_mat_double(char * varname, int n, int m, double * mat)
{
    printf("%s =\n", varname);
    const int space = 3;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < space; j++) {
            printf(" ");
        }
        printf("[");
        for (int j = 0; j < m; j++) {
            printf(" % 20.15g", mat[i * m + j]);
        }
        printf(" ]\n");
    }
}

void printf_mat_int(char * varname, int n, int m, int * mat)
{
    printf("%s =\n", varname);
    const int space = 3;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < space; j++) {
            printf(" ");
        }
        printf("[");
        for (int j = 0; j < m; j++) {
            printf(" % 10d", mat[i * m + j]);
        }
        printf(" ]\n");
    }
}

double * __simple_inv2(int n, double A[n][n])
{
    assert_int_equal(n, 2);
    const double a = A[0][0];
    const double b = A[0][1];
    const double c = A[1][0];
    const double d = A[1][1];
    const double detA = a * d - b * c;
    const double B[2][2] = {{+d, -b}, {-c, +a}};
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++) {
            A[i][j] = B[i][j] / detA;
        }
    return &(A[0][0]);
}

double * __simple_inv3(int n, double A[n][n])
{
    assert_int_equal(n, 3);
    const double a = A[0][0];
    const double b = A[0][1];
    const double c = A[0][2];
    const double d = A[1][0];
    const double e = A[1][1];
    const double f = A[1][2];
    const double g = A[2][0];
    const double h = A[2][1];
    const double i = A[2][2];
    const double B[3][3] = {
        {e * i - f * h, -(b * i - c * h), b * f - c * e},
        {-(d * i - f * g), a * i - c * g, -(a * f - c * d)},
        {d * h - e * g, -(a * h - b * g), a * e - b * d},
    };
    const double detA = a * B[0][0] + b * B[1][0] + c * B[2][0];
    for(int k = 0; k < n; k++)
        for(int j = 0; j < n; j++) {
            A[k][j] = B[k][j] / detA;
        }
    return &(A[0][0]);
}

double * mat_wrap(int n, double in[n][n], double out[n][n], double * (*func)(int n, double A[n][n]))
{
    memcpy(&(out[0][0]), &(in[0][0]), pow(n, 2)*sizeof(double));
    return func(n, out);
}

double * simple_inv2(int n, double in[n][n], double out[n][n])
{
    return mat_wrap(n, in, out, &(__simple_inv2)) ;
}

double * simple_inv3(int n, double in[n][n], double out[n][n])
{
    return mat_wrap(n, in, out, &(__simple_inv3)) ;
}

void test_inv2(void ** state)
{
    double I[2][2] = {
        {1, 0},
        {0, 1},
    };
    double A[2][2] = {
        {2, 0},
        {0, 4},
    };
    double Ai[2][2] = {
        {0.5, 0},
        {0, 0.25},
    };
    double B[2][2] = {
        {0.1, 5.3},
        {1.1, 2.1},
    };
    double tmp[2][2] = {{0}};
    double * expected = NULL;
    printf_mat_double("I^-1", 2, 2, simple_inv2(2, I, tmp));
    printf_mat_double("A^-1", 2, 2, simple_inv2(2, A, tmp));
    printf_mat_double("Ai", 2, 2, &(Ai[0][0]));
    printf_mat_double("B^-1", 2, 2, simple_inv2(2, B, tmp));
    expected = &(I[0][0]);
    assert_mat_equal( 2, 2, simple_inv2(2, I, tmp), expected, 0.0);
    expected = &(Ai[0][0]);
    assert_mat_equal( 2, 2, simple_inv2(2, A, tmp), expected, 0.0);
    double Bi[2][2] = {{0}};
    double BiB[2][2] = {{0}};
    simple_inv2(2, B, Bi);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                BiB[i][j] += Bi[i][k] * B[k][j];
            }
        }
    }
    expected = &(I[0][0]);
    double * in = &(BiB[0][0]);
    assert_mat_equal( 2, 2, in, expected, DBL_EPSILON);
    printf("checks I^-1 = I; A^-1 = Ai for diag A; B^-1 * B = I\n");
    double tt[2][2] = {{0}};
    assert_mat_equal( 2, 2, mat_wrap(2, I, tmp, inv), simple_inv2(2, I, tt), 0);
    assert_mat_equal( 2, 2, mat_wrap(2, A, tmp, inv), simple_inv2(2, A, tt), 0);
    assert_mat_equal( 2, 2, mat_wrap(2, Ai, tmp, inv), simple_inv2(2, Ai, tt), 0);
    assert_mat_equal( 2, 2, mat_wrap(2, B, tmp, inv), simple_inv2(2, B, tt), DBL_EPSILON);
}

void test_inv3(void ** state)
{
    double I[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };
    double A[3][3] = {
        {2, 0, 0},
        {0, 4, 0},
        {0, 0, 8},
    };
    double Ai[3][3] = {
        {1.0 / 2.0, 0, 0},
        {0, 1.0 / 4.0, 0},
        {0, 0, 1.0 / 8.0},
    };
    double B[3][3] = {
        {0.1, 5.3, 3.2},
        {1.1, 2.1, 4.1},
        {2.1, 3.5, 0.3},
    };
    double tmp[3][3] = {{0}};
    double * expected = NULL;
    printf_mat_double("I^-1", 3, 3, simple_inv3(3, I, tmp));
    printf_mat_double("A^-1", 3, 3, simple_inv3(3, A, tmp));
    printf_mat_double("Ai", 3, 3, &(Ai[0][0]));
    printf_mat_double("B^-1", 3, 3, simple_inv3(3, B, tmp));
    expected = &(I[0][0]);
    assert_mat_equal( 3, 3, simple_inv3(3, I, tmp), expected, 0.0);
    expected = &(Ai[0][0]);
    assert_mat_equal( 3, 3, simple_inv3(3, A, tmp), expected, 0.0);
    double Bi[3][3] = {{0}};
    double BiB[3][3] = {{0}};
    simple_inv3(3, B, Bi);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                BiB[i][j] += Bi[i][k] * B[k][j];
            }
        }
    }
    expected = &(I[0][0]);
    double * in = &(BiB[0][0]);
    printf_mat_double("B^-1*B = I? ...", 3, 3, &(BiB[0][0]));
    assert_mat_equal( 3, 3, in, expected, 5 * DBL_EPSILON);
    printf("checks I^-1 = I; A^-1 = Ai for diag A; B^-1 * B = I\n");
    double tt[3][3] = {{0}};
    assert_mat_equal( 3, 3, mat_wrap(3, I, tmp, inv), simple_inv3(3, I, tt), 0);
    assert_mat_equal( 3, 3, mat_wrap(3, A, tmp, inv), simple_inv3(3, A, tt), 0);
    assert_mat_equal( 3, 3, mat_wrap(3, Ai, tmp, inv), simple_inv3(3, Ai, tt), 0);
    assert_mat_equal( 3, 3, mat_wrap(3, B, tmp, inv), simple_inv3(3, B, tt), DBL_EPSILON);
}

void test_inv4(void ** state)
{
    double I[4][4] = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
    };
    double A[4][4] = {
        {2, 0, 0, 0},
        {0, 4, 0, 0},
        {0, 0, 8, 0},
        {0, 0, 0, 16},
    };
    double Ai[4][4] = {
        {1.0 / 2.0, 0, 0, 0},
        {0, 1.0 / 4.0, 0, 0},
        {0, 0, 1.0 / 8.0, 0},
        {0, 0, 0, 1.0 / 16.0},
    };
    double B[4][4] = {
        {0.1, 5.3, 3.2, 0.1},
        {1.1, 2.1, 4.1, 0.2},
        {2.1, 3.5, 0.3, 0.3},
        {5.1, 1.2, 3.7, 8.2},
    };
    double tmp[4][4] = {{0}};
    double * expected = NULL;
    printf_mat_double("I^-1", 4, 4, mat_wrap(4, I, tmp, inv));
    printf_mat_double("A^-1", 4, 4, mat_wrap(4, A, tmp, inv));
    printf_mat_double("Ai", 4, 4, &(Ai[0][0]));
    printf_mat_double("B^-1", 4, 4, mat_wrap(4, B, tmp, inv));
    expected = &(I[0][0]);
    assert_mat_equal( 4, 4, mat_wrap(4, I, tmp, inv), expected, 0.0);
    expected = &(Ai[0][0]);
    assert_mat_equal( 4, 4, mat_wrap(4, A, tmp, inv), expected, 0.0);
    double Bi[4][4] = {{0}};
    double BiB[4][4] = {{0}};
    mat_wrap(4, B, Bi, inv);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                BiB[i][j] += Bi[i][k] * B[k][j];
            }
        }
    }
    expected = &(I[0][0]);
    double * in = &(BiB[0][0]);
    printf_mat_double("B^-1*B = I? ...", 3, 3, &(BiB[0][0]));
    assert_mat_equal( 4, 4, in, expected, 5 * DBL_EPSILON);
    printf("checks I^-1 = I; A^-1 = Ai for diag A; B^-1 * B = I\n");
}

void test_bad_inv(void ** state)
{
    double singularA[2][2] = {
        {0, 0},
        {0, 1},
    };
    assert_true( inv(2, singularA) == NULL );
}

double * sym_to_full(int n, double * As, double * Bf)
{
    memset(Bf, 0, sizeof(double)*n * n);
    int m = 0;
    printf("sym = ");
    for (int i = 0; i < n; i++ ) { /* row-major */
        for (int j = i; j < n; j++ ) {
            Bf[i + j * n] = As[m];
            Bf[i * n + j] = As[m];
            printf(" %0.2f", As[m]);
            m++;
        }
    }
    printf(" (row-major, upper-triangular)\n");
    return Bf;
}

void test_shape_Se_v(void ** state)
{
    int ii [4 * 4];
    int jj [4 * 4];
    double ss[4 * 4];
    double out[4 * 4];
    int elems[4] = {1, 2, 3, 4};
    mesh m = { 0 };
    m.dim = 2;
    m.elems = &(elems[0]);
    m.n_elems = 1;
    m.n_nodes = 3;
    {
        double singularA[3][2] = {
            {0, 0},
            {0, 0},
            {1, 0},
        };
        printf("Se: check for singular matrix erroring out\n");
        m.nodes = &(singularA[0][0]);
        assert_int_equal( calc_Se(&m, ii, jj, ss), 1 );
    }
    /* matlab:
     *   E = inv([1 1 0; 1 0 1; 1 0 0])
     *   Sh = E(2:3,:)*sqrt(1/2/abs(det(E)))
     *   Se = Sh'*Sh
     */
    {
        double A[3][2] = {
            {1.0, 0.0},
            {0.0, 1.0},
            {0.0, 0.0},
        };
        m.nodes = &(A[0][0]);
        assert_int_equal( calc_Se(&m, ii, jj, ss), 0 );
        printf_mat_double("A", 3, 2, &(A[0][0]));
        printf_mat_double("A --> S", 3, 3, sym_to_full(3, ss, out));
        double Se2[3][3] = {
            { 0.5,  0.0, -0.5},
            { 0.0,  0.5, -0.5},
            {-0.5, -0.5,  1.0},
        };
        printf_mat_double("Se", 3, 3, &(Se2[0][0]));
        double * tmp2Se = &(Se2[0][0]);
        assert_mat_equal( 3, 3, out, tmp2Se, DBL_EPSILON);
    }
    /* matlab:
     *   E = inv([1 1 0 0; 1 0 1 0; 1 0 0 1; 1 0 0 0])
     *   Sh = E(2:4,:)*sqrt(1/6/abs(det(E)))
     *   S = Sh'*Sh
     */
    {
        double B[4][3] = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {0.0, 0.0, 0.0},
        };
        m.dim = 3;
        m.nodes = &(B[0][0]);
        assert_int_equal( calc_Se(&m, ii, jj, ss), 0 );
        printf_mat_double("B", 4, 3, &(B[0][0]));
        printf_mat_double("B --> S", 4, 4, sym_to_full(4, ss, out));
        double Se3[4][4] = {
            {1.0 / 6.0, 0, 0, -1.0 / 6.0},
            {0, 1.0 / 6.0, 0, -1.0 / 6.0},
            {0, 0, 1.0 / 6.0, -1.0 / 6.0},
            {-1.0 / 6.0, -1.0 / 6.0, -1.0 / 6.0, 1.0 / 2.0},
        };
        printf_mat_double("Se", 4, 4, &(Se3[0][0]));
        double * tmp3Se = &(Se3[0][0]);
        assert_mat_equal( 4, 4, out, tmp3Se, DBL_EPSILON);
    }
}

void test_shape_Se_n(void ** state)
{
    assert_true(calc_Se_n(2) == 6);
    assert_true(calc_Se_n(3) == 10);
}

void printf_ii(const char * name, int n, int * ii)
{
    printf("%5s = %d", name, ii[0]);
    for(int j = 1; j < n; j++) {
        printf(", %d", ii[j]);
    }
    printf("\n");
}

void test_shape_Se_ij(void ** state)
{
    const int n3 = calc_Se_n(3);
    int ii [n3 + 1];
    int jj [n3 + 1];
    double ss [n3 + 1];
    double nodes[4][3] = {
        {1, 0, 0},
        {1, 1, 0},
        {0, 1, 1},
        {0, 0, 0},
    };
    int elems[4] = {1, 2, 3, 4};
    mesh m = { 0 };
    m.elems = &(elems[0]);
    m.nodes = &(nodes[0][0]);
    m.n_elems = 1;
    m.n_nodes = 4;
    for(int nd = 2; nd <= 3; nd++) {
        m.dim = nd;
        const int n = calc_Se_n(nd);
        ii[n] = -1;
        jj[n] = -1;
        printf_ii("elems", nd + 1, elems);
        calc_Se(&m, ii, jj, ss);
        printf_ii("ii", n, ii);
        printf_ii("jj", n, jj);
        int k = 0;
        for(int i = 0; i < nd + 1; i++ ) {
            for(int j = i; j < nd + 1; j++ ) {
                assert_int_equal(ii[k], i);
                assert_int_equal(jj[k], j);
                k++;
            }
        }
        /* check guard values at end of array */
        assert_int_equal(ii[k], -1);
        assert_int_equal(jj[k], -1);
    }
}

int cmp_int_tests( const void * a, const void * b )
{
    return *(int *)a - *(int *)b;
}

void test_shape_2d(void ** state)
{
    /* 2D square:
     *    first four nodes at z=0
     *    first two elems */
    double nodes[5][2] = {
        {0, 0},
        {1, 0},
        {0, 1},
        {1, 1},
        {0, 0}, /* guard */
    };
    printf_mat_double("nodes", 4, 2, &(nodes[0][0]));
    /* TODO currently each row of elems must be sorted
     * or we won't get an upper triangular matrix */
    /* Our 'elems' starts at 1 which agrees with netgen */
    int elems[4][3] = { /* from netgen cube.geo */
        {0, 0, 0}, /* guard */
        {2, 1, 4},
        {1, 3, 4},
        {0, 0, 0}, /* guard */
    };
    int * e = &(elems[1][0]);
    printf_mat_int("elems", 2, 3, e);
    for(int i = 0; i < 2; i++) {
        qsort(&(e[i * 3]), 3, sizeof(int), &cmp_int_tests);
    }
    printf_mat_int("elems (sorted)", 2, 3, e);
    const int nnz = calc_Se_n(2) * 2;
    int * ii = malloc(sizeof(int) * nnz);
    int * jj = malloc(sizeof(int) * nnz);
    double * ss = malloc(sizeof(double) * nnz);
    mesh m = { 0 };
    m.dim = 2;
    m.elems = &(elems[1][0]);
    m.nodes = &(nodes[0][0]);
    m.n_elems = 2;
    m.n_nodes = 4;
    for(int gnd = 0; gnd < 7; gnd++ ) {
        int ret = calc_Se(&m, ii, jj, ss);
        assert_int_equal(ret, 0);
        int nnz_local = nnz;
        int ngnd = calc_gnd(gnd, &nnz_local, ii, jj, ss);
        printf("calc_gnd() gnd node #%d: deleted %d entries\n", gnd, ngnd);
        if (gnd > 0 && gnd <= 4 ) {
            assert_int_not_equal(ngnd, 0);
        }
        else {
            assert_int_equal(ngnd, 0);
        }
    }
    /* intentionally try building a bad mesh */
    {
        m.elems = &(elems[0][0]);
        int ret = calc_Se(&m, ii, jj, ss);
        assert_int_equal(ret, 1);
    }
    {
        m.elems = &(elems[2][0]);
        int ret = calc_Se(&m, ii, jj, ss);
        assert_int_equal(ret, 2);
    }
    /* clean up */
    free(ii);
    free(jj);
    free(ss);
}

void test_shape_3d(void ** state)
{
    /* 3D cube:
     *    all nodes used as in elems */
    double nodes[9][3] = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 1},
        {0, 0, 0}, /* guard */
    };
    printf_mat_double("nodes", 8, 3, &(nodes[0][0]));
    /* TODO currently each row of elems must be sorted
     * or we won't get an upper triangular matrix */
    /* Our 'elems' starts at 1 which agrees with netgen */
    int elems[8][4] = { /* from netgen cube.geo */
        {0, 0, 0, 0}, /* guard */
        {2, 5, 1, 4},
        {1, 5, 3, 4},
        {3, 5, 7, 8},
        {8, 4, 5, 6},
        {2, 6, 5, 4},
        {8, 4, 3, 5},
        {0, 0, 0, 0}, /* guard */
    };
    int * e = &(elems[1][0]);
    printf_mat_int("elems", 6, 4, e);
    for(int i = 0; i < 6; i++) {
        qsort(&(e[i * 4]), 4, sizeof(int), &cmp_int_tests);
    }
    printf_mat_int("elems (sorted)", 6, 4, e);
    const int nnz = calc_Se_n(3) * 6;
    int * ii = malloc(sizeof(int) * nnz);
    int * jj = malloc(sizeof(int) * nnz);
    double * ss = malloc(sizeof(double) * nnz);
    mesh m = { 0 };
    m.dim = 3;
    m.elems = &(elems[1][0]);
    m.nodes = &(nodes[0][0]);
    m.n_elems = 6;
    m.n_nodes = 8;
    for(int gnd = 0; gnd < 11; gnd++ ) {
        int ret = calc_Se(&m, ii, jj, ss);
        assert_int_equal(ret, 0);
        int nnz_local = nnz;
        int ngnd = calc_gnd(gnd, &nnz_local, ii, jj, ss);
        printf("calc_gnd() gnd node #%d: deleted %d entries\n", gnd, ngnd);
        if (gnd > 0 && gnd <= 8 ) {
            assert_int_not_equal(ngnd, 0);
        }
        else {
            assert_int_equal(ngnd, 0);
        }
    }
    /* intentionally try building a bad mesh */
    {
        m.elems = &(elems[0][0]);
        int ret = calc_Se(&m, ii, jj, ss);
        assert_int_equal(ret, 1);
    }
    {
        m.elems = &(elems[2][0]);
        int ret = calc_Se(&m, ii, jj, ss);
        assert_int_equal(ret, 6);
    }
    /* clean up */
    free(ii);
    free(jj);
    free(ss);
}

void test_bc_2d_neumann (void ** state)
{
    int bc[4] = {0, 2, 0, 1};
    int se[4][2] = {
        {1, 2},
        {2, 4},
        {4, 3},
        {3, 1},
    };
    mesh m = {0};
    m.dim = 2;
    m.bc = &(bc[0]);
    m.surfaceelems = &(se[0][0]);
    m.n_se = 4;
    m.n_nodes = 4;
    for (int gnd = 0; gnd <= 4; gnd++) {
        double b[4] = {0}; /* nodes */
        int ret = calc_stim_neumann(&m, +1, 1, gnd, &(b[0]));
        ret += calc_stim_neumann(&m, -1, 2, gnd, &(b[0]));
        double expect[4] = { +0.5,
                             -0.5,
                             +0.5,
                             -0.5,
                           };
        if(gnd > 0) {
            expect[gnd - 1] = 0;
            assert_int_equal(ret, 3);
        }
        else {
            assert_int_equal(ret, 4);
        }
        printf("gnd = %d\n", gnd);
        printf_mat_double("b", 4, 1, &(b[0]));
        assert_mat_equal( 4, 1, b, expect, DBL_EPSILON);
    }
}

void test_bc_3d_neumann (void ** state)
{
    int bc[12] = {1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6};
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
    mesh m = {0};
    m.dim = 3;
    m.bc = &(bc[0]);
    m.surfaceelems = &(se[0][0]);
    m.n_se = 12;
    m.n_nodes = 8;
    for (int gnd = 0; gnd <= 8; gnd++) {
        double b[8] = {0}; /* nodes */
        int ret = calc_stim_neumann(&m, +1, 1, gnd, &(b[0]));
        ret += calc_stim_neumann(&m, -1, 2, gnd, &(b[0]));
        /* TODO these are not right... I is the same for all boundary edge nodes */
        double expect[8] = { +1.0 / 6.0,
                             +2.0 / 6.0,
                             -2.0 / 6.0,
                             -1.0 / 6.0,
                             +2.0 / 6.0,
                             +1.0 / 6.0,
                             -1.0 / 6.0,
                             -2.0 / 6.0,
                           };
        if(gnd > 0) {
            expect[gnd - 1] = 0;
        }
        int expect_ret;
        switch(gnd) {
        case 0:
            expect_ret = 12;
            break;
        case 1:
        case 4:
        case 6:
        case 7:
            expect_ret = 11;
            break;
        default:
            expect_ret = 10;
        }
        printf("gnd = %d, expect %d local nodes\n", gnd, expect_ret);
        printf_mat_double("b", 8, 1, &(b[0]));
        assert_mat_equal( 8, 1, b, expect, 10 * DBL_EPSILON);
        assert_int_equal(ret, expect_ret);
    }
}

void test_2d_resistor (void ** state)
{
    /* 2D square:
     *    first four nodes at z=0
     *    first two elems */
    double nodes[5][2] = {
        {0, 0},
        {1, 0},
        {0, 1},
        {1, 1},
        {0, 0}, /* guard */
    };
    int elems[4][3] = { /* from netgen cube.geo */
        {0, 0, 0}, /* guard */
        {1, 2, 4},
        {1, 3, 4},
        {0, 0, 0}, /* guard */
    };
    int bc[4] = {0, 1, 0, 2};
    int se[4][2] = {
        {1, 2},
        {2, 4},
        {4, 3},
        {3, 1},
    };
    mesh m = { 0 };
    m.dim = 2;
    m.elems = &(elems[1][0]);
    m.nodes = &(nodes[0][0]);
    m.n_elems = 2;
    m.n_nodes = 4;
    m.bc = &(bc[0]);
    m.surfaceelems = &(se[0][0]);
    m.n_se = 4;
    const int gnd = 1;
    int nnz = calc_Se_n(2) * m.n_elems;
    int * ii = malloc(sizeof(int) * nnz);
    int * jj = malloc(sizeof(int) * nnz);
    double * ss = malloc(sizeof(double) * nnz);
    int ret = calc_Se(&m, ii, jj, ss);
    assert_int_equal(ret, 0);
    int ngnd = calc_gnd(gnd, &nnz, ii, jj, ss);
    assert_int_not_equal(ngnd, 0);
    double bb[4] = {0}; /* nodes */
    ret = calc_stim_neumann(&m, +1, 1, gnd, &(bb[0]));
    ret += calc_stim_neumann(&m, -1, 2, gnd, &(bb[0]));
    printf_mat_double("bb", 4, 1, &(bb[0]));
    /* fwd_solve */
    cholmod_dense * x, *b, *r ;
    cholmod_factor * L ;
    double one [2] = {1, 0}, m1 [2] = {-1, 0} ;     /* basic scalars */
    cholmod_common c ;
    cholmod_start (&c) ;                /* start CHOLMOD */
    cholmod_triplet Acoo = {
        m.n_nodes, m.n_nodes, nnz, nnz,
        ii, jj, ss, NULL, +1,
        CHOLMOD_INT, CHOLMOD_REAL, CHOLMOD_DOUBLE
    };
    cholmod_sparse * A = cholmod_triplet_to_sparse(&Acoo, nnz, &c);
    assert_ptr_not_equal(A, NULL); /* if allocate failed */
    assert_int_equal(cholmod_check_sparse (A, &c), 1); /* okay matrix? TRUE=1/FALSE=0 */
    cholmod_print_sparse (A, "A", &c) ;         /* print the matrix */
    assert_int_not_equal(A->stype, 0); /* A must be symmetric */
    b = cholmod_ones (A->nrow, 1, A->xtype, &c) ;   /* b = ones(n,1) */
    /* set b->x = bb */
    for (int i = 0; i < m.n_nodes; i++) {
        ((double *) b->x)[i] = bb[i];
    }
    L = cholmod_analyze (A, &c) ;           /* analyze */
    cholmod_factorize (A, L, &c) ;          /* factorize */
    x = cholmod_solve (CHOLMOD_A, L, b, &c) ;       /* solve Ax=b */
    r = cholmod_copy_dense (b, &c) ;            /* r = b */
    cholmod_sdmult (A, 0, m1, one, x, r, &c) ;      /* r = r-Ax */
    double norm = cholmod_norm_dense (r, 0, &c);
    printf ("norm(b-Ax) = %8.1e\n",
            cholmod_norm_dense (r, 0, &c)) ;        /* print norm(r) */
    assert_float_equal(norm, 0.0, DBL_EPSILON);
    double * soln = x->x;
    printf_mat_double("x", 4, 1, soln);
    /* check x */
    double expect[4] = { 0, +1, 0, +1 };
    assert_mat_equal( m.n_nodes, 1, soln, expect, 10 * DBL_EPSILON);
    /* clean up */
    cholmod_free_factor (&L, &c) ;          /* free matrices */
    cholmod_free_sparse (&A, &c) ;
    cholmod_free_dense (&r, &c) ;
    cholmod_free_dense (&x, &c) ;
    cholmod_free_dense (&b, &c) ;
    cholmod_finish (&c) ;               /* finish CHOLMOD */
    /* clean up */
    free(ii);
    free(jj);
    free(ss);
}

void test_3d_resistor (void ** state)
{
    /* 3D cube:
     *    all nodes used as in elems */
    double nodes[9][3] = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 1},
        {0, 0, 0}, /* guard */
    };
    printf_mat_double("nodes", 8, 3, &(nodes[0][0]));
    /* TODO currently each row of elems must be sorted
     * or we won't get an upper triangular matrix */
    /* Our 'elems' starts at 1 which agrees with netgen */
    int elems[8][4] = { /* from netgen cube.geo */
        {0, 0, 0, 0}, /* guard */
        {2, 5, 1, 4},
        {1, 5, 3, 4},
        {3, 5, 7, 8},
        {8, 4, 5, 6},
        {2, 6, 5, 4},
        {8, 4, 3, 5},
        {0, 0, 0, 0}, /* guard */
    };
    int * e = &(elems[1][0]);
    printf_mat_int("elems", 6, 4, e);
    for(int i = 0; i < 6; i++) {
        qsort(&(e[i * 4]), 4, sizeof(int), &cmp_int_tests);
    }
    printf_mat_int("elems (sorted)", 6, 4, e);
    int bc[12] = {2, 2, 1, 1, 3, 3, 4, 4, 5, 5, 6, 6};
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
    mesh m = { 0 };
    m.dim = 3;
    m.elems = &(elems[1][0]);
    m.nodes = &(nodes[0][0]);
    m.n_elems = 6;
    m.n_nodes = 8;
    m.bc = &(bc[0]);
    m.surfaceelems = &(se[0][0]);
    m.n_se = 12;
    const int gnd = 1;
    int nnz = calc_Se_n(3) * m.n_elems;
    int * ii = malloc(sizeof(int) * nnz);
    int * jj = malloc(sizeof(int) * nnz);
    double * ss = malloc(sizeof(double) * nnz);
    int ret = calc_Se(&m, ii, jj, ss);
    assert_int_equal(ret, 0);
    int ngnd = calc_gnd(gnd, &nnz, ii, jj, ss);
    assert_int_not_equal(ngnd, 0);
    double bb[8] = {0}; /* nodes */
    ret = calc_stim_neumann(&m, +1, 1, gnd, &(bb[0]));
    ret += calc_stim_neumann(&m, -1, 2, gnd, &(bb[0]));
    printf_mat_double("bb", 8, 1, &(bb[0]));
    /* fwd_solve */
    cholmod_dense * x, *b, *r ;
    cholmod_factor * L ;
    double one [2] = {1, 0}, m1 [2] = {-1, 0} ;     /* basic scalars */
    cholmod_common c ;
    cholmod_start (&c) ;                /* start CHOLMOD */
    cholmod_triplet Acoo = {
        m.n_nodes, m.n_nodes, nnz, nnz,
        ii, jj, ss, NULL, +1,
        CHOLMOD_INT, CHOLMOD_REAL, CHOLMOD_DOUBLE
    };
    cholmod_sparse * A = cholmod_triplet_to_sparse(&Acoo, nnz, &c);
    assert_ptr_not_equal(A, NULL); /* if allocate failed */
    assert_int_equal(cholmod_check_sparse (A, &c), 1); /* okay matrix? TRUE=1/FALSE=0 */
    cholmod_print_sparse (A, "A", &c) ;         /* print the matrix */
    assert_int_not_equal(A->stype, 0); /* A must be symmetric */
    b = cholmod_ones (A->nrow, 1, A->xtype, &c) ;   /* b = ones(n,1) */
    /* set b->x = bb */
    for (int i = 0; i < m.n_nodes; i++) {
        ((double *) b->x)[i] = bb[i];
    }
    L = cholmod_analyze (A, &c) ;           /* analyze */
    cholmod_factorize (A, L, &c) ;          /* factorize */
    x = cholmod_solve (CHOLMOD_A, L, b, &c) ;       /* solve Ax=b */
    r = cholmod_copy_dense (b, &c) ;            /* r = b */
    cholmod_sdmult (A, 0, m1, one, x, r, &c) ;      /* r = r-Ax */
    double norm = cholmod_norm_dense (r, 0, &c);
    printf ("norm(b-Ax) = %8.1e\n",
            norm) ;        /* print norm(r) */
    assert_float_equal(norm, 0.0, 2 * DBL_EPSILON);
    double * soln = x->x;
    printf_mat_double("x", m.n_nodes, 1, soln);
    /* check x */
    double expect[8] = { 0, 0, +1, +1, 0, 0, +1, +1 };
    assert_mat_equal( m.n_nodes, 1, soln, expect, 10 * DBL_EPSILON);
    /* clean up */
    cholmod_free_factor (&L, &c) ;          /* free matrices */
    cholmod_free_sparse (&A, &c) ;
    cholmod_free_dense (&r, &c) ;
    cholmod_free_dense (&x, &c) ;
    cholmod_free_dense (&b, &c) ;
    cholmod_finish (&c) ;               /* finish CHOLMOD */
    /* clean up */
    free(ii);
    free(jj);
    free(ss);
}

int main(void)
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_model),
        cmocka_unit_test(test_det2),
        cmocka_unit_test(test_det3),
        cmocka_unit_test(test_det4),
        cmocka_unit_test(test_inv2),
        cmocka_unit_test(test_inv3),
        cmocka_unit_test(test_inv4),
        cmocka_unit_test(test_bad_inv),
        cmocka_unit_test(test_shape_Se_n),
        cmocka_unit_test(test_shape_Se_v),
        cmocka_unit_test(test_shape_Se_ij),
        cmocka_unit_test(test_shape_2d),
        cmocka_unit_test(test_shape_3d),
        cmocka_unit_test(test_bc_2d_neumann),
        cmocka_unit_test(test_bc_3d_neumann),
        cmocka_unit_test(test_2d_resistor),
        cmocka_unit_test(test_3d_resistor),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
