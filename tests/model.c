/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#include <stdarg.h> /* cmocka.h */
#include <stddef.h> /* cmocka.h */
#include <setjmp.h> /* cmocka.h */

#include <stdio.h> /* cmocka.h */
#include <string.h> /* cmocka.h, memcpy */
#include <float.h> /* DBL_EPSILON */
#include <math.h> /* fabs, pow */

#include "cmocka.h"
#include "model.h"

#define assert_double_equal(a, b, delta) do { \
   if(!(fabs(a - b) <= delta)) { \
      printf("|a:%0.16g - b:%0.16g| = %0.16g > %0.16g\n", a, b, fabs(a-b), delta); \
   } \
} while(0)

#define assert_mat_equal(n, a, b, delta) do { \
   int i; \
   for (i=0; i<n*n; i++) { \
         assert_true(fabs(a[i] - b[i]) <= delta); \
   }\
} while(0)

double simple_det3(int n, double A[3][3])
{
    assert_int_equal(n, 3);
    int i;
    double d = 0;
    for (i = 0; i < 3; i++) {
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
    int i, j, k;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            for (k = 0; k < 2; k++) {
                AB[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    printf("|I| = %g\n", det(2, I));
    printf("|A| = %g\n", det(2, A));
    printf("|At| = %g\n", det(2, At));
    printf("|B| = %g\n", det(2, B));
    printf("|AB| = %g\n", det(2, AB));
    assert_double_equal( det(2, I), 1.0, 0.0);
    assert_double_equal( det(2, A), det(2, At), 0.0);
    assert_double_equal( det(2, AB), det(2, A)*det(2, B), 2e4 * DBL_EPSILON);
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
    int i, j, k;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                AB[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    printf("|I| = %g\n", det(3, I));
    printf("|A| = %g\n", det(3, A));
    printf("|At| = %g\n", det(3, At));
    printf("|B| = %g\n", det(3, B));
    printf("|AB| = %g\n", det(3, AB));
    assert_double_equal( det(3, I), 1.0, 0.0);
    assert_double_equal( det(3, A), det(3, At), 0.0);
    assert_double_equal( det(3, AB), det(3, A)*det(3, B), 2e4 * DBL_EPSILON);
    printf("checks |I| = 1; |A| = |At|; |A*B| = |A|*|B|\n");
    /* compare to simple_det3 when n == 3 */
    assert_double_equal( det(3, I), simple_det3(3, I), 0.0);
    assert_double_equal( det(3, A), simple_det3(3, A), 0.0);
    assert_double_equal( det(3, B), simple_det3(3, B), 0.0);
    assert_double_equal( det(3, At), simple_det3(3, At), 0.0);
    assert_double_equal( det(3, AB), simple_det3(3, AB), 0.0);
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
    int i, j, k;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {
                AB[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    printf("|I| = %g\n", det(4, I));
    printf("|A| = %0.16g\n", det(4, A));
    printf("|At| = %g\n", det(4, At));
    printf("|B| = %0.16g\n", det(4, B));
    printf("|AB| = %0.16g\n", det(4, AB));
    assert_double_equal( det(4, I), 1.0, 0.0);
    assert_double_equal( det(4, A), det(4, At), 0.0);
    assert_double_equal( det(4, AB), det(4, A)*det(4, B), 2e5 * DBL_EPSILON);
    printf("checks |I| = 1; |A| = |At|; |A*B| = |A|*|B|\n");
}

void printf_mat(char * varname, int n, double * m)
{
    printf("%s =\n", varname);
    const int space = 3;
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < space; j++) {
            printf(" ");
        }
        printf("[");
        for (j = 0; j < n; j++) {
            printf(" % 20.15g", m[i * n + j]);
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
    int i, j;
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++) {
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
    int k, j;
    for(k = 0; k < n; k++)
        for(j = 0; j < n; j++) {
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
    printf_mat("I^-1", 2, simple_inv2(2, I, tmp));
    printf_mat("A^-1", 2, simple_inv2(2, A, tmp));
    printf_mat("Ai", 2, &(Ai[0][0]));
    printf_mat("B^-1", 2, simple_inv2(2, B, tmp));
    expected = &(I[0][0]);
    assert_mat_equal( 2, simple_inv2(2, I, tmp), expected, 0.0);
    expected = &(Ai[0][0]);
    assert_mat_equal( 2, simple_inv2(2, A, tmp), expected, 0.0);
    double Bi[2][2] = {{0}};
    double BiB[2][2] = {{0}};
    simple_inv2(2, B, Bi);
    int i, j, k;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            for (k = 0; k < 2; k++) {
                BiB[i][j] += Bi[i][k] * B[k][j];
            }
        }
    }
    expected = &(I[0][0]);
    double * in = &(BiB[0][0]);
    assert_mat_equal( 2, in, expected, DBL_EPSILON);
    printf("checks I^-1 = I; A^-1 = Ai for diag A; B^-1 * B = I\n");
    double tt[2][2] = {{0}};
    assert_mat_equal( 2, mat_wrap(2, I, tmp, inv), simple_inv2(2, I, tt), 0);
    assert_mat_equal( 2, mat_wrap(2, A, tmp, inv), simple_inv2(2, A, tt), 0);
    assert_mat_equal( 2, mat_wrap(2, Ai, tmp, inv), simple_inv2(2, Ai, tt), 0);
    assert_mat_equal( 2, mat_wrap(2, B, tmp, inv), simple_inv2(2, B, tt), DBL_EPSILON);
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
    printf_mat("I^-1", 3, simple_inv3(3, I, tmp));
    printf_mat("A^-1", 3, simple_inv3(3, A, tmp));
    printf_mat("Ai", 3, &(Ai[0][0]));
    printf_mat("B^-1", 3, simple_inv3(3, B, tmp));
    expected = &(I[0][0]);
    assert_mat_equal( 3, simple_inv3(3, I, tmp), expected, 0.0);
    expected = &(Ai[0][0]);
    assert_mat_equal( 3, simple_inv3(3, A, tmp), expected, 0.0);
    double Bi[3][3] = {{0}};
    double BiB[3][3] = {{0}};
    simple_inv3(3, B, Bi);
    int i, j, k;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                BiB[i][j] += Bi[i][k] * B[k][j];
            }
        }
    }
    expected = &(I[0][0]);
    double * in = &(BiB[0][0]);
    printf_mat("B^-1*B = I? ...", 3, &(BiB[0][0]));
    assert_mat_equal( 3, in, expected, 5 * DBL_EPSILON);
    printf("checks I^-1 = I; A^-1 = Ai for diag A; B^-1 * B = I\n");
    double tt[3][3] = {{0}};
    assert_mat_equal( 3, mat_wrap(3, I, tmp, inv), simple_inv3(3, I, tt), 0);
    assert_mat_equal( 3, mat_wrap(3, A, tmp, inv), simple_inv3(3, A, tt), 0);
    assert_mat_equal( 3, mat_wrap(3, Ai, tmp, inv), simple_inv3(3, Ai, tt), 0);
    assert_mat_equal( 3, mat_wrap(3, B, tmp, inv), simple_inv3(3, B, tt), DBL_EPSILON);
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
    printf_mat("I^-1", 4, mat_wrap(4, I, tmp, inv));
    printf_mat("A^-1", 4, mat_wrap(4, A, tmp, inv));
    printf_mat("Ai", 4, &(Ai[0][0]));
    printf_mat("B^-1", 4, mat_wrap(4, B, tmp, inv));
    expected = &(I[0][0]);
    assert_mat_equal( 4, mat_wrap(4, I, tmp, inv), expected, 0.0);
    expected = &(Ai[0][0]);
    assert_mat_equal( 4, mat_wrap(4, A, tmp, inv), expected, 0.0);
    double Bi[4][4] = {{0}};
    double BiB[4][4] = {{0}};
    mat_wrap(4, B, Bi, inv);
    int i, j, k;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {
                BiB[i][j] += Bi[i][k] * B[k][j];
            }
        }
    }
    expected = &(I[0][0]);
    double * in = &(BiB[0][0]);
    printf_mat("B^-1*B = I? ...", 3, &(BiB[0][0]));
    assert_mat_equal( 4, in, expected, 5 * DBL_EPSILON);
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

int main(void)
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_det2),
        cmocka_unit_test(test_det3),
        cmocka_unit_test(test_det4),
        cmocka_unit_test(test_inv2),
        cmocka_unit_test(test_inv3),
        cmocka_unit_test(test_inv4),
        cmocka_unit_test(test_bad_inv),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
