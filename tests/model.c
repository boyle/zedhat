/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#include <stdarg.h> /* cmocka.h */
#include <stddef.h> /* cmocka.h */
#include <setjmp.h> /* cmocka.h */

#include <stdio.h> /* cmocka.h */
#include <string.h> /* cmocka.h */
#include <float.h> /* DBL_EPSILON */
#include <math.h> /* fabs */

#include "cmocka.h"
#include "model.h"

#define assert_double_equal(a, b, delta) do { \
   if(!(fabs(a - b) <= delta)) { \
      printf("|a:%0.16g - b:%0.16g| = %0.16g > %0.16g\n", a, b, fabs(a-b), delta); \
   } \
   assert_true(fabs(a - b) <= delta); \
} while(0)

double simple_det3(int n, double A[3][3])
{
    assert_int_equal(n,3);

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

int main(void)
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_det2),
        cmocka_unit_test(test_det3),
        cmocka_unit_test(test_det4),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
