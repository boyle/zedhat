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
      printf("|a:%g - b:%g| = %g > %g\n", a, b, fabs(a-b), delta); \
   } \
   assert_true(fabs(a - b) <= delta); \
} while(0)

void test_det3( void ** state)
{
    double I[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
    double A[3][3] = {
        {1, 9, 3},
        {4, 4, 6},
        {7, 8, 9}
    };
    double At[3][3] = {
        {1, 4, 7},
        {9, 4, 8},
        {3, 6, 9}
    };
    double B[3][3] = {
        {0.1, 5.3, 3.3},
        {1.1, 2.1, 7.2},
        {4.1, 0.8, 9.1}
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
}

int main(void)
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_det3),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
