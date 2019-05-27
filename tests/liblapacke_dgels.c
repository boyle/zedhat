/* from lapack-netlib, 3-clause BSD licensed */
/* modified (C) 2016, A. Boyle */

#include "config.h"
#include <stdio.h>
#include <math.h> /* for fabs */
#include <lapacke.h>

int main (int argc, const char * argv[])
{
    double A[] = {1, 1, 1,
                  2, 3, 4,
                  3, 5, 2,
                  4, 2, 5,
                  5, 4, 3
                 };
    double b[] = { -10, -3,
                   12, 14,
                   14, 12,
                   16, 16,
                   18, 16
                 };
    int i, j;
    /* Solve least squares problem: DGELS using row-major layout */
    int info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', 5, 3, 2, A, 3, b, 2);
    if (info) {
        printf("info=%d\n", info);
        exit(EXIT_FAILURE);
    }
    /* matlab: expect = A \ b */
    double const expect[] = {2, 1,
                             1, 1,
                             1, 2
                            };
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 2; j++) {
            printf(" %g", b[i * 2 + j]);
        }
        printf("\n");
    }
    for(i = 0; i < 6; i++) {
        if(fabs(expect[i] - b[i]) > 4e-15) {
            printf("%d %g\n", i, fabs(expect[i] - b[i]));
            printf("FAIL\n");
            exit(EXIT_FAILURE);
        }
    }
    printf("PASS\n");
    exit(EXIT_SUCCESS);
}
