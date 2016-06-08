/* from lapack-netlib, 3-clause BSD licensed */
/* modified 2016, A. Boyle */

#include <stdio.h>
#include <lapacke.h>

int main (int argc, const char * argv[])
{
    double A[] = {1, 1, 1, 2, 3, 4, 3, 5, 2, 4, 2, 5, 5, 4, 3}; /* 5 x 3 */
    double b[] = { -10, -3, 12, 14, 14, 12, 16, 16, 18, 16}; /* 5 x 2 */
    int i;
    /* Solve least squares problem: DGELS using row-major layout */
    int info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', 5, 3, 2, A, 3, b, 2);
    if (info) {
        printf("info=%d\n", info);
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < 3; i++) {
        printf(" %g", b[i]);
    } printf( "\n" );
    for(i = 3; i < 6; i++) {
        printf(" %g", b[i]);
    } printf( "\n" );
    exit( EXIT_SUCCESS );
}
