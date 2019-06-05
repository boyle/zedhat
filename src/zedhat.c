/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdio.h>
#include <stdlib.h> /* malloc, free */
#include <string.h> /* memcpy */
#include <math.h> /* sqrt */
#include <assert.h> /* assert */
//#include <cblas.h>
//#include <lapacke.h>

#include "argv.h"
#include "file.h"
#include "fwd.h"
#include "model.h"
#include "matrix.h"

int main(int argc, char ** argv)
{
    args_t args = {0};
    int ret = parse_argv(argc, argv, &args);
    if(ret != 0) {
        return 1;
    }
    if(args.mode == FORWARD_SOLVER) {
        ret = 1;
        model * mdl = malloc_model();
        printf("-- loading from file --\n");
        int nret = readfile(args.file[0], mdl);
        const int n_meas = mdl->n_stimmeas;
        const int frames = mdl->n_params[1];
        double meas[frames * n_meas];
        if(!nret) {
            fprintf(stderr, "error: read %s failed\n", args.file[0]);
            goto fwd_quit;
        }
        if(mdl->n_data[0] > 0) {
            assert(mdl->n_data[0] == mdl->n_stimmeas);
            assert(mdl->n_data[1] == mdl->n_params[1]);
        }
        printf("-- forward solutions --\n");
        nret = fwd_solve(mdl, meas);
        if(!nret) {
            fprintf(stderr, "error: bad forward solve\n");
            goto fwd_quit;
        }
        for(int idx = 0; idx < frames; idx++) {
            /* calculated measurements to expected */
            if(mdl->n_data[0] > 0) {
                printf("frame#%d\n", idx + 1);
                printf("meas#  %18s %18s\n", "calculated", "from file");
                double rmse = 0; /* root mean squared error from expected */
                for(int i = 0; i < mdl->n_stimmeas; i++) {
                    const double data_calc = meas[idx * n_meas + i];
                    const double data_expect = mdl->data[i * mdl->n_data[1] + idx];
                    printf("%4d  %18.8g %18.8g\n", i + 1, data_calc, data_expect);
                    const double dm = data_calc - data_expect;
                    rmse += dm * dm;
                }
                rmse = sqrt(rmse / (double) mdl->n_data[0]);
                printf("                      RMSE = %g\n", rmse);
                if (rmse > args.tol) {
                    fprintf(stderr, "fail: calculated measurements do not match expected\n");
                    ret = 1;
                    goto fwd_quit;
                }
            }
            else {
                printf("frame#%d\n", idx + 1);
                printf("meas#  %18s\n", "calculated");
                for(int i = 0; i < mdl->n_stimmeas; i++) {
                    const double data_calc = meas[idx * n_meas + i];
                    printf("%4d  %18.04g\n", i + 1, data_calc);
                }
            }
        }
        // /* Find: || A X - B ||_2 */
        // /* Compute  C = a A B + c C */
        // cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
        //              A->m, X->n, A->n, /* m, n, k, */
        //              +1.0, A->x.dense, A->m, X->x.dense, X->m, /* a, A, m, B, n, */
        //              -1.0, BB, B->m); /* c, C, k */
        // /* Compute  || B ||_2 */
        // const double err_mul = cblas_dnrm2( Bmn, BB, 1); /* n, X, incX */
        // if (err_mul > args.tol) {
        //     printf("fail: || A X - B ||_2 = %g\n", err_mul);
        //     ret = 1;
        //     goto fwd_quit;
        // }
        // /* Find BB = A\B; || BB - X ||_2 */
        // memcpy(AA, A->x.dense, Amn * sizeof(double));
        // memcpy(BB, B->x.dense, Bmn * sizeof(double));
        // // memset(&BB[Bmn], 0, (Xmn-Bmn)*sizeof(double));
        // lapack_int err = LAPACKE_dgels( LAPACK_COL_MAJOR, 'N', /* row/col major, trans, */
        //                                 A->m, A->n, B->n, /* m, n, n_rhs, */
        //                                 AA, A->m, BB, B->m ); /* A, m, B, n */
        // if (err) {
        //     printf("fail: DGELS (X=A\\B) info=%d\n", err);
        //     ret = 1;
        //     goto fwd_quit;
        // }
        // /* Compute  || BB - X ||_2 */
        // for(int i = 0; i < Xmn; i++) {
        //     BB[i] -= X->x.dense[i];
        // }
        // const double err_div = cblas_dnrm2( Xmn, BB, 1); /* n, X, incX */
        // if (err_div > args.tol) {
        //     printf("fail: || A\\B - X ||_2 = %g\n", err_div);
        //     ret = 1;
        //     goto fwd_quit;
        // }
        printf("-- completed --\n");
        ret = 0;
fwd_quit:
        free_model(mdl);
    }
    return ret;
}
