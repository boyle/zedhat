/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdio.h>
#include <stdlib.h> /* malloc, free */
#include <string.h> /* memcpy, memset */
#include <math.h> /* sqrt */
#include <assert.h> /* assert */
#include <time.h> /* clock */
//#include <cblas.h>
//#include <lapacke.h>

#include "argv.h"
#include "file.h"
#include "fwd.h"
#include "inv.h"
#include "model.h"
#include "matrix.h"

static int do_fwd_solve(const char * file, double tol)
{
    int ret = 1;
    model  * mdl = malloc_model();
    printf("-- loading from file %s --\n", file);
    int nret = readfile(file, mdl);
    const int n_meas = mdl->n_stimmeas;
    const int frames = mdl->n_params[1];
    double meas[frames * n_meas];
    if(!nret) {
        fprintf(stderr, "error: read %s failed\n", file);
        goto fwd_quit;
    }
    printf("-- forward solutions --\n");
    printf("%d parameter frame%s\n", frames, frames == 1 ? "" : "s");
    if(frames < 1) {
        fprintf(stderr, "error: nothing to do\n");
        goto fwd_quit;
    }
    if(mdl->n_data[0] > 0) {
        assert(mdl->n_data[0] == mdl->n_stimmeas);
        assert(mdl->n_data[1] == mdl->n_params[1]);
    }
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
                const double data_expect = mdl->data[i + mdl->n_data[0] * idx];
                printf("%4d  %18.8g %18.8g\n", i + 1, data_calc, data_expect);
                const double dm = data_calc - data_expect;
                rmse += dm * dm;
            }
            rmse = sqrt(rmse / (double) mdl->n_data[0]);
            printf("                      RMSE = %g\n", rmse);
            printf("                      tolerance = %g\n", tol);
            if (rmse > tol) {
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
    ret = 0;
    printf("-- completed --\n");
fwd_quit:
    free_model(mdl);
    return ret;
}

#define max(a,b) (((a) > (b)) ? (a) : (b))
static int do_inv_solve(const char * file, const double tol)
{
    int ret = 1;
    model  * mdl = malloc_model();
    printf("-- loading from file %s --\n", file);
    int nret = readfile(file, mdl);
    const int frames = mdl->n_data[1] > 1 ? mdl->n_data[1] - 1 : 0;
    if(!nret) {
        fprintf(stderr, "error: read %s failed\n", file);
        free_model(mdl);
        return 1;
    }
    printf("-- inverse solutions --\n");
    printf("%d measurement frame%s\n", frames, frames == 1 ? "" : "s");
    if(frames < 1) {
        fprintf(stderr, "error: nothing to do\n");
        free_model(mdl);
        return 1;
    }
    printf_model(mdl, 1);
    assert(mdl->n_data[0] > 0);
    assert(mdl->n_data[0] == mdl->n_stimmeas);
    if(mdl->n_params[1] > 1) {
        assert(mdl->n_data[1] == mdl->n_params[1]);
    }
    const int is_pmap = (mdl->fwd.n_pmap > 0);
    int n_params = 0;
    if( mdl->n_params[1] > 0 ) {
        n_params = mdl->n_params[0];
    }
    else if( !is_pmap ) {
        n_params = mdl->fwd.n_elems;
    }
    else {
        for(int i = 0; i < mdl->fwd.n_pmap; i++) {
            n_params = max(n_params, mdl->fwd.pmap_param[i]);
        }
        n_params++;
    }
    mdl->n_params[0] = n_params;
    double params[frames * n_params];
    memset(params, 0, sizeof(double) * frames * n_params);
    nret = inv_solve(mdl, params);
    if(!nret) {
        fprintf(stderr, "error: bad inverse solve\n");
        goto inv_quit;
    }
    for(int idx = 0; idx < frames; idx++) {
        /* calculated measurements to expected */
        if(mdl->n_params[1] > 1) {
            printf("frame#%d\n", idx + 1);
            printf("param#%18s %18s\n", "calculated", "from file");
            const int cols = mdl->n_params[1];
            const int rows = mdl->n_params[0];
            assert(cols == 3); /* TODO currently: params = [ J_background, ignore, expected ] */
            double rmse = 0; /* root mean squared error from expected */
            for(int i = 0; i < rows; i++) {
                const double params_calc = params[i + 0 * rows];
                const double params_expect = mdl->params[i + 2 * rows];
                printf("%4d  %18.8e %18.8e\n", i + 1, params_calc, params_expect);
                const double dp = params_calc - params_expect;
                rmse += dp * dp;
            }
            rmse = sqrt(rmse / (double) n_params);
            printf("                      RMSE = %g\n", rmse);
            if (rmse > tol) {
                fprintf(stderr, "fail: calculated parameters do not match expected\n");
                goto inv_quit;
            }
            break; /* TODO only handle one frame in checking mode currently */
        }
        else {
            printf("frame#%d\n", idx + 1);
            printf("param#%18s\n", "calculated");
            const int rows = mdl->n_params[0];
            for(int i = 0; i < n_params; i++) {
                const double params_calc = params[i + idx * rows];
                printf("%4d  %18.04g\n", i + 1, params_calc);
            }
        }
    }
    ret = 0;
    printf("-- completed --\n");
inv_quit:
    free_model(mdl);
    return ret;
}

int main(int argc, char ** argv)
{
    int ret = 0;
    args_t args;
    const int N = 1;
    long double dt [N];
    long double mean_dt = 0;
    long double std_dt = 0;
    if(parse_argv(argc, argv, &args) != 0) {
        return 1;
    }
    for(int i = 0; i < N; i++) {
        clock_t start = clock();
        printf("%s %s %d\n", args.file[0], args.file[1], args.mode);
        switch(args.mode) {
        case FORWARD_SOLVER: ret = do_fwd_solve(args.file[0], args.tol); break;
        case INVERSE_SOLVER: ret = do_inv_solve(args.file[1], args.tol); break;
        default: return 0;
        }
        clock_t stop = clock();
        dt[i] = (stop - start) / CLOCKS_PER_SEC;
        mean_dt += dt[i] / N;
        if(ret) {
            break;
        }
    }
    if(ret) {
        printf("error: solve failed\n");
        return ret;
    }
    for(int i = 0; i < N; i++) {
        const double long tmp = dt[i] - mean_dt;
        std_dt += tmp * tmp;
    }
    std_dt = sqrt(std_dt / N);
    printf("%Lf\n", mean_dt);
    printf("runtime: %0.0f:%02.0f:%02.4Lfs ± %2.4Lfs\n", 0.0, 0.0, mean_dt, std_dt);
    printf("total:   %0.0f:%02.0f:%02.4Lfs\n", 0.0, 0.0, mean_dt*N);
    return 0;
}
