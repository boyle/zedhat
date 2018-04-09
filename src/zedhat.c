/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include <stdio.h>
#include "config.h"
#include "argv.h"
#include "matrix.h"

matrix_t * get_fwd_matrix(const args_t args, const char * name)
{
    matrix_t * A = matrix_malloc(name);
    int ret = matrix_load(args.file[0], A);
    printf("fwd %s: ret = %d (%s)\n", name, ret, args.file[0]);
    if(ret) {
        matrix_free(A);
        return NULL;
    }
    else {
        return A;
    }
}

int main(int argc, char ** argv)
{
    args_t args = {0};
    int ret = parse_argv(argc, argv, &args);
    if(ret != 0) {
        goto quit;
    }
    if(args.mode == FORWARD_SOLVER) {
        matrix_t * A = get_fwd_matrix(args, "A");
        matrix_t * B = get_fwd_matrix(args, "B");
        matrix_t * X = get_fwd_matrix(args, "X");
        ret = (!A || !B || !X); /* all ptrs good */
        matrix_free(A);
        matrix_free(B);
        matrix_free(X);
        goto quit;
    }
quit:
    return ret;
}
