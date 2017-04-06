/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include <stdio.h>
#include "config.h"
#include "argv.h"

int main(int argc, char ** argv)
{
    args_t args = {0};
    int ret = parse_argv(argc, argv, &args);
    if(ret != 0) {
        goto quit;
    }
    if(args.mode == FORWARD_SOLVER) {
        printf("fwd %s %s\n", args.file[0], args.file[1]);
    }
quit:
    return ret;
}
