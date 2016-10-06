/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include <stdio.h>
#include "config.h"
#include "argv.h"

int main(int argc, char ** argv)
{
    int ret = parse_argv(argc, argv);
//    if(ret != 0) {
    ret = (ret == 1);
    goto quit;
//    }
//    printf("test %s %s\n", "hello", "world");
quit:
    return ret;
}
