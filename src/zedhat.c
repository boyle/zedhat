/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include <stdio.h>
#include "config.h"
#include "argv.h"

int main(int argc, char ** argv)
{
    parse_argv(argc, argv);
    printf("test %s %s\n", "hello", "world");
    return 0;
}
