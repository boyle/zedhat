/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include <stdio.h>     /* for printf */
#include <getopt.h>
#include "config.h"
#include "argv.h"

int parse_argv(int argc, char ** argv, args_t * args)
{
    int i;
    int files = 0;
    int ret = 0;
    while (1) {
        static struct option long_options[] = {
            {"help",    no_argument, 0, 'h'},
            {"version", no_argument, 0, 'V'},
            {"fwd", no_argument, 0, 'f'},
            {"forward-solver", no_argument, 0, 'f'},
            {"unhandled-arg", no_argument, 0, 'x'},
            {0,         0,           0,  0 }
        };
        int c = getopt_long(argc, argv, "?hV", long_options, NULL);
        if (argc == 1 || argc - optind < files) {
            c = '!';    /* if no args, give help */
        }
        if (c == -1) {
            break;    /* getopt is done */
        }
        switch (c) {
        case 'V':
            printf("%s\n", PACKAGE_STRING);
            return 0;
        case '!':
            ret = 1;
        case '?':
        case 'h':
            printf("%s [options]\n", PACKAGE_NAME);
            printf(" --help -h     this help\n");
            printf(" --version -V  version info\n");
            if (argc == 1) {
                return 1;
            }
            else {
                return ret;
            }
        case 'f':
            files = 2;
            args->mode = FORWARD_SOLVER;
            break;
        default:
            fprintf(stderr, "%s error: ??? oops, unhandled arg\n", PACKAGE_NAME);
            return 1;
        }
    }
    if (argc - optind  > files) { // too many args
        fprintf(stderr, "%s error:", PACKAGE_NAME);
        do {
            fprintf(stderr, " %s", argv[optind]);
        }
        while(++optind < argc);
        fprintf(stderr, ": extra options\n");
        return 1;
    }
    for(i = 0; i < (argc - optind) && i < files; i++) {
        args->file[i] = argv[optind + i];
    }
    return 0;
}
