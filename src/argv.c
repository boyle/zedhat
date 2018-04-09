/* Copyright 2016, 2018, Alistair Boyle, 3-clause BSD License */
#include <stdio.h> /* for printf */
#include <stdlib.h> /* strtof */
#include <errno.h> /* strtof -> errno */
#include <getopt.h>
#include "config.h"
#include "argv.h"

int parse_argv(int argc, char ** argv, args_t * args)
{
    int err = (argc <= 1);
    while (1) {
        static struct option long_options[] = {
            {"help",    no_argument, 0, 'h'},
            {"version", no_argument, 0, 'V'},
            {"fwd", required_argument, 0, 'f'},
            {"forward-solver", required_argument, 0, 'f'},
            {"tolerance", required_argument, 0, 't'},
            {"unhandled-arg", no_argument, 0, 'x'},
            {0,         0,           0,  0 }
        };
        int c = getopt_long(argc, argv, "?hVf:t:", long_options, NULL);
        if (c == -1) {
            break;    /* getopt is done */
        }
        switch (c) {
        case 'V':
            printf("%s\n", PACKAGE_STRING);
            break;
        case '?':
            err = 1;
        case 'h':
            goto _help;
        case 'f':
            args->mode = FORWARD_SOLVER;
            args->file[0] = optarg;
            break;
        case 't':
            errno = 0;
            args->tol = strtod(optarg, NULL);
            if (errno != 0) {
                args->tol = -1e10;
            }
            break;
        default:
            fprintf(stderr, "error: %s: unhandled argument\n", PACKAGE_NAME);
            err = 1;
        }
    }
    if (optind != argc) {
        fprintf(stderr, "error: %s: extra option\n", argv[optind]);
        err = 1;
    }
    if (!err) {
        return 0;
    }
_help:
    printf("%s [options]\n", PACKAGE_NAME);
    printf(" --help -h     this help\n");
    printf(" --version -V  version info\n");
    printf(" --forward-solver --fwd -f <file.mat>\n");
    printf("               solves X=A\\B\n");
    printf("               with A B from file.mat\n");
    return err;
}
