/* Copyright 2016, 2018, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdio.h> /* for printf */
#include <stdlib.h> /* strtof */
#include <errno.h> /* strtof -> errno */
#include <getopt.h>
#include <assert.h> /* assert */
#include <math.h> /* INFINITY */
#include <string.h> /* bzero */
#include "argv.h"

int parse_argv(int argc, char ** argv, args_t * args)
{
    int err = (argc <= 1);
    assert(args != NULL);
    bzero(args, sizeof(args_t));
    args->tol = INFINITY;
    while (1) {
        static struct option long_options[] = {
            {"help",    no_argument, 0, 'h'},
            {"version", no_argument, 0, 'V'},
            {"fwd", required_argument, 0, 'f'},
            {"forward-solver", required_argument, 0, 'f'},
            {"inv", required_argument, 0, 'i'},
            {"inverse-solver", required_argument, 0, 'i'},
            {"tolerance", required_argument, 0, 't'},
            {"unhandled-arg", no_argument, 0, 'x'},
            {0,         0,           0,  0 }
        };
        int c = getopt_long(argc, argv, "?hVf:i:t:", long_options, NULL);
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
        case 'i':
            args->mode = INVERSE_SOLVER;
            args->file[1] = optarg;
            break;
        case 't':
            errno = 0;
            char * endptr;
            args->tol = strtod(optarg, &endptr);
            if((errno != 0) || (optarg == endptr)) {
                fprintf(stderr, "error: --tolerance %s: not a number\n", optarg);
                err = 1;
            }
            if(args->tol < 0) {
                fprintf(stderr, "error: --tolerance %s: must be non-negative\n", optarg);
                err = 1;
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
    printf(" --forward-solver --fwd -f <fwd.zh>\n");
    printf("       simulate measurements for model in fwd.zh\n");
    printf(" --inverse-solver --inv -f <inv.zh>\n");
    printf("       solve for model parameters using model and data from inv.zh\n");
    printf(" --tolerance --tol -t <#.##e#>\n");
    printf("       for checking --fwd and --inv solutions when available\n");
    return err;
}
