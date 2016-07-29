/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#include <stdio.h>     /* for printf */
#include <stdlib.h>    /* for exit */
#include <getopt.h>
#include "argv.h"
#include "config.h"

int parse_argv(int argc, char ** argv)
{
    int c;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = {
            {"help",    no_argument, 0, 'h'},
            {"version", no_argument, 0, 'V'},
            {"unhandled-arg", no_argument, 0, 'x'},
            {0,         0,           0,  0 }
        };
        c = getopt_long(argc, argv, "?hV", long_options, &option_index);
        if (argc == 1) {
            c = 'h';    /* if no args, give help */
        }
        if (c == -1) {
            break;
        }
        switch (c) {
        case 'V':
            printf("%s\n", PACKAGE_STRING);
            exit(EXIT_SUCCESS);
        case '?':
        case 'h':
            printf("%s [options]\n", PACKAGE_NAME);
            printf(" --help -h     this help\n");
            printf(" --version -V  version info\n");
            if (argc == 1) {
                exit(EXIT_FAILURE);
            }
            else {
                exit(EXIT_SUCCESS);
            }
        default:
            fprintf(stderr, "%s error: ??? oops, unhandled arg\n", PACKAGE_NAME);
            exit(EXIT_FAILURE);
        }
    }
    if (optind  != argc - 1) {
        fprintf(stderr, "%s error:", PACKAGE_NAME);
        do {
            fprintf(stderr, " %s", argv[optind]);
        }
        while(++optind < argc);
        fprintf(stderr, ": extra options\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "%s: TODO parse a config file\n", PACKAGE_NAME);
    exit(EXIT_FAILURE);
    return 0;
}
