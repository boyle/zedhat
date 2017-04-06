/* Copyright 2016, Alistair Boyle, 3-clause BSD License */
#ifndef __ARGV_H__
#define __ARGV_H__

typedef enum {IDLE, FORWARD_SOLVER, INVERSE_SOLVER} args_mode_t;

typedef struct args {
    args_mode_t mode;
    const char * file[2];
} args_t;

/* returns: 0 on success, 1 on failure, 2 on success-but-exit-now */
int parse_argv(int argc, char ** argv, args_t * args);
#endif /* __ARGV_H__ */

