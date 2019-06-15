/* Copyright 2019, Alistair Boyle, 3-clause BSD License */
#ifndef __FWD_H__
#define __FWD_H__

#include "model.h"
#include "matrix.h"

int fwd_solve(model * mdl, double * meas);
int calc_jacobian(model * mdl, double * J);

#endif
