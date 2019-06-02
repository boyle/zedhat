/* Copyright 2018, Alistair Boyle, 3-clause BSD License */
#ifndef __FILE_H__
#define __FILE_H__

#include "model.h"
int readfile_ngvol(const char filename[], model * m);
int readfile(const char filename[], model * m);

#endif
