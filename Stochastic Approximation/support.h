//
//  supporting.h
//  Stochastic Approximation
//
//  Created by jjxu on 2019/3/31.
//  Copyright © 2019 jjxu. All rights reserved.
//



#include <stdio.h>
typedef        double                *vector;

int solveSub(vector x,vector lambda, long long *SEED);
FILE *open_ofile(char *name);

double get_CI(double *CI, double* x);
float scalit(float lower, float upper, long long *RUN_SEED);
float randUniform(long long *SEED);
