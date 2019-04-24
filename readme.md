README

This project contains the Stochatic Approximation (SA) algorithm and Robust Stochatic Approximation (RSA) algorithm 
for Multi-dimensional Newsvendor (MDNV) problem. The code is written in C language. 

To run the SA, please select SA_small.c and support.h files for compile.
To run the RSA, please select Robust_SA_small.c and support.h files for compile.

subprob_small.lp is the model file for the subproblem, which is needed for both algorithm. 

CI.c and CI.h calculate the confidence interval of the objective value, given the first stage variables.
