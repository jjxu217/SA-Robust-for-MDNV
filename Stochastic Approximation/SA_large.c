//
//  main.c
//  Stochastic Approximation
//
//  Created by jjxu on 2019/3/31.
//  Copyright © 2019 jjxu. All rights reserved.
//


#include <stdlib.h>
#include <ilcplex/cplex.h>
#include <math.h>
#include "support.h"

//#define DUAL_check
#define NUMofRESOURCE 20
#define NUMofPRODUCTION 30
#define NUMofROW (NUMofRESOURCE + NUMofPRODUCTION)
#define N1 1000
#define N2 1
#define Tolerance 0.00001
#define seed 17
#define alpha 2
#define beta 20

int main(int argc, char * argv[]) {
    double lambda[NUMofRESOURCE], x[NUMofRESOURCE];
    //double d[NUMofPRODUCTION];
    double step_size, improve, a = 0.01;
    int n,i;
    srand(seed);
    
    //Initialization for x
    for(i = 0; i < NUMofRESOURCE; i++){
        x[i] = 0.2;
    }
    
    for(n = 1; n < N1; n++){
        
        /*Inverse Transform Sampling
        *F(x) = 1 - (1+x^alpha)^(-beta)
        *F^-1(x) = ((1-x)^(-1/beta)-1)^(1/alpha)*/
        
        printf("\n Iteration %d, \n  ", n);
        
        solveSub(x, lambda);
        step_size = a / n ;
        improve = 0;
        for(i = 0; i < NUMofRESOURCE; i++){
            x[i] = x[i] - step_size * lambda[i];
            if(x[i] < 0){
                x[i] = 0;
            }
            improve += step_size * step_size * lambda[i] * lambda[i];
        }
        if(improve < Tolerance * Tolerance){
            break;
        }
    }
    printf("Optimal production vector: ");
    for(i = 0; i < NUMofRESOURCE; i++){
        printf("%f, ", x[i]);
    }
}


/*
 *   minimize  -c*x
 *   subject to Ay <= x
 *              y <= d
 *              l <= y <= u
 *   where

 *
 *
 */
int solveSub(double* x, double* lambda){
    int status,i, j;
    int indices[NUMofROW];
    double rhs[NUMofROW], y[NUMofPRODUCTION], one_lambda[NUMofRESOURCE];
    double d[NUMofPRODUCTION];
    double   objval, mean_objval;
    int      solnstat, solnmethod, solntype;
    int cur_numcols, cur_numrows;
    CPXENVptr     env = NULL;
    CPXLPptr      lp = NULL;
    
    for(i = 0; i < NUMofRESOURCE; i++){
        lambda[i] = 0;
    }
    mean_objval = 0.0;
    
    env = CPXopenCPLEX (&status);
    if ( env == NULL ) {
        char  errmsg[CPXMESSAGEBUFSIZE];
        fprintf (stderr, "Could not open CPLEX environment.\n");
        CPXgeterrorstring (env, status, errmsg);
        fprintf (stderr, "%s", errmsg);
        goto TERMINATE;
    }
    
    status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_OFF);
    if ( status ) {
        fprintf (stderr,
                 "Failure to turn on screen indicator, error %d.\n", status);
        goto TERMINATE;
    }
    
    lp = CPXcreateprob (env, &status, "subprob");
    if ( lp == NULL ) {
        fprintf (stderr,"Failed to create subproblem\n");
        status = 1;
        goto TERMINATE;
    }
    
    status = CPXreadcopyprob (env, lp, "subprob.lp", NULL);
    if ( status ) {
        fprintf (stderr, "Failed to read and copy the problem data.\n");
        goto TERMINATE;
    }
    
    status = CPXsetintparam (env, CPXPARAM_LPMethod, CPX_ALG_AUTOMATIC);
    if ( status ) {
        fprintf (stderr,
                 "Failed to set the optimization method, error %d.\n", status);
        goto TERMINATE;
    }
    
    for(i = 0; i < NUMofROW; i++){
        indices[i] = i;
    }
    for(i = 0; i < NUMofRESOURCE; i++){
        rhs[i] = x[i];
    }
    
    for(j = 0; j < N2; j++){
        for(i = 0; i < NUMofPRODUCTION; i++){
            d[i] = (rand() % 1000) / 1000.0;
            d[i] = pow(pow(1 - d[i], -1.0/beta) - 1.0, 1.0/alpha);
        }
    
        for(i = 0; i < NUMofPRODUCTION; i++){
            rhs[NUMofRESOURCE + i] = d[i];
        }
    
        status = CPXchgrhs (env, lp, NUMofROW, indices, rhs);
        if ( status ) {
            fprintf (stderr, "CPXchgrhs failed, error code %d.\n", status);
            goto TERMINATE;
        }

        status = CPXlpopt (env, lp);
        if ( status ) {
            fprintf (stderr, "Failed to optimize LP.\n");
            goto TERMINATE;
        }
    
        solnstat = CPXgetstat (env, lp);
    
        if      ( solnstat == CPX_STAT_UNBOUNDED ) {
            printf ("Model is unbounded\n");
            goto TERMINATE;
        }
        else if ( solnstat == CPX_STAT_INFEASIBLE ) {
            printf ("Model is infeasible\n");
            goto TERMINATE;
        }
        else if ( solnstat == CPX_STAT_INForUNBD ) {
            printf ("Model is infeasible or unbounded\n");
            goto TERMINATE;
        }
    
        status = CPXsolninfo (env, lp, &solnmethod, &solntype, NULL, NULL);
        if ( status ) {
            fprintf (stderr, "Failed to obtain solution info.\n");
            goto TERMINATE;
        }
        //printf ("Solution status %d, solution method %d\n", solnstat, solnmethod);
    
        if ( solntype == CPX_NO_SOLN ) {
            fprintf (stderr, "Solution not available.\n");
            goto TERMINATE;
        }
    
        status = CPXgetobjval (env, lp, &objval);
        mean_objval += objval;
        if ( status ) {
            fprintf (stderr, "Failed to obtain objective value.\n");
            goto TERMINATE;
        }
        
    
        status = CPXgetx (env, lp, y, 0, NUMofPRODUCTION-1);
        if ( status ) {
            fprintf (stderr, "Failed to get y.\n");
            goto TERMINATE;
        }
        cur_numcols = CPXgetnumcols (env, lp);
        cur_numrows = CPXgetnumrows (env, lp);

        status = CPXgetpi(env, lp, one_lambda, 0, NUMofRESOURCE-1);
        
        for(i = 0; i < NUMofRESOURCE; i++){
            lambda[i] += one_lambda[i];
        }
        
        if ( status ) {
            fprintf (stderr, "Failed to get pi.\n");
            goto TERMINATE;
        }
    
    #ifdef DUAL_check
        status = CPXgetpi(env, lp, pi, 0, NUMofROW-1);
        printf("Dual variables: ");
        for(i=0;i<NUMofROW;i++){
            printf("%f ", pi[i]);
        }
    #endif
    }
    
    for(i = 0; i < NUMofRESOURCE; i++){
        lambda[i] = lambda[i] / N2;
    }
    mean_objval = mean_objval / N2;
    printf ("Objective value %.10g.\n", mean_objval);
    
    
TERMINATE:
    /* Free up the problem as allocated by CPXcreateprob, if necessary */
    if ( lp != NULL ) {
        status = CPXfreeprob (env, &lp);
        if ( status ) {
            fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
        }
    }
    
    /* Free up the CPLEX environment, if necessary */
    
    if ( env != NULL ) {
        status = CPXcloseCPLEX (&env);
        
        /* Note that CPXcloseCPLEX produces no output,
         so the only way to see the cause of the error is to use
         CPXgeterrorstring.  For other CPLEX routines, the errors will
         be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON. */
        
        if ( status ) {
            char  errmsg[CPXMESSAGEBUFSIZE];
            fprintf (stderr, "Could not close CPLEX environment.\n");
            CPXgeterrorstring (env, status, errmsg);
            fprintf (stderr, "%s", errmsg);
        }
    }
    
    return 0;
}
