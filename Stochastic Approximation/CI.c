//
//  CI.c
//  Stochastic Approximation
//
//  Created by jjxu on 2019/4/10.
//  Copyright © 2019 jjxu. All rights reserved.
//

#include "CI.h"
#include <stdlib.h>
#include <stdio.h>
#include <ilcplex/cplex.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "support.h"
#define NUMofRESOURCE 2
#define NUMofPRODUCTION 3
#define NUMofROW (NUMofRESOURCE + NUMofPRODUCTION)
#define N1 10000
#define N2 10
#define Tolerance 0.0001

#define alpha 2
#define beta 20

#define NUMofREPLICATION 10


int main(int argc, char * argv[]) {
    double x[2] = { 0.1734,0.2296 };
    double CI[2];
    double temp;
    temp = get_CI(CI, x);
    printf("Objective functions Expectation %f \n", temp);
    printf("Objective functions 0.95 CI: [%.4f,%.4f] \n", CI[0], CI[1]);
}

double get_CI(double *CI, double *x){
    int status,i, count;
    int indices[NUMofROW];
    double rhs[NUMofROW];
    double d[NUMofPRODUCTION];
    double   objval, mean=0, stdev,temp=0,vari,ans;
    long long seed=2668655841019641;
    CPXENVptr     env = NULL;
    CPXLPptr      lp = NULL;
    
    
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
    
    lp = CPXcreateprob (env, &status, "subprob_small");
    if ( lp == NULL ) {
        fprintf (stderr,"Failed to create subproblem\n");
        status = 1;
        goto TERMINATE;
    }
    
    status = CPXreadcopyprob (env, lp, "subprob_small.lp", NULL);
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
    
    count = 0;
    mean = 0.0;
    vari = 0.0;
    stdev = 10000000.0;
    temp = 0.0;
    while (3.92 * stdev > 0.01 * fabs(mean) || count < 100){
        for(i = 0; i < NUMofPRODUCTION; i++){
            d[i] = scalit(0.0, 1.0, &seed);
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
        
        objval = 0.0;
        status = CPXgetobjval (env, lp, &objval);
        
        ans =  5.0 * x[0] + 5.0 * x[1] + objval;
        
        if ( status ) {
            fprintf (stderr, "Failed to obtain objective value.\n");
            goto TERMINATE;
        }
        
        
        if (count == 0)
        {
            mean = ans;
        }
        else
        {
            temp = mean;
            mean = mean + (ans - mean) / (double) (count + 1);
            vari = (1 - 1 / (double) count) * vari
            + (count + 1) * (mean - temp) * (mean - temp);
            stdev = sqrt(vari / (double) count);
        }
        
        count++;
        
        /* Print the results every once in a while for long runs */
        if (!(count % 250))
        {
            printf(".");
        }
        if (!(count % 10000))
        {
            printf("\n\nobs:%d mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n",
                   count, mean, 3.92 * stdev / mean,
                   mean - 1.96 * stdev, mean + 1.96 * stdev);
        }
    }
    
    CI[0] = mean - 1.96*stdev;
    CI[1] = mean + 1.96*stdev;
    
    return mean;
    
        
        
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
        return 0.0;
}


float scalit(float lower, float upper, long long *RUN_SEED)
{
    float val, wide;
    
    wide = upper - lower;
    
    val = randUniform(RUN_SEED);
    
    /* Yifan 06/25/2012 batch mean */
    //val = randUniform(&config.RUN_SEED1);
    return ((wide * val) + lower);
}

float randUniform(long long *SEED)
{
    /* static int to static long int: modified by Yifan 2013.02.18 */
    static int lo_bits, hi_bits;
    
    lo_bits = ((*SEED) & 0xFFFFL) * 16807;
    hi_bits = (int) (((*SEED) >> 16) * 16807) + (lo_bits >> 16);
    *SEED = ((lo_bits & 0xFFFFL) - 0x7FFFFFFFL) + ((hi_bits & 0x7FFFL) << 16)
    + (hi_bits >> 15);
    return ((*SEED) < 0 ? ((*SEED) += 0x7FFFFFFFL) : (*SEED)) * 4.656612875E-10;
}
