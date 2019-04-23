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
#include <time.h>
#include <string.h>
#include "support.h"
//#include "get_final_obj.c"

//#define DUAL_check 1
#define NUMofRESOURCE 2
#define NUMofPRODUCTION 3
#define NUMofROW (NUMofRESOURCE + NUMofPRODUCTION)
#define N1 10000
#define N2 10

#define alpha 2
#define beta 20

#define NUMofREPLICATION 10





int main(int argc, char * argv[]) {
    double lambda[NUMofRESOURCE], x[NUMofRESOURCE], x_mean[NUMofRESOURCE];
    //double d[NUMofPRODUCTION];
    double step_size, a = 0.0001;
    int n,i,r;
    double mean_obj=0.0, duration;
    long long seed[30] = {4650175399072632,3554548844580680,6070772756632709,5451675876709589,5285327724846206,5588857889468088,1098833779416153,6192593982049265,4756774140130874,6784592265109609,9728429908537680,1163479388309571,3279282318700126,8773753208032360,9337302665697748,4415169667296773,4220432037464045,3554548844580680,1814300451929103,5339672949292608,5638710736762732,3154245808720589,2414929536171258,7998609999427572,7080145164625719,3612848862740490586,7772725003305823,5982768791029230,1395182510837913,3735836402047426};
    FILE *oFile;
    clock_t start, finish;
    double CI[2];
    
    oFile = open_ofile(argv[1]);
    fprintf(oFile, "RobustSA Algo: \n");
    fclose(oFile);
    
    
    for(r = 0; r < NUMofREPLICATION;r++){
        srand((unsigned)seed[r]);
        
        //Initialization for x
        for(i = 0; i < NUMofRESOURCE; i++){
            x[i] = 0.1;
        }
        
        start = clock();
        for(n = 1; n < N1; n++){
            
            /*Inverse Transform Sampling
            *F(x) = 1 - (1+x^alpha)^(-beta)
            *F^-1(x) = ((1-x)^(-1/beta)-1)^(1/alpha)*/
            
            printf("\n Iteration %d, \n  ", n);
            
            solveSub(x, lambda, &seed[r]);
            step_size = a ;
    //        improve = 0;
            for(i = 0; i < NUMofRESOURCE; i++){
                x[i] = x[i] - step_size * (5 + lambda[i]);
                if(x[i] < 0){
                    x[i] = 0;
                }
                if(x[i] > 5){
                    x[i] = 5;
                }
                x_mean[i] += x[i];
    //            improve += step_size * step_size * lambda[i] * lambda[i];
            }
    //        if(improve < Tolerance * Tolerance){
    //            break;
    //        }
        }
        finish = clock();
        duration = (double)(finish - start)/ CLOCKS_PER_SEC;
        
        printf("Optimal Solution: ");
        for(i = 0; i < NUMofRESOURCE; i++){
            x_mean[i] = x_mean[i]/N1;
            printf("%f ", x_mean[i]);
        }
        
        mean_obj = get_CI(CI,x);
        
        
        
        /********Output file ***************/
        oFile = open_ofile(argv[1]);
        
        fprintf(oFile, "%d-th Repelication: \n", r);
        fprintf(oFile, "Solution, 0.95 CI, time : &[%.4f,%.4f] &[%.4f,%.4f] &%.2f \n", x[0], x[1], CI[0], CI[1], duration);
        fprintf(oFile, "Objective functions = %f\n", mean_obj);
        
        fprintf(oFile, "NUMofRESOURCE, NUMofPRODUCTION, NUMofROW = %d, %d, %d\n", NUMofRESOURCE, NUMofPRODUCTION, NUMofROW);
        fprintf(oFile, "N1, N2 = %d, %d\n", N1, N2);
        fprintf(oFile, "seed =  %lld\n", seed[r]);
        fprintf(oFile, "constant step size =  %f\n\n", a);
        
        
        fclose(oFile);
    }

    
    /********Output file end***************/
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
int solveSub(double* x, double* lambda, long long *SEED){
    int status,i, j;
    int indices[NUMofROW];
    double rhs[NUMofROW], y[NUMofPRODUCTION], one_lambda[NUMofRESOURCE], pi[NUMofROW];
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
    
    for(j = 0; j < N2; j++){
        for(i = 0; i < NUMofPRODUCTION; i++){
            d[i] = scalit(0.0, 1.0, SEED);
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
        mean_objval = mean_objval + objval + 5*x[0]+5*x[1];
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
//        printf("Dual variables: ");
//        for(i=0;i<NUMofROW;i++){
//            printf("%f ", pi[i]);
//        }
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

FILE *open_ofile(char* name){
    FILE *fp = NULL;
    char probName[100];
    
    strcpy(probName, name);
    strcat(probName, ".out");
    fp = fopen(probName, "a+");
    if ( fp == NULL ) {
        fprintf(stderr, "failed to open file %s", probName);
        return NULL;
    }
    return fp;
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
