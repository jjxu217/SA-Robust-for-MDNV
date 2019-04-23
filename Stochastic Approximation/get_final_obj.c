//
//  get_final_obj.c
//  Stochastic Approximation
//
//  Created by jjxu on 2019/4/1.
//  Copyright © 2019 jjxu. All rights reserved.
//



#include <stdlib.h>
#include <ilcplex/cplex.h>
#include <math.h>



#define NUMofRESOURCE 2
#define NUMofPRODUCTION 3
#define NUMofROW (NUMofRESOURCE + NUMofPRODUCTION)
#define alpha 2
#define beta 20

//int main(int argc, char * argv[]) {
//    double x[NUMofRESOURCE] = {0.5483855, 1.143005, 0.2710045,  0.6335799,  0.2652503, 6.307027e-01, 2.594962e-01  , 6.278255e-01  , 2.537421e-01, 6.249483e-01, 2.479880e-01   ,6.220712e-01  ,2.422339e-01  , 6.191940e-01 ,2.364798e-01  ,6.163168e-01  ,2.307258e-01  ,6.134396e-01  , 2.249718e-01, 6.105624e-01 };
//    get_final_obj(x, 100);
//}

double get_final_obj(double* x, int runlength){
    int status,i, j;
    int indices[NUMofROW];
    double rhs[NUMofROW], y[NUMofPRODUCTION], one_lambda[NUMofRESOURCE];
    double d[NUMofPRODUCTION];
    double   objval, mean_objval;
    int      solnstat, solnmethod, solntype;
    int cur_numcols, cur_numrows;
    CPXENVptr     env = NULL;
    CPXLPptr      lp = NULL;
    
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
    
    for(j = 0; j < runlength; j++){
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
        CPXwriteprob(env, lp, "test.lp", NULL);
        
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
        
        if ( status ) {
            fprintf (stderr, "Failed to get pi.\n");
            goto TERMINATE;
        }
        
    }
    
    mean_objval = 5.0 * x[0] + 5.0 * x[1] + mean_objval / runlength;
    printf ("Final Objective value %.10g.\n", mean_objval);
    return mean_objval;
    
    
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
