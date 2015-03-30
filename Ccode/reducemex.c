/*Copyright 2015 Jeff Shamma
 Address:
 School of Electrical and Computer Engineering
 Georgia Institute of Technology
 777 Atlantic Dr NW
 Atlanta, GA 30332-0250
*/
#include "cplex.h"
#include <string.h>
#include <math.h>
#include "mex.h"

/* mex-file gateway

   [costs,stats,ikeep] = reducemex(M,matbeg,matcnt,matind,matval,
                          rhs,ncon,nvar,numnz,tol)
*/

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])

{

/* Declare variables with "natural" cplex data type */
int       *matbeg; /* columnwise beginning index    */
int       *matcnt; /* columnwise number of elements */
int       *matind; /* columnwise row indices        */
double    *matval; /* values */
char      *sense;
int       objsen;
double    *lb;
double    *ub;
int       numrows;
int       numcols;
int       numnz;

double    *M,*m,*costs,*stats,*ikeep,tol;

/* Declare matlab versions for non-double variables */

double    *d_matbeg;
double    *d_matcnt;
double    *d_matind;
double    d_numrows;
double    d_numcols;
double    d_numnz;

/* Declare usual counters */
int kc;

/* Load rhs arguments */
M          =  mxGetPr(prhs[0]);
d_matbeg   =  mxGetPr(prhs[1]);
d_matcnt   =  mxGetPr(prhs[2]);
d_matind   =  mxGetPr(prhs[3]);
matval     =  mxGetPr(prhs[4]);
m          =  mxGetPr(prhs[5]);
d_numrows  = *mxGetPr(prhs[6]);
d_numcols  = *mxGetPr(prhs[7]);
d_numnz    = *mxGetPr(prhs[8]);
tol        = *mxGetPr(prhs[9]);

/* Convert dimensions to integers */
numrows = (int) (d_numrows + 0.1);
numcols = (int) (d_numcols + 0.1);
numnz   = (int) (d_numnz + 0.1);

/* Allocate storage for vectors */
matbeg = mxCalloc(numcols,sizeof(int));
matcnt = mxCalloc(numcols,sizeof(int));
matind = mxCalloc(numnz,sizeof(int));
sense  = mxCalloc(numrows,sizeof(char));
lb     = mxCalloc(numcols,sizeof(double));
ub     = mxCalloc(numcols,sizeof(double));

/* Set CPLEX data */
objsen =  CPX_MAX;

for(kc = 0; kc <= numnz-1; kc++)
  matind[kc] = (int) d_matind[kc] + 0.1;

for(kc = 0; kc <= numcols-1; kc++)
  {
    matbeg[kc] = (int) (d_matbeg[kc] + 0.1);
    matcnt[kc] = (int) (d_matcnt[kc] + 0.1);
	lb[kc] = -CPX_INFBOUND;
    ub[kc] =  CPX_INFBOUND;
  }

for(kc = 0; kc <= numrows-1; kc++)
  {
    sense[kc] = 'L';
  }

/* Create lhs matrices */

plhs[0] = mxCreateDoubleMatrix(numrows, 1, mxREAL);
plhs[1] = mxCreateDoubleMatrix(numrows, 1, mxREAL);
plhs[2] = mxCreateDoubleMatrix(numrows, 1, mxREAL);

costs = mxGetPr(plhs[0]);
stats = mxGetPr(plhs[1]);
ikeep = mxGetPr(plhs[2]);

/* Solve lp's*/

reducemex(M,costs,stats,ikeep,
          m,
          matbeg,matcnt,matind,matval,
          sense,objsen,lb,ub,
          numrows,numcols,tol);
}

reducemex(double M[], double costs[], double stats[], double ikeep[],
          double m[],
		  int matbeg[], int matcnt[], int matind[], double matval[],
          char sense[], int objsen, double lb[], double ub[],
          int numrows, int numcols, double tol)

{
CPXENVptr     env = NULL;
CPXLPptr      lp = NULL;
int           status;
char          probname[16];
double        *pi;
double        *slack;
double        *dj;


double        *obj;
int           solstat;
double        objval;
double        *x;

double        newval;
int           *indices;

int kc, i;

obj = mxCalloc(numcols,sizeof(double));
x   = mxCalloc(numcols,sizeof(double));

indices = mxCalloc(numcols,sizeof(int));

pi = mxCalloc(numrows,sizeof(double));
slack = mxCalloc(numrows,sizeof(double));
dj = mxCalloc(numcols,sizeof(double));

   /* Initialize the CPLEX environment */

   env = CPXopenCPLEX (&status);

   if ( env == NULL ) {
   char  errmsg[1024];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      goto TERMINATE;
   }

   /* Turn off output to the screen */

   status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
   if ( status ) {
      fprintf (stderr,
               "Failure to turn on screen indicator, error %d.\n", status);
      goto TERMINATE;
   }

   /* Turn off pre-processor

   status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);
   if ( status ) {
      fprintf (stderr,
               "Failure to turn off pre-processor, error %d.\n", status);
      goto TERMINATE;
   }

   /* Turn off aggregator

   status = CPXsetintparam (env, CPX_PARAM_AGGIND, CPX_OFF);
   if ( status ) {
      fprintf (stderr,
               "Failure to turn of aggregator, error %d.\n", status);
      goto TERMINATE;
   }*/

   /* Create the problem. */

   strcpy(probname,"LPnoname");
   lp = CPXcreateprob (env, &status, probname);

   if ( lp == NULL ) {
      fprintf (stderr, "Failed to create LP.\n");
      goto TERMINATE;
   }

   /* Initialize data */
   for (kc = 0; kc <= numcols-1; kc++)
     {
       indices[kc] = kc;
     }

   /* Now copy the problem data into the lp */

   status = CPXcopylp (env, lp, numcols, numrows, objsen, obj, m,
                       sense, matbeg, matcnt, matind, matval,
                       lb, ub, NULL);

   if ( status ) {
      fprintf (stderr, "Failed to copy problem data.\n");
      goto TERMINATE;
   }

/**********MAIN LOOP************/

   for(i = 0; i <= numrows-1; i++)
   {
   for (kc = 0; kc <= numcols-1; kc++)
     {
       obj[kc] = M[i+kc*numrows];
     }

   CPXchgobj(env, lp, numcols, indices, obj);

 
   newval = m[i] + 10;
   CPXchgrhs(env, lp, 1, &i, &newval);
 
   /* Optimize the problem and obtain solution. */

   status = CPXlpopt (env, lp);   /*status = CPX[prim/dual]opt (env, lp);*/
   if ( status ) {
      fprintf (stderr, "Failed to optimize LP.\n");
      goto TERMINATE;
   }

   status = CPXsolution (env, lp, &solstat, &objval, x, pi, slack, dj);
   if ( status ) {
      fprintf (stderr, "Failed to obtain solution.\n");
      goto TERMINATE;
   }

   ikeep[i] = 0;
   if (objval > m[i] + tol)  /* keep row */
     {
       ikeep[i] = 1;
       newval = m[i];
       CPXchgrhs(env, lp, 1, &i, &newval);
     }


/*   Write the output to the screen.  */
/*
   printf ("\nSolution status = %d\n", solstat);
   printf ("Solution value  = %f\n\n", objval);
*/

   costs[i] = objval;
   stats[i] = solstat;

}

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

      if ( status ) {
      char  errmsg[1024];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
      }
   }

}

