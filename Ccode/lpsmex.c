/*Copyright 2015 Jeff Shamma
 Address:
 School of Electrical and Computer Engineering
 Georgia Institute of Technology
 777 Atlantic Dr NW
 Atlanta, GA 30332-0250
 */
#include "cplex.h"
#include "mex.h"
#include <string.h>
#include <math.h>

/* mex-file gateway

   [v,x,solstat] = lpsmex(c,matbeg,matcnt,matind,matval,
                          b,e,objsen,lb,ub,ncon,nvar,numnz)
*/


void
free_and_null (char **ptr);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])

{

/* Declare variables with "natural" cplex data type */
double    *objval;
double    *x;
int       solstat;
double    *obj;
int       *matbeg; /* columnwise beginning index    */
int       *matcnt; /* columnwise number of elements */
int       *matind; /* columnwise row indices        */
double    *matval; /* values */
double    *rhs;
char      *sense;
int       objsen;
double    *lb;
double    *ub;
int       numrows;
int       numcols;
int       numnz;

/* Declare matlab versions for non-double variables */
double    *d_matbeg;
double    *d_matcnt;
double    *d_matind;
double    *d_sense;
double    d_objsen;
double    d_numrows;
double    d_numcols;
double    d_numnz;

/* Declare usual counters */
int kc;

/* Load rhs arguments */
obj        =  mxGetPr(prhs[0]);
d_matbeg   =  mxGetPr(prhs[1]);
d_matcnt   =  mxGetPr(prhs[2]);
d_matind   =  mxGetPr(prhs[3]);
matval     =  mxGetPr(prhs[4]);
rhs        =  mxGetPr(prhs[5]);
d_sense    =  mxGetPr(prhs[6]);
d_objsen   = *mxGetPr(prhs[7]);
lb         =  mxGetPr(prhs[8]);
ub         =  mxGetPr(prhs[9]);
d_numrows  = *mxGetPr(prhs[10]);
d_numcols  = *mxGetPr(prhs[11]);
d_numnz    = *mxGetPr(prhs[12]);

/* Convert dimensions to integers */
numrows = (int) (d_numrows + 0.1);
numcols = (int) (d_numcols + 0.1);
numnz   = (int) (d_numnz + 0.1);

/* Allocate storage for vectors */
matbeg = mxCalloc(numcols,sizeof(int));
matcnt = mxCalloc(numcols,sizeof(int));
matind = mxCalloc(numnz,sizeof(int));
sense  = mxCalloc(numrows,sizeof(char));

/* Convert rhs arguments to cplex data type */
if (d_objsen < 0)
  objsen =  CPX_MAX;
else
  objsen =  CPX_MIN;


for(kc = 0; kc <= numnz-1; kc++)
  matind[kc] = d_matind[kc] + 0.1;

for(kc = 0; kc <= numcols-1; kc++)
  {
    matbeg[kc] = (int) (d_matbeg[kc] + 0.1);
    matcnt[kc] = (int) (d_matcnt[kc] + 0.1);
  }

for(kc = 0; kc <= numrows-1; kc++)
  {
    if (d_sense[kc] < -0.5)
      sense[kc] = 'L';
    else if (d_sense[kc] < 0.5)
      sense[kc] = 'E';
    else if (d_sense[kc] < 1.5)
      sense[kc] = 'G';
    else
      sense[kc] = 'R';

/* Gurdal Arslan added on 7/30/2002 */
}

for(kc = 0; kc < numcols; kc++)
  {
/***********************************/
    if ( mxIsInf(-lb[kc]) )
      lb[kc] = -CPX_INFBOUND;
                                 /* CHANGES VALUES*/
    if ( mxIsInf(ub[kc]) )
      ub[kc] =  CPX_INFBOUND;
  }

/* Create lhs matrices */

plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
plhs[1] = mxCreateDoubleMatrix(numcols, 1, mxREAL);
plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);


/* Point to lhs */

objval = mxGetPr(plhs[0]);
x = mxGetPr(plhs[1]);

/* Solve lp*/

lpsmex(objval,x,&solstat,
         obj,matbeg,matcnt,matind,matval,
         rhs,sense,objsen,lb,ub,
         numrows,numcols);

*(mxGetPr(plhs[2])) = solstat;   /*  double <- int  */

}

lpsmex(
  double *objval, double x[], int *solstat,
  double  obj[], int matbeg[], int matcnt[], int matind[], double matval[],
  double rhs[], char sense[], int objsen, double lb[], double ub[],
  int numrows, int numcols)

{
CPXENVptr     env = NULL;
CPXLPptr      lp = NULL;
int           status;
char          probname[16];
double        *pi;
double        *slack;
double        *dj;

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
   
   /*fprintf(stderr, "TP 1");*/

   /* Turn off output to the screen */

   status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
   if ( status ) {
      printf ("Failure to turn off screen indicator, error %d.\n", status);
      goto TERMINATE;
   }
   /* Turn off pre-processor

   status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);
   if ( status ) {
      fprintf (stderr,
               "Failure to turn off pre-processor, error %d.\n", status);
      goto TERMINATE;
   }

   **** Turn off aggregator

   status = CPXsetintparam (env, CPX_PARAM_AGGIND, CPX_OFF);
   if ( status ) {
      fprintf (stderr,
               "Failure to turn of aggregator, error %d.\n", status);
      goto TERMINATE;
   

   **** Create the problem. */
   
   strcpy(probname,"LPnoname");
   lp = CPXcreateprob (env, &status, probname);

   if ( lp == NULL ) {
      printf ("Failed to create LP.\n");
      goto TERMINATE;
   }

   /* Now copy the problem data into the lp */

   /*
   status = CPXcheckcopylp (env, lp, numcols, numrows, objsen, obj,
                          rhs, sense, matbeg, matcnt, matind,
                          matval, lb, ub, NULL);

   if ( status ) {
      printf ("Failed before copying problem data.\n");
      goto TERMINATE;
   }
   
   
   printf("Check point 4.5\n");
   */
   status = CPXcopylp (env, lp, numcols, numrows, objsen, obj, rhs,
                       sense, matbeg, matcnt, matind, matval,
                       lb, ub, NULL);

   if ( status ) {
      printf ("Failed to copy problem data.\n");
      goto TERMINATE;
   }

   /* Optimize the problem and obtain solution. */

   
  
   /*CPXsavwrite(env, lp, "test.sav");*/
   status = CPXlpopt (env, lp);
/*   status = CPX[prim/dual]opt (env, lp); */
   if ( status ) {
       printf ("Failed to optimize LP.\n");
      fprintf (stderr, "Failed to optimize LP.\n");
      goto TERMINATE;
   }

   
   status = CPXsolution (env, lp, solstat, objval, x, pi, slack, dj);
   if ( status ) {
      fprintf (stderr, "Failed to obtain solution.\n");
      goto TERMINATE;
   }

/*
   Write the output to the screen.

   printf ("\nSolution status = %d\n", *solstat);
   printf ("Solution value  = %f\n\n", *objval);
*/

TERMINATE:

   /* Free up the problem as allocated by CPXcreateprob, if necessary */

   if ( lp != NULL ) {
      int  frstatus;
      frstatus = CPXfreeprob (env, &lp);
      if ( frstatus ) {
         fprintf (stderr, "CPXfreeprob failed, error code %d.\n", frstatus);
         if (( !status ) && frstatus )  status = frstatus;
      }
   }

   /* Free up the CPLEX environment, if necessary */

   if ( env != NULL ) {
      int  clstatus;
      clstatus = CPXcloseCPLEX (&env);

      if ( clstatus ) {
         fprintf (stderr, "CPXcloseCPLEX failed, error code %d.\n", clstatus);
         if (( !status ) && clstatus )  status = clstatus;
      }
   }

   if ( status ) {
      char  errmsg[1024];

      /* Note that since we have turned off the CPLEX screen indicator,
         we'll need to print the error message ourselves. */

      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
   }

    
    return (status);
}



/* This simple routine frees up the pointer *ptr, and sets *ptr to NULL*/

void
free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
}  /**/  
