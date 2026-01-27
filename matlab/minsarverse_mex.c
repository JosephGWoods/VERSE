/* Adapted from Brian Hargreaves' mintverse code      */
/* by Joseph G. Woods, May 2022 */

#define MEX_COMPILE

/*
#define DEBUG
*/

#include <mex.h>
#include "verse.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* function [bv, gv] = mintverse_mex(b, g, dt, gmax, smax) */
{
    double *br;  /* real part of B1 pulse (arbitrary units) */
    double *bi;  /* imag part of B1 pulse (arbitrary units) */
    double *g;   /* gradient pulse (arbitrary units)  */
    double dt;   /* time increment (arbitrary units) */
    int n;       /* number of points */
    double gmax; /* maximum gradient amplitude, (units of g) */
    double smax; /* maximum slew rate, (units of g per unit of dt) */

    double *brwork; /* output B1 pulse (arbitrary units) */
    double *biwork; /* output B1 pulse (arbitrary units) */
    double *gwork;  /* output Gradient pulse (arbitrary units) */

    double *brout; /* Matlab output B1 pulse (real) */
    double *biout; /* Matlab output B1 pulse (imag) */
    double *gout;  /* Matlab output gradient pulse. */

    int count;        /* For loops         */
    int bcomplex = 0; /* 1 if b is complex */

    /* ================ Get input parameters ================ */

    br = mxGetPr(prhs[0]);
    n = (int) (mxGetM(prhs[0]) * mxGetN(prhs[0]));
    if (mxIsComplex(prhs[0])) {
        bi = mxGetPi(prhs[0]);
        bcomplex = 1;
    } else {
        bi = (double *) malloc(n*sizeof(double));
        for (count=0; count <n; count++)
            bi[count] = 0.0;
        bcomplex = 0;
    }

    g = mxGetPr(prhs[1]);
    if ( (int) (mxGetM(prhs[1])* mxGetN(prhs[1])) != n)
        mexErrMsgTxt("B1 and Gradient should be the same length.");

    dt   = *mxGetPr(prhs[2]);
    gmax = *mxGetPr(prhs[3]);
    smax = *mxGetPr(prhs[4]);

    /* Allocate outputs */
    brwork = (double *) malloc(n*sizeof(double));
    biwork = (double *) malloc(n*sizeof(double));
    gwork  = (double *) malloc(n*sizeof(double));

    /* Call C-function here */
    minsarverse(br,bi,g,dt,n,gmax,smax,brwork,biwork,gwork);

    /* Allocate Matlab outputs, and copy function outputs to them. */
    if (bcomplex==0)
        plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    else
        plhs[0] = mxCreateDoubleMatrix(n,1,mxCOMPLEX);

    plhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
    brout = mxGetPr(plhs[0]);
    if (bcomplex==1)
        biout = mxGetPi(plhs[0]);
    gout = mxGetPr(plhs[1]);

    #ifdef DEBUG
    printf("mex function freeing memory. \n");
    #endif

    memcpy(brout,brwork,n*sizeof(double));
    if (bcomplex==1)
        memcpy(biout,biwork,n*sizeof(double));
    memcpy(gout,gwork,n*sizeof(double));

    /*	Free up space that was allocated here and in mintverse(). */
    if (bcomplex==0)
        free(bi);
    free(brwork);
    free(biwork);
    free(gwork);
}
