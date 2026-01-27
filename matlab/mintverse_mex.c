/* Adapted from Brian Hargreaves' mintverse code       */
/* by Joseph G. Woods, University of Oxford, June 2022 */

#define MEX_COMPILE

/*
#define DEBUG
*/

#include <mex.h>
#include "verse.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/*  function [bv, gv] = mintverse_mex(b1, g, dt, bmax, gmax, smax, emax)*/
{
    double *br;   /* real part of B1 pulse (arbitrary units)   */
    double *bi;   /* imag part of B1 pulse (arbitrary units)   */
    double *g;    /* gradient pulse (arbitrary units)          */
    double dt;    /* time increments (arbitrary units)         */
    long n;       /* number of points                          */
    double bmax;  /* maximum b1-field amplitude (units of b1)  */
    double gmax;  /* maximum gradient amplitude, (units of g)  */
    double smax;  /* maximum slew rate, (units of g per unit of dt) */
    double emax;  /* maximum energy (optional)                 */

    double *brwork;	/* output B1 pulse (arbitrary units)       */
    double *biwork;	/* output B1 pulse (arbitrary units)       */
    double *gwork;	/* output Gradient pulse (arbitrary units) */

    long nout;      /* number of output points       */
    double *brout;	/* matlab Output B1 pulse (real) */
    double *biout;	/* matlab Output B1 pulse (imag) */
    double *gout;	/* matlab Output gradient pulse  */

    long count;       /* for loops         */
    int bcomplex = 0; /* 1 if b is complex */

    /*================	Get required input parameters================	*/

    br = mxGetPr(prhs[0]);
    n = (long) (mxGetM(prhs[0]) * mxGetN(prhs[0]));
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
    if ( (long) (mxGetM(prhs[1])* mxGetN(prhs[1])) != n) {
        printf("Error: B1 and Gradient should be the same length. \n");
        return;
    }

    dt   = *mxGetPr(prhs[2]);
    bmax = *mxGetPr(prhs[3]);
    gmax = *mxGetPr(prhs[4]);
    smax = *mxGetPr(prhs[5]);

    /*================ Get optional input parameters  ===============	*/
    /* If a default is used, a message is displayed. */

    if (nrhs > 6) {
        emax = *mxGetPr(prhs[6]);
    } else {
        emax = -1.0; /* Maximum pulse energy */
        printf("No max B1 energy constraint given.\n");
    }

    /* Call C-function here */
    mintverse(br,bi,g,dt,n,bmax,gmax,smax,emax,&nout,&brwork,&biwork,&gwork);

    /*	Allocate Matlab outputs, and copy function outputs to them. */
    if (bcomplex==0) plhs[0] = mxCreateDoubleMatrix(nout,1,mxREAL);
    else             plhs[0] = mxCreateDoubleMatrix(nout,1,mxCOMPLEX);
    brout = mxGetPr(plhs[0]);
    if (bcomplex==1)
        biout = mxGetPi(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix(nout,1,mxREAL);
    gout = mxGetPr(plhs[1]);

    #ifdef DEBUG
    printf("mex function freeing memory. \n");
    #endif

    memcpy(brout,brwork,nout*sizeof(double));
    if (bcomplex == 1)
        memcpy(biout,biwork,nout*sizeof(double));
    memcpy(gout,gwork,nout*sizeof(double));

    /*	Free up space that was allocated here and in mintverse(). */
    if (bcomplex==0)
        free(bi);
    free(brwork);
    free(biwork);
    free(gwork);
}




