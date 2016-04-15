/* function which generates the MH chain */
/* [draws, acceptance] = try_mh(theta, lnw, u) */

#include "mex.h"
#include <math.h>

/* The computational routine */
void try_mh(double *theta, double *draws, double *lnw, double *u, mwSize Nk, mwSize N)
// void try_mh(double *theta, double *draws, double *lnw, double *u, mwSize N, mwSize k, double *a)

// in matrices come as ROW vectors 
{
    //mwSize i; 
    // mwSize is a type that represents size values, such as array dimensions.
    // Use this function for cross-platform flexibility. 
    // By default, mwSize is equivalent to int in C. 
    int i, j, m, s;
    double e;
    const mwSize k=Nk/N;
    double tmp[k];
    
    /* initialize chain */
    for (j=0; j<k; j++) {
        tmp[j] = theta[j];
        draws[j] = theta[j];
    }
    
    /* iterate over the chain*/
    s = 0;
    for (i=1; i<N; i++){
        e = exp(lnw[i] - lnw[s]);
        if (e > 1) // min between e and 1
        { 
            e = 1;
        }
        
        m = i*k;
        // compare with random uniform 
        if (u[i] <= e)
        {
           for (j=0; j<k; j++)
           {
               tmp[j] = theta[m+j];
           }
           s = i; // move the pointer
//            *a++; // increase the counter fo the acceptance rate;
        }
        // assign tmp to draws 
        for (j=0; j<k; j++)
        {
            draws[m+j] = tmp[j];
        }
    }
}



//////////////////////////////////////////////////////////////////////////


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
//     double multiplier;                       /* input scalar */
    double *inDraws, *inWeights, *inUnif;    /* 1 x inlen, 1 x nrows, 1 x nrows, input matrices (vectors) */
    size_t inlen, nrows;              /* sizes, ncols = inlen/nrows  */
//     size_t inlen, nrows, ncols;              /* sizes, ncols = inlen/nrows  */
    double *outDraws;                       /* output matrix */
//     double *outAcc;                           /* output acceptance rate */
                                     
                                     
//     /* check for proper number of arguments */
//     if(nrhs!=2) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
//     }
//     if(nlhs!=1) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
//     }
//     /* make sure the first input argument is scalar */
//     if( !mxIsDouble(prhs[0]) || 
//          mxIsComplex(prhs[0]) ||
//          mxGetNumberOfElements(prhs[0])!=1 ) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
//     }
//     
//     /* make sure the second input argument is type double */
//     if( !mxIsDouble(prhs[1]) || 
//          mxIsComplex(prhs[1])) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
//     }
//     
//     /* check that number of rows in second input argument is 1 */
//     if(mxGetM(prhs[1])!=1) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
//     }
    
    
    // step 1/3: in
    /* get the value of the scalar input  */
//     multiplier = mxGetScalar(prhs[4]);

    /* create a pointer to the real data in the input matrices */
    inDraws = mxGetPr(prhs[0]);
    inWeights = mxGetPr(prhs[1]);
    inUnif = mxGetPr(prhs[2]);
    
    /* get dimensions of the input matrices */
    inlen = mxGetN(prhs[0]);
    nrows = mxGetN(prhs[1]);
//     ncols = inlen/nrows; 
    
    // step 2/3: out
    /* create the output matrix of draws */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)inlen,mxREAL);
    /* get a pointer to the real data in the output matrix */
    outDraws = mxGetPr(plhs[0]);

    /* create the pointer to the scalar for acceptance rate*/ //????
//     plhs[1] = mxCreateDoubleScalar;
//     outAcc =  mxGetPr(plhs[1]);
    
    // step 3/3
    /* call the computational routine */
   // arrayProduct(multiplier,inMatrix,outMatrix,(mwSize)ncols);
    try_mh(inDraws, outDraws, inWeights, inUnif, (mwSize)inlen, (mwSize)nrows);
//     try_mh(inDraws, outDraws, inWeights, inUnif, (mwSize)nrows, (mwSize)ncols, outAcc)

}
