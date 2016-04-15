#include "mex.h"

/* The computational routine */
// double *theta, int *N, int *k, double *lnw, double *u, double *draws, int *ns)
// {
//   int m, s;
//   double e;
//   double tmp[*k];
void play_mh(double *x, double *y, double *u, double *z, int n)
{
    int i;
    /* "inverse" x to y ad create z=max(x,u)*/
    for (i=0; i<n; i++) {
        y[i] = x[n-i-1];
        if (x[i] < u[i]) {
            z[i] = u[i];
        }
        else {
            z[i] = x[i];
        }
    }
}

//#include "matrix.h"
// size_t mxGetM(const mxArray *pm); pm - Pointer to an mxArray, returns number of rows in the mxArray to which pm points.
// size_t mxGetN(const mxArray *pm); pm - Pointer to an mxArray, returns number of columns in the mxArray to which pm points.

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *inMatrix1, *inMatrix2;               /* 1xN input matrix */
    size_t ncols;                   /* size of matrix */
    //size_t is an unsigned data type defined by several C/C++ standards; this type is used to represent the size of an object
    double *outMatrix1, *outMatrix2;              /* output matrix */

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
    
//     /* get the value of the scalar input  */
//     multiplier = mxGetScalar(prhs[0]);

    // step 1/3: in
    /* create a pointer to the real data in the input matrix  */
    inMatrix1 = mxGetPr(prhs[0]);
    inMatrix2 = mxGetPr(prhs[1]);

    /* get dimensions of the input matrix */ 
    ncols = mxGetN(prhs[0]); //lenghtof theinput vector

    // step 2/3: out
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL); // allocate memory
    plhs[1] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL); // allocate memory

    /* get a pointer to the real data in the output matrix */
    outMatrix1 = mxGetPr(plhs[0]);
    outMatrix2 = mxGetPr(plhs[1]);

    // step 3/3
    /* call the computational routine */
    play(inMatrix1,outMatrix1,inMatrix2,outMatrix2,(mwSize)ncols);
}
