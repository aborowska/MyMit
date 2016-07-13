#include "mex.h"
#include "math.h"
#include "matrix.h"

 
/******************************************************************* */
 void predict_t_garch_noS_mex(double *theta, double *y_T, double *h_T, double *eps, mwSignedIndex N, mwSignedIndex hp, double *y_hp)
{
    double *rho, rhoh ; 
    double *h; 
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    rho = mxMalloc((N)*sizeof(double));
    h = mxMalloc((N)*sizeof(double));
     
    /* Initialise */
    for (i=0; i<N; i++)
    {
        rho[i] = (theta[4*N+i]-2)/theta[4*N+i];
        h[i] = h_T[i]; 
     }    
    
    /* predict: h and y*/
    for (i=0; i<N; i++) 
    {
        for (j=0; j<hp; j++)
        {   
            if (j==0)
            {
                h[i] = theta[2*N+i]*h[i] + theta[i] + theta[N+i]*(y_T[0]-theta[3*N+i])*(y_T[0]-theta[3*N+i]); /* forecast based on the last observation */
            }
            else
            {
                h[i] = theta[2*N+i]*h[i] + theta[i] + theta[N+i]*(y_hp[i+N*(j-1)]-theta[3*N+i])*(y_hp[i+N*(j-1)]-theta[3*N+i]);  /* forcast based on the previous forecast */                      
            }

            rhoh = rho[i]*h[i];
            y_hp[i+N*j] = theta[3*N+i] + pow(rhoh,0.5)*eps[i+N*j]; 
         } 

      }
    
    /* Free allocated memory */
    mxFree(rho);
    mxFree(h);
}



/******************************************************************* */
 

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, hp;                         /* size of matrix */
    double *y_T, *theta, *h_T, *eps;             /* input */
    double *y_hp;                                /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y_T = mxGetPr(prhs[1]); // the last observation (double)
    h_T = mxGetPr(prhs[2]); // computed last volatility (array of doubles of length N)
    eps = mxGetPr(prhs[3]); // a matrix of future disturbances (N)x(hp)

    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    hp = mxGetN(prhs[3]); /* length of prediction horizon */
  
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,hp,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    y_hp = mxGetPr(plhs[0]);
    
    /* call the function */
    predict_t_garch_noS_mex(theta, y_T, h_T, eps, N, hp, y_hp);
  
}
