#include "mex.h"
#include "math.h"

const double M = -1e100;

void duvt_garch(y, theta[i+2*N], rho[i]*h[0], theta[i+3*N], GamMat, pdf[0])
{
}

void prior_t_garch(double *theta, mwSignedIndex N, int *r1, double *r2)
{
    mwSignedIndex i;
    
    /* Variable size arrays */
  /*  c1 = malloc((N)*sizeof(double));             */
    
    for (i=0; i<N; i++)
    {
        r1[i] = 1;
        if ((theta[i] <0 ) || (theta[i] >= 1))
        {
            r1[i] = 0;
        }
        if ((theta[i+N] < 0) || (theta[i+N] >= 1))
        {
            r1[i] = 0;
        }
        if (theta[i] + theta[i+N] >= 1)
        {
            r1[i] = 0;
        }
        if (theta[i+3*N] <= 2)
        {
            r1[i] = 0;
        }
        
        if (r1[i] == 1)
        {
            r2[i] = 2 - theta[i+3*N];
        }
    }
}

void posterior_t_garch(double *y, mwSignedIndex N, mwSignedIndex T, double *S,
        double *theta, double *GamMat, double *d)
{
    int *r1;
    double *r2; 
    double *omega, *rho;
    double *h, *pdf; 
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = malloc((N)*sizeof(int));              
    r2 = malloc((N)*sizeof(double));              
    omega = malloc((N)*sizeof(double));   
    rho = malloc((N)*sizeof((double));
    
    prior_t_garch(theta, N, r1, r2);
    
    /* Initialise */
    for (i=0; i<N; i++)
    {
        omega[i] = S[0]*(1-theta[i]-theta[N+i]);  
        rho[i] = (theta[i+3*N]-2)/theta(i+3*N);
    }    

         
    /* PDF */
    for (i=0; i<N; i++) 
    {
        if (r1[i]==1)
        {
            d[i] = r2[i];
            h[0] = S[0]; 
            for (j=1; j<T; j++)
            {   
    /*            h[i] = beta[i]*h[i] + omega[i] + alpha[i]*(y[j-1]-mu[i])*(y[j-1]-mu[i]);     */    
                h[0] = theta[N+i]*h[0] + omega[i] + theta[i]*(y[j-1]-theta[2*N+i])*(y[j-1]-theta[2*N+i]);         
                duvt_garch(y[j], theta[i+2*N], rho[i]*h[0], theta[i+3*N], GamMat, pdf[0]);
                d[i] = d[i] + log(pdf[0]);
            }         
        }
        else
        {
            d[i] = M;  
        }    
     }
    
    /* Free allocated memory */
    free(omega);  
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, T;                               /* size of matrix */
    double *y, *S, *theta;                            /* input*/
    double *GamMat;                                   /* output */
    
    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    S = mxGetPr(prhs[2]);
    GamMat = mxGetPr(prhs[3])
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    
    /* call the function */
    posterior_t_garch(y, N, T, S, theta, GamMat, d);
  
}
