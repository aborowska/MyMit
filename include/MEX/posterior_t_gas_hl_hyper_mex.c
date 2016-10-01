#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI  3.1415926535897932;

/******************************************************************* */

void predict_t_gas(double *theta, double *y, mwSignedIndex N, mwSignedIndex hp, mwSignedIndex T, double *PL_hp)
{
    double *rho, rhoh ; 
    double *h, *y_hp ; 
    double *A, *nu_con;
    double tmp;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    rho = mxMalloc((N)*sizeof(double));
    A = mxMalloc((N)*sizeof(double));
    nu_con = mxMalloc((N)*sizeof(double));
    
    h = mxMalloc((N)*sizeof(double));
    y_hp = mxMalloc((N)*sizeof(double));
    
    /* Initialise */
    for (i=0; i<N; i++)
    {
        rho[i] = (theta[4*N+i]-2)/theta[4*N+i];
        A[i] = theta[i+2*N]*(theta[i+4*N]+3)/theta[i+4*N];
        nu_con[i] = (theta[i+4*N]+1)/(theta[i+4*N]-2);
        PL_hp[i] = 0;
        y_hp[i] = y[T-1];
    }    
     
    
    /* volatility_t_garch to get h_T*/   
    for (i=0; i<N; i++) 
    {
        for (j=1; j<T; j++)
        {
            tmp = (y[j-1]-theta[i])*(y[j-1]-theta[i]);
            tmp = tmp/(h[i]*(theta[i+4*N]-2));
            tmp = 1 + tmp;
            tmp = nu_con[i]/tmp;
            h[i] = theta[3*N+i]*h[i] + theta[i+N] + A[i]*(tmp*(y[j-1]-theta[i])*(y[j-1]-theta[i]) - h[i]);
        }                  
     }
    
    /* predict: h, y and pl */
    for (i=0; i<N; i++) 
    {
        for (j=1; j<(hp+1); j++)
        {   
             /* if j==1 then forecast based on the last observation */
             /* else forcast based on the previous forecast */   
            tmp = (y_hp[i]-theta[i])*(y_hp[i]-theta[i]);
            tmp = tmp/(h[i]*(theta[i+4*N]-2));
            tmp = 1 + tmp;
            tmp = nu_con[i]/tmp;
            h[i] = theta[3*N+i]*h[i] + theta[i+N] + A[i]*(tmp*(y_hp[i]-theta[i])*(y_hp[i]-theta[i]) - h[i]);                      
            rhoh = rho[i]*h[i];
            y_hp[i] = theta[i] + pow(rhoh,0.5)*theta[N*(j+4)+i]; 
            PL_hp[i] = PL_hp[i] + y_hp[i];
        } 
        PL_hp[i] = 100*(exp(0.01*PL_hp[i])-1);
     }
    
    /* Free allocated memory */
    mxFree(rho);
    mxFree(h);
    mxFree(y_hp);
    mxFree(A); 
    mxFree(nu_con); 
}

/******************************************************************* */

void duvt_garch(double x, double mu, double sigma, double df, double *GamMat, mwSignedIndex G, double *pdf)
{
    double c0, c1, c2, c, e, tmp, etmp, df5;
    int ind;
    
    c0 = df+1;
    
    if (c0 <= 100) /* gamma((nu+1)/2) */
    {
        ind = floor(c0*50000) - 1;
        c1 = GamMat[ind];
    }
    else
    {
        c1 = GamMat[G-1];
    }

    if (df <= 100) /* gamma(nu/2) */
    {
        df5 = df*50000;
        if ((df5 - floor(df5)) < (floor(df5+1) - df5))
        {
            ind = floor(df5);
        }
        else
        {
            ind = floor(df5 + 1);
        }
        ind = ind -  1;
        c2 = GamMat[ind];
    }
    else
    {
        c2 = GamMat[G-1];
    } 

    c = df*PI;
    c = pow(c,0.5);
    c2 = c2*c;

    c2 = c2*pow(sigma,0.5);

    c = c1/c2;

    e = -0.5*c0; 
    
    tmp = (x-mu)*(x-mu)/sigma;
    tmp = 1 + tmp/df;
    etmp = pow(tmp,e); 
    pdf[0] = exp(log(c) + log(etmp));
}

/******************************************************************* */

void prior_t_gas_hl_hyper(double *theta, double *hyper, double *PL_hp, double *VaR,
        mwSignedIndex N, mwSignedIndex *r1, double *r2)
{
    mwSignedIndex i;
    
    /* Variable size arrays */
    
    for (i=0; i<N; i++)
    {
        r1[i] = 1;        
        r2[i] = -mxGetInf();

        if (theta[i+N] <= 0) // omega>0
        {
            r1[i] = 0;
        }   
        if ((theta[i+3*N] < 0) || (theta[i+3*N] >= 1)) // 0<=B<1
        {
            r1[i] = 0;
        }
        if (theta[i+4*N] <= 2) //nu>2
        {
            r1[i] = 0;
        }     
        if (PL_hp[i] > VaR[0])
        {
            r1[i] = 0;
        }  
        if (r1[i] == 1)
        {
            r2[i] = log(hyper[0]) - hyper[0]*(theta[i+4*N] - 2);
        }
    }
}

/*******************************************************************  */

void posterior_t_gas_hl_hyper_mex(double *y, mwSignedIndex N, mwSignedIndex hp, mwSignedIndex T,
        double *theta, double *hyper, double *VaR, double *GamMat, mwSignedIndex G, double *d)
{
    mwSignedIndex *r1;
    double *r2; 
    double *PL_hp;
    double *rho, *A, *nu_con;
    double h, *pdf; 
    double rhoh, tmp;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));              
    rho = mxMalloc((N)*sizeof(double));
    A = mxMalloc((N)*sizeof(double));
    nu_con = mxMalloc((N)*sizeof(double));
    pdf = mxMalloc((1)*sizeof(double));
    PL_hp = mxMalloc((N)*sizeof(double)); 
    
    
    predict_t_gas(theta, y, N, hp, T, PL_hp);  
    prior_t_gas_hl_hyper(theta, hyper, PL_hp, VaR, N, r1, r2);

    /* Initialise */
    for (i=0; i<N; i++)
    {   
        rho[i] = (theta[i+4*N]-2)/theta[i+4*N];
        A[i] = theta[i+2*N]*(theta[i+4*N]+3)/theta[i+4*N];
        nu_con[i] = (theta[i+4*N]+1)/(theta[i+4*N]-2);
    }
        
    /* PDF */
    for (i=0; i<N; i++) 
    {          
      
        if (r1[i]==1)
        {
            d[i] = r2[i];
            h = theta[i+N]/(1-theta[i+3*N]); 
            
            for (j=1; j<T; j++)
            {   
                tmp = (y[j-1]-theta[i])*(y[j-1]-theta[i]);
                tmp = tmp/(h*(theta[i+4*N]-2));
                tmp = 1 + tmp;
                tmp = nu_con[i]/tmp;
                h = theta[3*N+i]*h + theta[i+N] + A[i]*(tmp*(y[j-1]-theta[i])*(y[j-1]-theta[i]) - h);   
                rhoh = rho[i]*h;  
                duvt_garch(y[j], theta[i], rhoh, theta[i+4*N], GamMat, G, pdf);
                d[i] = d[i] + log(pdf[0]);
            } 

            for (j=1; j<hp+1; j++) /* pdf of future disturbances */
            {
                duvt_garch(theta[i+N*(j+4)], 0, 1, theta[i+4*N], GamMat, G, pdf);
                pdf[0] = log(pdf[0]);
                d[i] = d[i] + pdf[0];
            }
        }
        else
        {
            d[i] = -mxGetInf();
        }    
     }
    
    /* Free allocated memory */
    mxFree(r1); 
    mxFree(r2); 
    mxFree(rho); 
    mxFree(A); 
    mxFree(nu_con); 
    mxFree(pdf);
    mxFree(PL_hp);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, hp, T, G;                          /* size of matrix */
    double *y, *theta, *hyper, *VaR;                /* input*/
    double *GamMat;                                     /* input */
    double *d;                                          /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    hyper = mxGetPr(prhs[2]);
    VaR = mxGetPr(prhs[3]);
    GamMat = mxGetPr(prhs[4]);
   
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    hp = mxGetN(prhs[0]); /* length of prediction horizon */
    hp = hp - 5;
    T = mxGetM(prhs[1]); /* no. of observations */
    G = mxGetM(prhs[4]);
        
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    
    /* call the function */
    posterior_t_gas_hl_hyper_mex(y, N, hp, T, theta, hyper, VaR, GamMat, G, d);
  
}
