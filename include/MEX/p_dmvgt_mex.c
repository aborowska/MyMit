#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "blas.h"
#include "lapack.h"

#define  PI  3.1415926535897932;

void dmvt_mex(double *x, double *mu, double *Sigma, double *df,
        double *GamMat,  mwSignedIndex d,  mwSignedIndex G, mwSignedIndex N, 
        double *dens)
/* density of the multivariate Student's t distribution*/
{
    mwSignedIndex i, j, info;
    mwSignedIndex one = 1;
    mwSignedIndex  *IPIV; /* pointer to ints */
    int ind;
    double *xmu;
    double c0, c1, c2, c, e, tmp, etmp, df5;
    double sqrt_det_Sigma;
    char *cht = "T";
    
    /* variable size arrays */  
    xmu = mxMalloc((d)*sizeof(double));    
    IPIV = mxMalloc((d)*sizeof(mwSignedIndex)); 
       
//     mexPrintf("**** IN ****\n");

    
    /* compute the constants */
    c0 = df[0] + d;
    
    if (c0 <= 100)
    {
        ind = floor(c0*50000) - 1;
        c1 = GamMat[ind];
    }
    else
    {
        c1 = GamMat[G-1];
    }
    
    if (df[0] <= 100)
    {
//         ind = floor(df[0]*50000) - 1; /* useround insteadas follows: */
        df5 = df[0]*50000;
        if ((df5 - floor(df5)) < (floor(df5+1) - df5))
        {
            ind = floor(df5);
        }
        else
        {
            ind = floor(df5 + 1);
        }
//     	mexPrintf("ind df = %i\n",ind);       
        ind = ind -  1;
        c2 = GamMat[ind];
    }
    else
    {
        c2 = GamMat[G-1];
    }      
    
    /* Cholesky decomposition of Sigma */
    /* & changes a value into a pointer*/ 
    dgetrf(&d, &d, Sigma, &d, IPIV, &info); 
    
    /* determinant of Sigma from LU factorisation*/
    sqrt_det_Sigma = 1;
    for (i=0; i<d; i++)
    {
        sqrt_det_Sigma = sqrt_det_Sigma*Sigma[i+d*i];
    }
    sqrt_det_Sigma = pow(fabs(sqrt_det_Sigma),0.5);
//     mexPrintf("sqrt_det_Sigma = %6.4f\n",sqrt_det_Sigma);

    c = df[0]*PI;
    c = pow(c,0.5*d);
//     mexPrintf("c = %6.4f\n",c);    
    c2 = c2*c;
    c2 = c2*sqrt_det_Sigma;
    c = c1/c2;
    e = -0.5*c0;    
//     mexPrintf("e = %6.4f\n",e);

    for (j=0; j<N; j++)
    {
        for (i=0; i<d; i++)
        {
            xmu[i] = x[j+N*i] - mu[i];
        }
        /* inv(Sigma)*xmu'  */
        /* on exit, xmu is the solution X to Sigma*X = xmu i.e. is equal to xmu*inv(Sigma)*/
         dgetrs(cht, &d, &one, Sigma, &d, IPIV, xmu, &d, &info); /* Sigma is the lu factorisation of Sigma */

        tmp = 0; /* tmp = (x-mu)*inv(Sigma)*(x-mu)' */
        for (i=0; i<d; i++)
        {
            tmp = tmp + xmu[i]*(x[j+N*i] - mu[i]); /* xmu = inv(Sigma)*(x-mu)' */
//             mexPrintf("tmp[%i] = %6.4f\n",i,tmp);
        } 
        tmp = 1 + tmp/df[0];
//         mexPrintf("tmp = %16.14f\n",tmp);        
        etmp = pow(tmp,e);
//         mexPrintf("etmp = %16.14f\n",etmp);
        dens[j] = exp(log(c) + log(etmp));
//         mexPrintf("dens[%i] = %16.14f\n",j,dens[j]);

    }

//     mexPrintf("**** OUT ****\n");

    mxFree(IPIV); 
    mxFree(xmu);
}

void p_dmvgt_mex(double *x, double *beta, double *Sigma, double *df, double *p,
        double *X,
        double *GamMat, double *L,
        mwSignedIndex d, mwSignedIndex r, mwSignedIndex G, 
        mwSignedIndex N, mwSignedIndex H,  
        double *dens)
/* density of the mixture of the multivariate Studenet's t distributions */        

{
    mwSignedIndex i, j, k, m, one;
    double *mu_h, *Sigma_h, *df_h, *dens_h, *x_tmp;    

    df_h = mxMalloc((1)*sizeof(double));    
    mu_h = mxMalloc((d)*sizeof(double));    
    Sigma_h = mxMalloc((d*d)*sizeof(double));    
    dens_h = mxMalloc((1)*sizeof(double));
    x_tmp = mxMalloc((d)*sizeof(double));
    
    one = 1;
    /* initialitsation */
    for (i=0; i<N; i++)
    {
        dens[i] = 0;
    }
    
    /* loop over components */
     for (i=0; i<H; i++)
     {
//           mexPrintf("\n*** h = %i ***\n",i+1);

          df_h[0] = df[i];
//           mexPrintf("df_h = %4.2f\n", df_h[0]);

         for (j=0; j<(d*d); j++)
         {
             Sigma_h[j] = Sigma[i+H*j];
//              mexPrintf("Sigma_h[%i] = %6.4f\n",j,Sigma_h[j]);
         }
              
         for (j=0; j<N; j++)
         {
            /* determine the mean of the j-th draw*/
            for (k=0; k<d; k++)
            {
                x_tmp[k] = x[j+N*k];
                mu_h[k] = 0;
                for (m=0; m<r; m++)
                {
                  mu_h[k] =  mu_h[k] + X[j+N*m]*beta[i+H*(k*r+m)];
                }
//                 mexPrintf("x_tmp[%i,%i] = %6.4f\n",j,k,x_tmp[k]);   
//                 mexPrintf("mu_h[%i,%i] = %6.4f\n",j,k,mu_h[k]);   
            }
            dmvt_mex(x_tmp, mu_h, Sigma_h, df_h, GamMat, d, G, one, dens_h);  /* density on the i-th component*/

//              mexPrintf("dens_h[%i] = %16.14f\n",j,dens_h[j]);
//              mexPrintf("p[%i] = %16.14f\n",i,p[i]);
             dens[j] = dens[j] + exp(log(p[i]) + log(dens_h[0]));
//              mexPrintf("dens[%i] = %16.14f\n",j,dens[j]);

         }
     }

//     mexPrintf("L = %4.2f\n", L[0]);    
    if (L[0]==1)
    {
//         mexPrintf("YO!\n");
        for (i=0; i<N; i++)
        {
            dens[i] = log(dens[i]);
        }
    }
        
    mxFree(df_h); mxFree(mu_h); mxFree(Sigma_h); mxFree(dens_h); mxFree(x_tmp);
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *L; /*log or not*/
    mwSignedIndex G, d, N, H, r;                          /* size of matrix */
    double *x, *beta, *Sigma, *X, *GamMat, *df, *p;           /* input*/
    double *dens;                                      /* output */
    
    /* Getting the inputs */
    x = mxGetPr(prhs[0]);
    beta = mxGetPr(prhs[1]); /* mu = beta*X;    beta: (r)x(d), X: (N)x(r)*/
    Sigma = mxGetPr(prhs[2]);
    df = mxGetPr(prhs[3]);
    p = mxGetPr(prhs[4]);
    X = mxGetPr(prhs[5]);
    GamMat = mxGetPr(prhs[6]); 
    L = mxGetPr(prhs[7]); 
    
    N = mxGetM(prhs[0]); /* number of draws */
    d = mxGetN(prhs[1]); /* the dimension of t-distribution */
    H = mxGetM(prhs[1]); /* number of components */
    r = mxGetN(prhs[5]); /* number of regressors in X */
    G = mxGetM(prhs[6]); /* the length of GamMat */

    d = d/r; /* beta is stored as (1)x(d*r) vector */

//     mexPrintf("N = %i\n", N);    
//     mexPrintf("d = %i\n", d);    
//     mexPrintf("H = %i\n", H);    
//     mexPrintf("r = %i\n", r);    
//     mexPrintf("G = %i\n", G);    
      
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    
    /* get a pointer to the real data in the output matrix */
    dens = mxGetPr(plhs[0]);
           
    /* call the function */
    p_dmvgt_mex(x, beta, Sigma, df, p, X, GamMat, L, d, r, G, N, H, dens);  
    
}
