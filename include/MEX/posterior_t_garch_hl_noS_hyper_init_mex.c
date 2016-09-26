#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI  3.1415926535897932;
// const double M = -1e100;

/******************************************************************* */

void predict_t_garch(double *theta, double *y, double *S, mwSignedIndex N, mwSignedIndex hp, mwSignedIndex T, double *PL_hp)
{
    double *rho, rhoh ; 
    double *h, *y_hp ; 
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    rho = mxMalloc((N)*sizeof(double));
    h = mxMalloc((N)*sizeof(double));
    y_hp = mxMalloc((N)*sizeof(double));
    
    /* Initialise */
    for (i=0; i<N; i++)
    {
        rho[i] = (theta[4*N+i]-2)/theta[4*N+i];
        h[i] = S[0]; 
//         mexPrintf("h[%i] = %6.4f \n",0,h[i]);
//         mexPrintf("omega[%i] = %6.4f\n",i,omega[i] );       
//         mexPrintf("rho[%i] = %6.4f\n",i,rho[i]);  
        PL_hp[i] = 0;
//         h[i] = h_T[i];  
    }    
     
    
    /* volatility_t_garch to get h_T*/   
    for (i=0; i<N; i++) 
    {
        for (j=1; j<T; j++)
        {   
            h[i] = theta[2*N+i]*h[i] + theta[i] + theta[N+i]*(y[j-1]-theta[3*N+i])*(y[j-1]-theta[3*N+i]);         
//             mexPrintf("y[%i] = %6.4f \n",j,y[j-1]);
//             mexPrintf("h[%i] = %6.4f \n",j,h[i]);
        }                  
     }
    
    /* predict: h, y and pl */
    for (i=0; i<N; i++) 
    {
//         for (j=0; j<4+hp; j++)
//         {
//             mexPrintf("theta[%i] = %6.4f \n",j,theta[N*i+j]);    
//         }
        for (j=1; j<(hp+1); j++)
        {   
            if (j==1)
            {
                h[i] = theta[2*N+i]*h[i] + theta[i] + theta[N+i]*(y[T-1]-theta[3*N+i])*(y[T-1]-theta[3*N+i]); /* forecast based on the last observation */
            }
            else
            {
                h[i] = theta[2*N+i]*h[i] + theta[i] + theta[N+i]*(y_hp[i]-theta[3*N+i])*(y_hp[i]-theta[3*N+i]);  /* forcast based on the previous forecast */                      
            }
/*            y_hp[i] = theta[2*N+i] + sqrt(rho[i]*h[i])*eps[N*j+i]; */
//             mexPrintf("h_hp[%i] = %6.4f \n",j,h[i]);
            rhoh = rho[i]*h[i];
//             mexPrintf("eps[%i]= %6.4f\n",i,theta[N*(j+3)+i]);
            y_hp[i] = theta[3*N+i] + pow(rhoh,0.5)*theta[N*(j+4)+i]; 
            PL_hp[i] = PL_hp[i] + y_hp[i];
        } 
//         mexPrintf("y_hp[%i] = %6.4f \n",i,y_hp[i]);
//         mexPrintf("PL_hp[%i] = %6.4f \n",i,PL_hp[i]);
        PL_hp[i] = 100*(exp(0.01*PL_hp[i])-1);
//         mexPrintf("PL_hp[%i] = %6.4f \n",i,PL_hp[i]);
     }
    
    /* Free allocated memory */
    mxFree(rho);
    mxFree(h);
    mxFree(y_hp);
}

/******************************************************************* */

void duvt_garch(double x, double mu, double sigma, double df, double *GamMat, mwSignedIndex G, double *pdf)
{
    double c0, c1, c2, c, e, tmp, etmp, df5;
    int ind;
    
    c0 = df+1;
//     mexPrintf("c0 = %16.14f\n",c0);    
    
    if (c0 <= 100) /* gamma((nu+1)/2) */
    {
        ind = floor(c0*50000) - 1;
        c1 = GamMat[ind];
    }
    else
    {
     /*   c1 = tgamma(0.5*c0); */
        c1 = GamMat[G-1];
    }
//     mexPrintf("c1 = %16.14f\n",c1);    

    if (df <= 100) /* gamma(nu/2) */
    {
//         ind = floor(df[0]*50000) - 1; /* useround insteadas follows: */
        df5 = df*50000;
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
//     mexPrintf("c2 = %16.14f\n",c2);    

    c = df*PI;
    c = pow(c,0.5);
//     mexPrintf("sqrt(df*pi) = %16.14f\n",c);    
    c2 = c2*c;
//     mexPrintf("c2*sqrt(df*pi) = %16.14f\n",c2);    

    c2 = c2*pow(sigma,0.5);
//     mexPrintf("sqrt(sigma) = %16.14f\n",pow(sigma,0.5));    

    c = c1/c2;
//     mexPrintf("c1/c2 = %16.14f\n",c);    

    e = -0.5*c0; 
//     mexPrintf("e = %16.14f\n",e);    
    
    tmp = (x-mu)*(x-mu)/sigma;
    tmp = 1 + tmp/df;
    etmp = pow(tmp,e); 
    pdf[0] = exp(log(c) + log(etmp));
}

/******************************************************************* */

void prior_t_garch_hl_hyper(double *theta, double *hyper, double *PL_hp, double *VaR,
        mwSignedIndex N, mwSignedIndex hp, mwSignedIndex *r1, double *r2)
{
    mwSignedIndex i, j;
    
    /* Variable size arrays */
    /*  c1 = malloc((N)*sizeof(double));             */
    
    for (i=0; i<N; i++)
    {
        r1[i] = 1;        
        r2[i] = -mxGetInf();

        if (theta[i] <= 0) // omega>0
        {
            r1[i] = 0;
        }   
        if ((theta[i+N] <0 ) || (theta[i+N] >= 1)) // 0<=alpha<1 
        {
            r1[i] = 0;
        }
//         mexPrintf("r1[%i] = %i\n",i,r1[i]);   
        if ((theta[i+2*N] < 0) || (theta[i+2*N] >= 1)) // 0<=beta<1
        {
            r1[i] = 0;
        }
//         mexPrintf("r1[%i] = %i\n",i,r1[i]);   
        if (theta[i+N] + theta[i+2*N] >= 1) //alpha+beta<1
        {
            r1[i] = 0;
        }
//         mexPrintf("r1[%i] = %i\n",i,r1[i]);           

        if (theta[i+4*N] <= 2) //nu>2
        {
            r1[i] = 0;
        }     
//         mexPrintf("r1[%i] = %i\n",i,r1[i]);                  
        if (PL_hp[i] > VaR[0])
        {
            r1[i] = 0;
        }  
//         mexPrintf("r1[%i] = %i\n",i,r1[i]); 
        
        for (j=1; j<(hp+1); j++)
        {
            if (abs(theta[N*(j+4)+i]) > 10) /* for numerical stability */
            {
                r1[i] = 0;
            }
                
        }
        
        if (r1[i] == 1)
        {
            r2[i] = log(hyper[0]) - hyper[0]*(theta[i+4*N] - 2);
        }
    }
}

/*******************************************************************  */

void posterior_t_garch_hl_noS_hyper_mex(double *y, mwSignedIndex N, mwSignedIndex hp, mwSignedIndex T, double *S,
        double *theta, double *hyper, double *VaR, double *GamMat, mwSignedIndex G, double *d,
        double *Td)
{
    mwSignedIndex *r1;
    double *r2; 
    double *PL_hp;
    double *rho;
    double h, *pdf; 
    double rhoh;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));              
    rho = mxMalloc((N)*sizeof(double));
    pdf = mxMalloc((1)*sizeof(double));
    PL_hp = mxMalloc((N)*sizeof(double)); 
    
    Td[0] = (double) T;
    predict_t_garch(theta, y, S,N, hp, T, PL_hp);
    prior_t_garch_hl_hyper(theta, hyper, PL_hp, VaR, N, hp, r1, r2);

    /* Initialise */
    for (i=0; i<N; i++)
    {   
//         mexPrintf("PL_hp[%i] = %6.4f \n",i,PL_hp[i]);   
        rho[i] = (theta[i+4*N]-2)/theta[i+4*N];
//         mexPrintf("omega[%i] = %6.4f\n",i,omega[i] );       
//         mexPrintf("rho[%i] = %6.4f\n",i,rho[i]);  
//         mexPrintf("alpha[%i] = %6.4f\n", i, theta[i]); 
//         mexPrintf("beta[%i] = %6.4f\n", i, theta[i+N]);
//         mexPrintf("mu[%i] = %6.4f\n", i, theta[2*N+i]);
//         mexPrintf("nu[%i] = %6.4f\n", i, theta[i+3*N]);
    }
        
    /* PDF */
    for (i=0; i<N; i++) 
    {
//         mexPrintf("r1[%i] = %i\n",i,r1[i]);       
//         mexPrintf("r2[%i] = %6.4f\n",i,r2[i]);       
       
//         mexPrintf("alpha[%i] = %6.4f\n", i, theta[i]); 
//         mexPrintf("beta[%i] = %6.4f\n", i, theta[i+N]);
//         mexPrintf("mu[%i] = %6.4f\n", i, theta[2*N+i]);
//         mexPrintf("nu[%i] = %6.4f\n", i, theta[i+3*N]);           
//         mexPrintf("eps[%i] = %6.4f\n", i, theta[i+4*N]);           
      
        if (r1[i]==1)
        {
            d[i] = r2[i]; /* prior */
            h = S[0]; 
//             mexPrintf("d[%i] = %6.4f\n",i,d[i]);       
//             mexPrintf("h = %6.4f\n", h);   
            for (j=1; j<T; j++)
            {   
//                 mexPrintf("i = %i\n",i); 
//                 mexPrintf("j = %i\n",j); 
//                 mexPrintf("y[%i] = %6.4f\n", j-1, y[j-1]);  
//                 
//                 mexPrintf("alpha[%i] = %6.4f\n", i, theta[i]); 
//                 mexPrintf("beta[%i] = %6.4f\n", i, theta[i+N]);
//                 mexPrintf("mu[%i] = %6.4f\n", i, theta[2*N+i]);
//                 mexPrintf("nu[%i] = %6.4f\n", i, theta[i+3*N]);
// //     /*            h[i] = beta[i]*h[i] + omega[i] + alpha[i]*(y[j-1]-mu[i])*(y[j-1]-mu[i]);     */    
                h = theta[2*N+i]*h + theta[i] + theta[N+i]*(y[j-1]-theta[3*N+i])*(y[j-1]-theta[3*N+i]);   
//                 mexPrintf("h = %6.4f\n", h);   
                rhoh = rho[i]*h;
//                 mexPrintf("rhoh = %6.4f\n", rhoh);  
                duvt_garch(y[j], theta[i+3*N], rhoh, theta[i+4*N], GamMat, G, pdf);
                pdf[0] = log(pdf[0]);
//                 mexPrintf("pdf[%i] = %16.14f\n", j, pdf[0]);  
                d[i] = d[i] + pdf[0];
//                 mexPrintf("d[%i] = %6.4f\n",i,d[i]);  
            }
//             mexPrintf("d[%i] = %16.14f\n",i,d[i]);  
            for (j=1; j<hp+1; j++) /* pdf of future disturbances */
            {
                duvt_garch(theta[i+N*(j+4)], 0, 1, theta[i+4*N], GamMat, G, pdf);
                pdf[0] = log(pdf[0]);
//                 mexPrintf("eps[%i] = %16.14f\n", j, theta[i+N*(j+3)]); 
//                 mexPrintf("nu[%i] = %16.14f\n", j, theta[i+3*N]); 
//                 mexPrintf("pdf[%i] = %16.14f\n", j, pdf[0]);  
                d[i] = d[i] + pdf[0];
//                 mexPrintf("d[%i] = %6.4f\n",i,d[i]);  
            }
            d[i] = -d[i]/T;
        }
        else
        {
//             d[i] = M;  
            d[i] = mxGetInf();
        }    
//         mexPrintf("d[%i] = %6.4f\n",i,d[i]);  
     }
    
    /* Free allocated memory */
    mxFree(r1); 
    mxFree(r2); 
    mxFree(rho); 
    mxFree(pdf);
    mxFree(PL_hp);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, hp, T, G;                          /* size of matrix */
    double *y, *S, *theta, *hyper, *VaR;                /* input*/
    double *GamMat;                                     /* input */
    double *d;                                          /* output */
    double *Td;
    
    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    S = mxGetPr(prhs[2]);
    VaR = mxGetPr(prhs[3]);
    GamMat = mxGetPr(prhs[4]);
    hyper = mxGetPr(prhs[5]);
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    hp = mxGetN(prhs[0]); /* length of prediction horizon */
    hp = hp - 5;
    T = mxGetM(prhs[1]); /* no. of observations */
    G = mxGetM(prhs[4]);
    
//     mexPrintf("N = %i\n", N);  
//     mexPrintf("hp = %i\n", hp);  
//     mexPrintf("T = %i\n", T); 
//     mexPrintf("G = %i\n", G); 
        
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    Td = mxGetPr(plhs[1]);
    
    /* call the function */
    posterior_t_garch_hl_noS_hyper_mex(y, N, hp, T, S, theta, hyper, VaR, GamMat, G, d, Td);
  
}
