%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
 
x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

y = csvread('GSPC_ret_tgarch.csv');
y = 100*y;
T = length(y);

% theta = [mu, omega, A, B]
theta = [0, 0.01, 0.1, 0.89];

% TRANSFORMATIONS WHEN NO PRIORS:
% kernel_init = @(xx) - posterior_gas_init(xx, y);
% fn_delta = @(par, hessian) delta_gas(par, hessian, T);

% invlogsig = @(xx) -log((1-xx)./xx);
% mu_init = theta; 
% mu_init(2) = log(mu_init(2));
% mu_init(4) = invlogsig(mu_init(4));
% [mu, Sigma] = fn_initopt(kernel_init, mu_init, fn_delta);
% mu_opt = mu;
% mu_opt(2) = exp(mu_opt(2));
% mu_opt(4) = logsig(mu_opt(4));

mu_init = theta;
kernel_init = @(xx) -posterior_gas(xx, y, true);
kernel = @(xx) posterior_gas(xx, y, true);
[mu, Sigma] = fn_initopt(kernel_init, mu_init);    
try
    r = fn_testSigma(Sigma);
catch
    r = 1;
end
