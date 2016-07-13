clear all
close all
addpath(genpath('include/'));

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

y = csvread('GSPC_ret_tgarch.csv');
y = 100*y;
data = y;

T = size(data,1);
y_T = data(T);
S = var(data);

% theta = [omega, alpha, beta, mu, nu]
mu_init = [0.008, 0.07, 0.9, 0.01, 10];

cont = MitISEM_Control;
cont.mit.dfnc = 5;


% try
%     r = fn_testSigma(Sigma);
% catch
%     r = 1;
% end


%% fixed hyper = 1;
% L = true;
% kernel_init = @(a) - posterior_t_garch_noS(a, data, S, L, hyper, GamMat);
% kernel = @(a) posterior_t_garch_noS(a, data, S, L, hyper, GamMat);
kernel_init = @(a) - posterior_t_garch_noS_mex(a, data , S, GamMat);
kernel = @(a) posterior_t_garch_noS_mex(a, data, S, GamMat);
kernel(mu_init)

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
[mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
mit1

%% mean(nu) = 1
hyper = 1; 
kernel_init = @(a) - posterior_t_garch_noS_hyper_mex(a, data , S, GamMat, hyper);
kernel = @(a) posterior_t_garch_noS_hyper_mex(a, data, S, GamMat, hyper);
kernel(mu_init)

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
[mit1_hyper1, summary1_hyper1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
mit1_hyper1

%% mean(nu) = 10
hyper = 0.1; 
kernel_init = @(a) - posterior_t_garch_noS_hyper_mex(a, data , S, GamMat, hyper);
kernel = @(a) posterior_t_garch_noS_hyper_mex(a, data, S, GamMat, hyper);
mu_init = [0.009, 0.07, 0.9, 0.05, 11];

        [mu, Sigma] = fn_initopt(kernel_init, mu_init);

        try
            r = fn_testSigma(Sigma);
        catch
            r = 1;
        end


kernel(mu_init)
[mit1_hyper01, summary1_hyper01] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
mit1_hyper01

[mu, Sigma] = fn_initopt(kernel_init, mu_init);
r = fn_testSigma(Sigma);