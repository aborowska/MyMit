%% Initialisation
% clc
clear all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath('include/');

x_gam = (0:0.00001:50)' + 0.00001;
GamMat = gamma(x_gam);

T = 10000;
y = randn(T,1);


sigma_init = 0.9;
M = 10000;
N_sim = 20;

MitISEM_Control

% hyper params:
a = 1;
b = 1; 
% logkernel
kernel_init = @(x) - posterior_debug(x, y, a, b);
kernel = @(x) posterior_debug(x, y, a, b);
% kernel([-10:1:10]')
% posterior_debug([-1,0,1,1.5,0.2,-1]', [0,1,2,3,4]', 1, 1)
% mu_init = sigma_init;
[mit1, summary1] = MitISEM(kernel_init, kernel, sigma_init, cont, GamMat);


% for QERMit 2.: draws from the joint
[draw1, lnk1, ind_red1] = fn_rmvgt_robust(M, mit1, kernel);

a_post = a + T/2;
b_post = b + sum(y.^2)/2;
mean_post = b_post/(a_post);
var_post = b_post/(((a_post-1)^2)*(a_post-2));
std_post = sqrt(var_post);