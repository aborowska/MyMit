%% Initialisation
% clc
clear all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath('include/');

x = (0:0.00001:50)' + 0.00001;
GamMat = gamma(x);

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
kernel_init = @(x) ((a + 1 + T/2)*log(x) + (b + 0.5*sum(y.^2))./x);
kernel = @(x) (-(a + 1 + T/2)*log(x) - (b + 0.5*sum(y.^2))./x);

[mit1, summary1] = MitISEM(kernel_init, kernel, sigma_init, cont, GamMat);

% for QERMit 2.: draws from the joint
[draw1, lnk1, ind_red1] = fn_rmvgt_robust(M, mit1, kernel);
draw1 = [draw1,randn(M,1)];
