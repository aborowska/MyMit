
% clc
clear all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath('include/');

x = (0:0.00001:50)' + 0.00001;
GamMat = gamma(x);


data = csvread('GSPC_ret.csv');
data = 100*data;

% QERMit ARCH 
ind_arch = find(data<=-5.5, 1, 'last' );
data = data(1:ind_arch,1);
data = data - mean(data);

T = length(data);
y_T = data(T);
S = var(data);

mu_init = 0.03;
mu_hl = [0.15, -3];
M = 10000;
N_sim = 20;

plot_on = false;
print_on  = false;

MitISEM_Control

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1a.: 
L = true;
kernel_init = @(a) - posterior_arch(a, data, S, L);
kernel = @(a) posterior_arch(a, data, S, L);

VaR_prelim = zeros(N_sim,1);
for sim = 1:N_sim
    [mit1(sim), summary1(sim)] = MitISEM(kernel_init, kernel, mu_init, cont, GamMat);

    [draw1, lnk1, ind_red1] = fn_rmvgt_robust(M, mit1(sim), kernel);
    eps1 = randn(M,1);
    draw1 = [draw1,eps1];


    [alpha, accept] = Mit_MH(M+1000, kernel, mit1(sim), GamMat);
    fprintf('(MitISEM) MH acceptance rate: %4.2f. \n',accept);
    alpha = alpha(1001:M+1000);


    f_stdev = @(aa) sqrt(S+(y_T^2-S)*aa);
    stdev = f_stdev(alpha);
    y_T1 = stdev.*randn(M,1);

    % get the preliminary VaR estimate as the 100th of the ascendingly sorted percentage loss values
    p_bar = 0.01; % p_bar = 1-alpha, 100alpha% VaR
    PL_T1 = sort(fn_PL(y_T1));
    VaR_prelim(sim) = PL_T1(p_bar*M);

    fprintf('(MitISEM) Preliminary VAR estimate: %6.4f. \n',VaR_prelim(sim));
end