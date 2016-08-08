%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 't_gas';
algo = 'Direct';

y = csvread('GSPC_ret_tgarch.csv');
y = 100*y;
data = y;

T = size(data,1);
y_T = data(T);
S = var(data);
 
M = 10000;
BurnIn = 1000;
N_sim = 20;

plot_on = true;
save_on = true;

p_bar = 0.01;
H = 250;     % prediction horizon 

% Control parameters for MitISEM  
cont_direct = MitISEM_Control;
cont_direct.mit.dfnc = 3;

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
RNE_direct = zeros(N_sim,1);
accept_direct = zeros(N_sim,1);
time_direct = zeros(2,1);

% kernel_init = @(a) - posterior_t_garch_noS_mex(a, data, S, GamMat);

hyper = 0.01;
kernel_init = @(xx) - posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
% mu_init = [0, 0.01, 0.1, 0.89, 8];
mu_init = [0.047, 0.0095, 0.06, 0.96, 12];

%   0.0473    0.0098    0.0666    0.9931   11.8485

tic
[mu, Sigma] = fn_initopt(kernel_init, mu_init);
mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',cont_direct.mit.dfnc);
time_direct(1,1) = toc;

tic
for sim = 1:N_sim
    fprintf('\nDirect sim = %i.\n', sim);
    kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
    [theta_direct, accept_direct(sim,1)] = Mit_MH(M+BurnIn, kernel, mit_direct, GamMat);
%     [theta_direct, accept_direct(sim,1), lnw_direct, lnk_direct, lnd_diredct] = Mit_MH_new(M+BurnIn, kernel, mit_direct, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept_direct(sim,1), model, algo);
    theta_direct = theta_direct(BurnIn+1:M+BurnIn,:);

    f_direct = volatility_t_gas_mex(theta_direct, y);
    [y_direct, eps_direct] = predict_t_gas(theta_direct, y_T, f_direct, H);

    ind_real = find(sum(imag(y_direct),2)==0);
    M_real = length(ind_real); 
    fprintf('M_real = %i.\n',M_real)
    y_direct = y_direct(ind_real,:);
    theta_direct = theta_direct(ind_real,:);  
    eps_direct = eps_direct(ind_real,:);

    PL_direct_ind = fn_PL(y_direct);
    PL_direct = sort(PL_direct_ind);
    VaR_direct(sim,1) = PL_direct(round(p_bar*M_real));
    ES_direct(sim,1) = mean(PL_direct(round(1:p_bar*M_real)));   
     
    ind_direct = double((PL_direct_ind < VaR_direct(sim,1)));
    RNE_direct(sim,1) = fn_RNE(ind_direct, 'MH',[],'Q');
    fprintf('Direct 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_direct(sim,1), model, algo);
end
time_direct(2,1) = toc/N_sim;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','mit_direct','accept_direct','time_direct','RNE_direct')
end