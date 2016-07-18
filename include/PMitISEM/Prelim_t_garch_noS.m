%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 't_garch2_noS';
algo = 'Prelim';

y = csvread('GSPC_ret_tgarch.csv');
y = 100*y;
data = y;

T = size(data,1);
y_T = data(T);
S = var(data);
 
M = 10000;
BurnIn = 1000;
N_sim = 20;


% L = true;
% hyper = 1;
% theta = [omega, alpha, beta, mu, nu]
% mu_init = [0.008, 0.07, 0.9, 0.01, 10];
% % mu_init = [0.065 0.93 0.048 8.4];
mu_init = [0.009, 0.07, 0.9, 0.05, 11];
d = size(mu_init,2);

plot_on = true;
save_on = false;

% Control parameters for MitISEM
cont1 = MitISEM_Control;
cont1.mit.dfnc = 5;
cont1.mit.N = 10000;
cont1.resmpl_on = false;
 
p_bar = 0.01;
H = 10;     % prediction horizon 

VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
accept = zeros(N_sim,1);
time_prelim = zeros(2,1);

% % kernel_init = @(a) - posterior_t_garch_noS(a, data, S, L, hyper, GamMat);
% % kernel = @(a) posterior_t_garch_noS(a, data, S, L, hyper, GamMat);
% kernel_init = @(a) - posterior_t_garch_noS_mex(a, data , S, GamMat);
% kernel = @(a) posterior_t_garch_noS_mex(a, data, S, GamMat);
hyper = 0.1; 
kernel_init = @(a) - posterior_t_garch_noS_hyper_mex(a, data , S, GamMat, hyper);
kernel = @(a) posterior_t_garch_noS_hyper_mex(a, data, S, GamMat, hyper);


tic
[mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont1, GamMat);
time_prelim(1,1) = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'mit1','summary1')
end

%% QERMit 1b.: 
tic
for sim = 1:N_sim  
    fprintf('\nPrelim sim = %i.\n', sim);
%     kernel = @(a) posterior_t_garch_noS_mex(a, data, S, GamMat);
    kernel = @(a) posterior_t_garch_noS_hyper_mex(a, data, S, GamMat, hyper);
    [theta1, accept(sim,1)] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept(sim,1), model, algo);
    theta1 = theta1(BurnIn+1:M+BurnIn,:);

    h_T = volatility_t_garch_noS_mex(theta1, data, S);
    [y_H, eps_H] = predict_t_garch_noS(theta1, y_T, h_T, H);

    ind_real = (sum(imag(y_H),2)==0);
    M_real = sum(ind_real); 
    fprintf('M_real = %i.\n',M_real)
    y_H = y_H(ind_real,:);
    theta1 = theta1(ind_real,:);  
    eps_H = eps_H(ind_real,:);

    [PL_H, ind] = sort(fn_PL(y_H));
    VaR_prelim(sim,1) = PL_H(round(p_bar*M_real));
    ES_prelim(sim,1) = mean(PL_H(round(1:p_bar*M)));   
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
end    
time_prelim(2,1) = toc/N_sim;


if plot_on
    Plot_hor_direct(y_H,y_T,VaR_prelim(sim,1),model,true);
end

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','cont1','summary1','accept','time_prelim')
end

% If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
% kernel = @(a) posterior_t_garch_noS(a, data, S, L, hyper, GamMat);
% kernel = @(xx) posterior_t_garch_noS_mex(xx, data, S, GamMat);
kernel = @(a) posterior_t_garch_noS_hyper_mex(a, data, S, GamMat, hyper);
% y_predict = @(draw) predict_t_garch_new_noS(draw(:,1:d), data, S, H, draw(:,d+1:end));
y_predict = @(draw) predict_t_garch_new_noS(draw(:,1:d), data, S, draw(:,d+1:end));

tic
% [draw_hl, VaR_est, ~, ~] = BigDraw(cont1.mit.N, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat, d);
[draw_hl, VaR_est, ~, ~] = BigDraw(cont1.mit.N/10, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat, d);
time_bigdraw = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','cont1','summary1','accept','time_prelim','draw_hl','VaR_est','time_bigdraw')
end