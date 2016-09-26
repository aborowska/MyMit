%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
 
x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 't_gas';
algo = 'Prelim';

crisis = true;
recent = false;
old = false;
if crisis 
    y = csvread('GSPC_ret_updated.csv'); 
    results_path = 'results/PMitISEM/crisis/';
elseif recent
    y = csvread('GSPC_ret_updated_short.csv');
    results_path = 'results/PMitISEM/recent';
elseif old
    y = csvread('GSPC_ret_tgarch.csv');
    results_path = 'results/PMitISEM/old/';        
else
    y = csvread('GSPC_ret_updated_short_end.csv');
    results_path = 'results/PMitISEM/';    
end

y = 100*y;
T = size(y,1);
y_T = y(T,1);
 
M = 10000;
BurnIn = 1000;
N_sim = 20;
p_bar = 0.01;
H = 250;

plot_on = false;
save_on = true;

% Control parameters for MitISEM
cont1 = MitISEM_Control;

VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
RNE_prelim = zeros(N_sim,1);
RNE_ES_prelim =  zeros(N_sim,1);
accept = zeros(N_sim,1);
time_prelim = zeros(2,1);

hyper = 0.01;
% kernel_init = @(xx) - posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
kernel_init = @(xx) posterior_t_gas_hyper_init_mex(xx, y, hyper, GamMat);
kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);

% theta = [mu, omega, A, B, nu]
if crisis
    % mu_init = [0, 0.01, 0.1, 0.89, 8];
    mu_init = [0.07, 0.02, 0.1, 0.98, 6.9];    
elseif recent
    mu_init = [0.09, 0.05, 0.14, 0.95, 6.9];    
elseif old
    mu_init = [0, 0.01, 0.1, 0.89, 8];
else
    mu_init = [0.047, 0.0095, 0.06, 0.96, 12];
end
DD = size(mu_init,2);


cont1.mit.dfnc = 5;
cont1.mit.N = 10000;
cont1.resmpl_on = false;
cont1.mit.Hmax = 2;
cont1.mit.iter_max = 2;
tic
[mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont1, GamMat);
time_prelim(1,1) = toc;

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont1','mit1','summary1')
end


tic
for sim = 1:N_sim  
    fprintf('\nPrelim sim = %i.\n', sim);
    kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
    [theta1, accept(sim,1)] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept(sim,1), model, algo);
    theta1 = theta1(BurnIn+1:M+BurnIn,:);

    f_T = volatility_t_gas_mex(theta1, y);
    [y_H, eps_H] = predict_t_gas(theta1, y_T, f_T, H);

    ind_real = (sum(imag(y_H),2)==0);
    M_real = sum(ind_real); 
    fprintf('M_real = %i.\n',M_real)
    y_H = y_H(ind_real,:);
    theta1 = theta1(ind_real,:);  
    eps_H = eps_H(ind_real,:);

    PL_H_ind = fn_PL(y_H);
    PL_H = sort(PL_H_ind);
    VaR_prelim(sim,1) = PL_H(round(p_bar*M_real));
    ES_prelim(sim,1) = mean(PL_H(round(1:p_bar*M_real)));   
    
    ind_prelim = double((PL_H_ind <= VaR_prelim(sim,1)));
    RNE_prelim(sim,1) = fn_RNE(ind_prelim, 'MH',[],'Q');
    ind_prelim = PL_H_ind((PL_H_ind <= VaR_prelim(sim,1)));
    RNE_ES_prelim(sim,1) = fn_RNE(ind_prelim, 'MH',[],'Q');
    
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
end    
time_prelim(2,1) = toc/N_sim;

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','cont1','summary1','accept','time_prelim','RNE_prelim','RNE_ES_prelim')
end


%% RNE
for H = [10,20,40,100,250]
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    % load(name,'RNE_prelim')
    load(name,'mit1','VaR_prelim')
    accept2 = zeros(N_sim,1);
    
    for sim = 1:N_sim
        fprintf('\nPrelim sim = %i.\n', sim);
        kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
        [theta1, accept2(sim,1)] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
        fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept2(sim,1), model, algo);
        theta1 = theta1(BurnIn+1:M+BurnIn,:);

        f_T = volatility_t_gas_mex(theta1, y);
        y_H = predict_t_gas(theta1, y_T, f_T, H);

        ind_real = (sum(imag(y_H),2)==0);
        M_real = sum(ind_real); 
        fprintf('M_real = %i.\n',M_real)
        y_H = y_H(ind_real,:);

        PL_H_ind = fn_PL(y_H);

        ind_prelim = double((PL_H_ind < mean(VaR_prelim)));
        RNE_prelim(sim,1) = fn_RNE(ind_prelim, 'MH',[],'Q');
        ind_prelim = PL_H_ind((PL_H_ind < mean(VaR_prelim)));
        RNE_ES_prelim(sim,1) = fn_RNE(ind_prelim, 'MH',[],'Q');

        fprintf('sum ind_prelim = %i.\n',sum(ind_prelim))    
        fprintf('VaR RNE: %6.4f. \n',RNE_prelim(sim,1));
        fprintf('ES RNE: %6.4f. \n',RNE_ES_prelim(sim,1));            
    end

    if save_on
        name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_RNE_Nsim',num2str(N_sim),'.mat'];
        save(name,'RNE_prelim','RNE_ES_prelim')
    end
end

%% BIG DRAW
% If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
y_predict = @(draw) predict_t_gas_new(draw(:,1:DD), y, H, draw(:,DD+1:end));
tic
[draw_hl, VaR_est, ~, ~] = BigDraw(cont1.mit.N, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat, DD);
time_bigdraw = toc;
 
if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','cont1','summary1','accept','time_prelim','draw_hl','VaR_est','time_bigdraw','RNE_prelim')
end