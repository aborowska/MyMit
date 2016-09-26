%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 't_garch2_noS';
algo = 'Direct';

crisis = false;
recent = false;
old = true;
if crisis 
    y = csvread('GSPC_ret_updated.csv'); 
    results_path = 'results/PMitISEM/crisis/';
elseif recent
    % y = csvread('GSPC_ret_tgarch.csv');
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
y_T = y(T);
S = var(y);
 
M = 10000;
BurnIn = 1000;
N_sim = 20;
p_bar = 0.01;
H = 10;     % prediction horizon 

plot_on = false;
save_on = false;

% Control parameters 
cont_direct = MitISEM_Control;
cont_direct.mit.dfnc = 3;

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
RNE_direct = zeros(N_sim,1);
RNE_ES_direct = zeros(N_sim,1);
accept_direct = zeros(N_sim,1);
time_direct = zeros(2,1);

hyper = 0.01; % 0.01 - VERY UNINFORMATIVE PRIOR
if plot_on
    nu_prior = @(hh,xx) hh*exp(-hh*(xx-2));
    xx = 0.01:0.01:50;
    hold all
    plot(xx,nu_prior(0.01,xx));
    plot(xx,nu_prior(0.1,xx));
    % plot(xx,nu_prior(0.2,xx));
    plot(xx,nu_prior(1,xx));
    hold off
end

% kernel_init = @(a) - posterior_t_garch_noS_mex(a, data, S, GamMat);
% kernel_init = @(a) - posterior_t_garch_noS_hyper_mex(a, y , S, GamMat, hyper);
kernel_init = @(a) posterior_t_garch_noS_hyper_init_mex(a, y , S, GamMat, hyper);
kernel = @(a) posterior_t_garch_noS_hyper_mex(a, y, S, GamMat, hyper);

% theta = [omega, alpha, beta, mu, nu]
if crisis
    mu_init = [0.008, 0.07, 0.9, 0.01, 6.2];
    %     mu_init = [0.02, 0.12, 0.85, 0.075, 6.2];
    tic
    x = fminsearch(kernel_init,mu_init);
    [mu, Sigma] = fn_initopt(kernel_init, x);
    [mu, Sigma] = fn_initopt(kernel_init, mu);
    [mu, Sigma] = fn_initopt(kernel_init, mu);
    mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',cont_direct.mit.dfnc);
    time_direct(1,1) = toc;   
elseif recent
    % mu_init = [0.008, 0.07, 0.9, 0.01, 10];
%     mu_init = [0.009, 0.07, 0.9, 0.05, 11];    
    mu_init = [0.043, 0.17, 0.78, 0.08, 6.1];    
    tic
    [mu, Sigma] = fn_initopt(kernel_init, mu_init);
    mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',cont_direct.mit.dfnc);
    time_direct(1,1) = toc; 
elseif old
% % %     mu_init = [0.009, 0.07, 0.9, 0.05, 11];
% %     mu_init = [0.009, 0.06, 0.9, 0.05, 11];
%     mu_init = [0.006, 0.065, 0.92, 0.048, 10.0];    
    mu_init = [0.008, 0.07, 0.92, 0.048, 10.0];    
    tic
    [mu, Sigma] = fn_initopt(kernel_init, mu_init);
    mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',cont_direct.mit.dfnc);
    time_direct(1,1) = toc;
else
%     mu_init = [0.043, 0.17, 0.78, 0.08, 6.1];    
    mu_init = [0.05, 0.15, 0.78, 0.09, 6.1];    
    tic
    [mu, Sigma] = fn_initopt(kernel_init, mu_init);
    mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',cont_direct.mit.dfnc);
    time_direct(1,1) = toc;     
end


tic
for sim = 1:N_sim
    fprintf('\nDirect sim = %i.\n', sim);
%     kernel = @(a) posterior_t_garch_noS_mex(a, data, S, GamMat);
    kernel = @(a) posterior_t_garch_noS_hyper_mex(a, y , S, GamMat, hyper);
    [theta_direct, accept_direct(sim,1)] = Mit_MH(M+BurnIn, kernel, mit_direct, GamMat);
%     [theta_direct, accept_direct(sim,1), lnw_direct, lnk_direct, lnd_diredct] = Mit_MH_new(M+BurnIn, kernel, mit_direct, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept_direct(sim,1), model, algo);
    theta_direct = theta_direct(BurnIn+1:M+BurnIn,:);

    h_direct = volatility_t_garch_noS_mex(theta_direct, y, S);
    [y_direct, eps_direct] = predict_t_garch_noS(theta_direct, y_T, h_direct, H);

    ind_real = find(sum(imag(y_direct),2)==0);
    M_real = length(ind_real); 
    fprintf('M_real = %i.\n',M_real)
    y_direct = y_direct(ind_real,:);
    theta_direct = theta_direct(ind_real,:);  
    eps_direct = eps_direct(ind_real,:);

    PL_direct_ind = fn_PL(y_direct);
    PL_direct = sort(PL_direct_ind);
    VaR_direct(sim,1) = PL_direct(round(p_bar*M_real));
    ES_direct(sim,1) = mean(PL_direct(round(1:p_bar*M)));   

    ind_direct = double((PL_direct_ind <= VaR_direct(sim,1)));
    RNE_direct(sim,1) = fn_RNE(ind_direct, 'MH',[],'Q');
    ind_direct = PL_direct_ind((PL_direct_ind <= VaR_direct(sim,1)));
    RNE_ES_direct(sim,1) = fn_RNE(ind_direct, 'MH',[],'Q');
    
    fprintf('Direct 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_direct(sim,1), model, algo);
end
time_direct(2,1) = toc/N_sim;

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','mit_direct','accept_direct','time_direct','RNE_direct','RNE_ES_direct')
end

%% RNE
for H = [10,20,40,100,250]
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    % load(name,'RNE_direct')
    load(name,'mit_direct','VaR_direct')
    accept = zeros(N_sim,1);

    for sim = 1:N_sim
        fprintf('\nDirect sim = %i.\n', sim);
        kernel = @(a) posterior_t_garch_noS_hyper_mex(a, y , S, GamMat, hyper);
        [theta_direct, accept(sim,1)] = Mit_MH(M+BurnIn, kernel, mit_direct, GamMat);
        fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept(sim,1), model, algo);
        theta_direct = theta_direct(BurnIn+1:M+BurnIn,:);

        h_direct = volatility_t_garch_noS_mex(theta_direct, y, S);
        y_direct = predict_t_garch_noS(theta_direct, y_T, h_direct, H);

        ind_real = find(sum(imag(y_direct),2)==0);
        M_real = length(ind_real); 
        fprintf('M_real = %i.\n',M_real)
        y_direct = y_direct(ind_real,:);

        PL_direct_ind = fn_PL(y_direct);

        ind_direct = double((PL_direct_ind < mean(VaR_direct)));
        RNE_direct(sim,1) = fn_RNE(ind_direct, 'MH',[],'Q');
        ind_direct = PL_direct_ind((PL_direct_ind < mean(VaR_direct)));
        RNE_ES_direct(sim,1) = fn_RNE(ind_direct, 'MH',[],'Q');    

        fprintf('sum ind_direct = %i.\n',sum(ind_direct))    
        fprintf('VaR RNE: %6.4f. \n',RNE_direct(sim,1));
        fprintf('ES RNE: %6.4f. \n',RNE_ES_direct(sim,1));    
    end

    if save_on
        name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_RNE_Nsim',num2str(N_sim),'.mat'];
        save(name,'RNE_direct','RNE_ES_direct')
    end
end