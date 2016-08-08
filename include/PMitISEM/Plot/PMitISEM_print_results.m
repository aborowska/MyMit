clear all 
close all
addpath(genpath('include/'));

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

save_on = 1;
H = 250;
M = 10000; BurnIn = 1000;
N_sim = 20;
p_bar = 0.01;

% model = 't_gas';
% model = 't_garch2_noS';
% model = 'arch';
% model = 'WN';
% model = 't_gas_ML';
% model = 'WN_ML';

switch model
    case 't_gas'
        model_tex = '\\textbf{GAS(1,1)-$t$}';          
        ML = false;
    case 't_gas_ML'
        model_tex = '\\textbf{GAS(1,1)-$t$}';   
        ML = true;
    case 't_garch2_noS'
        model_tex = '\\textbf{GARCH(1,1)-$t$}';
        ML = false;        
    case 'arch'
        model_tex = '\\textbf{ARCH(1)}';
        ML = false;        
    case 'WN'
        model_tex = '\\textbf{White Noise}';
        ML = false;        
    case 'WN_ML'
        model_tex = '\\textbf{White Noise}';    
        ML = true;        
end

horizons = [10,20,40,100,250];

if ML
    algos = {'Direct','MitISEM','PMitISEM'};
else
    algos = {'Direct','Prelim','MitISEM','PMitISEM'};
end

%% Boxplot Combine
for h = horizons
    Boxplot_Combine(model, h, N_sim, p_bar, save_on)
end
Boxplot_Combine(model,100, N_sim, p_bar, save_on)

%% Boxplot PMitISEM and Prelim
% Bayesian case
for h = horizons
    name = ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_prelim','ES_prelim')
    name = ['results/PMitISEM/',model,'_PMitISEM_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_pmit','ES_pmit')
    Boxplot_PMitISEM(VaR_prelim,VaR_pmit,ES_prelim,ES_pmit,model,'PMitISEM',h,N_sim,true);
end

% ML case
estimation = '';
estimation = 'mle';

algo = 'PMitISEM';
for h = horizons
    name = ['results/PMitISEM/',model,'_Direct_',estimation,'_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_direct','ES_direct')
    name = ['results/PMitISEM/',model,'_PMitISEM_',estimation,'_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_pmit','ES_pmit')
    labels_in = {'naive','pmit'};
    if strcmp(model,'WN_ML')
        algo_est = [algo,'_',estimation];
    else
        algo_est = algo;
    end
    Boxplot_PMitISEM(VaR_direct,VaR_pmit,ES_direct,ES_pmit,model,algo_est,h,N_sim,save_on, labels_in);
end

%% VaR results diff algo 
Print_PMitISEM_alg_comb(model, horizons, algos, M, N_sim, p_bar)

%% Pmit print
Print_pmit(model,model_tex,10,p_bar,N_sim)


%% Pmit posterior
model = 't_gas'; 
model = 't_garch2_noS'; 
model = 'arch'; 
switch model
    case 't_gas'
        parameter = {'$\mu$','$\omega$','$A$','$B$','$\nu$'};
        y = csvread('GSPC_ret_tgarch.csv');
        y = 100*y; 
        hyper = 0.1; 
        kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat); 
    case 't_garch2_noS'   
        parameter = {'$\omega$','$\alpha$','$\beta$','$\mu$','$\nu$'};
        y = csvread('GSPC_ret_tgarch.csv');
        y = 100*y;
        S = var(y);  
        hyper = 0.1; 
        kernel = @(xx) posterior_t_garch_noS_hyper_mex(xx, y, S, GamMat, hyper);  
    case 'arch'
        parameter = {'$\alpha$'};
        y = csvread('GSPC_ret.csv');
        y = 100*y;
        ind_arch = find(y<=-5.5, 1, 'last' );
        y = y(1:ind_arch,1);
        y = y - mean(y);
        S = var(y); 
        kernel = @(xx) posterior_arch(xx, y, S, true);      
    case 'WN'     
        parameter = {'$\sigma^{2}$'};
        s = RandStream('mt19937ar','Seed',0);
        RandStream.setGlobalStream(s); 
        T = 10000;
        y = randn(T,1); 
        y = y - mean(y); 
        a = 1;
        b = 1;
        kernel = @(x) posterior_WN(x, y, a, b, true);         
end  
results_arch = Print_posterior([], model, parameter, kernel, save_on, GamMat);
results_t_garch = Print_posterior([], model, parameter, kernel, save_on, GamMat);
results_t_gas = Print_posterior([], model, parameter, kernel, save_on, GamMat);
 

%% Time-Precision Combined
save_on = true;
H = 250;
for h = horizons
    Plot_time_precision(model, save_on, h, p_bar, N_sim, M)
end

results_WN = Print_time_precision('WN',horizons,p_bar,N_sim,M);
results_arch = Print_time_precision('arch',horizons,p_bar,N_sim,M);
Plot_time_precision2(results_arch, 'arch', save_on, horizons) %, p_bar, N_sim, M)

results_t_garch = Print_time_precision('t_garch2_noS',horizons,p_bar,N_sim,M);
Plot_time_precision2(results_t_garch,'t_garch2_noS', save_on, horizons) %, p_bar, N_sim, M)

results_t_gas = Print_time_precision('t_gas',horizons,p_bar,N_sim,M);
Plot_time_precision2(results_t_gas,'t_gas', save_on, horizons) %, p_bar, N_sim, M)


results_t_gas_ML = Print_time_precision('t_gas_ML',horizons,p_bar,N_sim,M);
Plot_time_precision2(results_t_gas_ML,'t_gas_ML', save_on, horizons) %, p_bar, N_sim, M)



results_WN_ML_true = Print_time_precision('WN_ML',horizons,p_bar,N_sim,M,'true');
Plot_time_precision2(results_WN_ML_true,'WN_ML_true', save_on, horizons,'true') %, p_bar, N_sim, M)

results_WN_ML_mle = Print_time_precision('WN_ML',horizons,p_bar,N_sim,M,'mle');
Plot_time_precision2(results_t_gas_ML,'WN_ML', save_on, horizons,'mle') %, p_bar, N_sim, M)


% Efficiency

 Print_eff(model, horizons)
 
 
%% Horizons
H = 10;
% model = 't_gas';
% model = 't_garch2_noS';
model = 'arch';
switch model
    case 't_gas'
        y = csvread('GSPC_ret_tgarch.csv');
        y = 100*y;
        S = var(y); 
        y_T = y(end);
        hyper = 0.1;
        fn_vol = @(xx) volatility_t_gas_mex(xx, y);
        fn_predict = @(xx, vv) predict_t_gas(xx, y_T, vv, H);
        fn_predict2 = @(xx, vv, zz) predict_t_gas(xx, y_T, vv, H, zz);

        fn_const_X = @(xx) t_gas_hyper_const_X2(xx, y);
        fn_input_X = @(xx) t_gas_hyper_input_X(xx, y);
        kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
        DD = 5;
    case 't_garch2_noS'         
        y = csvread('GSPC_ret_tgarch.csv');
        y = 100*y;
        S = var(y); 
        y_T = y(end);
        hyper = 0.1; 
        fn_vol = @(xx) volatility_t_garch_noS_mex(xx, y, S);
        fn_predict = @(xx,vv) predict_t_garch_noS(xx, y_T, S, vv, H);
        fn_predict2 = @(xx,vv,zz) predict_t_garch_noS(xx, y_T, S, vv, H, zz);
        fn_const_X = @(xx) t_garch_noS_const_X2(xx, y, S);
        fn_input_X = @(xx) t_garch_noS_input_X(xx, y, S);          
        kernel = @(xx) posterior_t_garch_noS_hyper_mex(xx, y, S, GamMat, hyper);
        DD = 5;
    case 'arch'
        y = csvread('GSPC_ret.csv');
        y = 100*y;
        ind_arch = find(y<=-5.5, 1, 'last' );
        y = y(1:ind_arch,1);
        y = y - mean(y);
        S = var(y);
        y_T = y(end);
        fn_vol = @(xx) randn(M,H); 
        fn_predict = @(xx,vv) predict_arch(xx, y_T, S, H, vv);
        fn_predict2 = @(xx,aa,zz) predict_arch(xx, y_T, S, H, zz);

        fn_const_X = @(xx) arch_const_X(xx, y_T, S);
        fn_input_X = @(xx) xx;
        kernel = @(xx) posterior_arch(xx, y, S, true);
        DD = 1;
    case 'WN'     
        s = RandStream('mt19937ar','Seed',0);
        RandStream.setGlobalStream(s); 
        T = 10000;
        y = randn(T,1); 
        y = y - mean(y);
        fn_vol = @(xx) randn(M,H); 
        fn_predict = @(xx,vv) bsxfun(@times,vv,sqrt(xx)); 
        fn_predict2 = @(xx,aa,vv) bsxfun(@times,vv,sqrt(xx)); 
       
        fn_const_X = @(xx) WN_const_X(xx);
        fn_input_X = @(xx) xx;
        a = 1;
        b = 1;
        kernel = @(x) posterior_WN(x, y, a, b, true); 
        DD = 1;
end 
save_on = true;
for H = horizons
    fn_vol = @(xx) randn(M,H);
    partition = [1,DD+2:H+DD];
    Plot_hor(y, model, DD, H, p_bar, N_sim, M, BurnIn, save_on, kernel, fn_vol, fn_predict, fn_predict2, partition, fn_const_X, fn_input_X, GamMat)
end