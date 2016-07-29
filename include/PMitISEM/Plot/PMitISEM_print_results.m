clear all 
close all
addpath(genpath('include/'));

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
        model_tex = 'GAS(1,1)-$t$';          
        ML = false;
    case 't_gas_ML'
        model_tex = 'GAS(1,1)-$t$';   
        ML = true;
    case 't_garch2_noS'
        model_tex = 'GARCH(1,1)-$t$';
        ML = false;        
    case 'arch'
        model_tex = 'ARCH(1)';
        ML = false;        
    case 'WN'
        model_tex = 'White Noise';
        ML = false;        
    case 'WN_ML'
        model_tex = 'White Noise';    
        ML = true;        
end

horizons = [10,20,40,100,250];

if ML
    algos = {'Direct','MitISEM','PMitISEM'};
else
    algos = {'Direct','Prelim','MitISEM','PMitISEM'};
end
VaR_mat_prelim = zeros(2,5);
VaR_mat_pmit = zeros(2,5);
ES_mat_prelim = zeros(2,5);
ES_mat_pmit = zeros(2,5);

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
estimation = 'mle';
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

%% VaR results PMitISEM all H
fname = ['results/PMitISEM/results_',model,'_pmitisem_comp.tex'];
FID = fopen(fname, 'w+');
fprintf(FID, '\\begin{table}[h] \n');
fprintf(FID, '\\centering \n');

caption = ['\\caption{Results for the $99\\%%$ VaR and ES, in the ',...
    model_tex,' model, based on $N=',int2str(M),'$ candidate draws and $',...
    int2str(N_sim),'$ replications to obtain NSEs.} \n'];
fprintf(FID, caption);

label = ['\\label{tab:res_pmit_',model,'} \n'];
fprintf(FID, label);
fprintf(FID, '\\begin{tabular}{ccccccc}  \n');
fprintf(FID, ' Horizon & & $VaR_{adapt}$ & $VaR_{pmit}$ & & $ES_{adapt}$ & $ES_{pmit}$ \\\\ \\hline \n');

ii = 0;
for h = horizons
    ii = ii+1;
    VaR_prelim = NaN;
    VaR_pmit = NaN;
    ES_prelim = NaN;
    ES_pmit = NaN;
    fprintf(FID, '$%s$ & & ',num2str(h));
    name = ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_prelim','ES_prelim')
    name = ['results/PMitISEM/',model,'_PMitISEM_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    try
        load(name,'VaR_pmit','ES_pmit')
    catch
        
    end
    fprintf(FID, '%6.4f & ' ,mean(VaR_prelim)); VaR_mat_prelim(1,ii) = mean(VaR_prelim);
    fprintf(FID, '%6.4f & & ',mean(VaR_pmit));  VaR_mat_pmit(1,ii) = mean(VaR_pmit);
    fprintf(FID, '%6.4f & ' ,mean(ES_prelim));  ES_mat_prelim(1,ii) = mean(ES_prelim);
    fprintf(FID, '%6.4f  \\\\ \n',mean(ES_pmit));  ES_mat_pmit(1,ii) = mean(ES_pmit);

    fprintf(FID, ' & & ');
    fprintf(FID, '(%6.4f) & ' ,std(VaR_prelim)); VaR_mat_prelim(2,ii) = std(VaR_prelim);
    fprintf(FID, '(%6.4f) & & ',std(VaR_pmit)); VaR_mat_pmit(2,ii) = std(VaR_pmit);
    fprintf(FID, '(%6.4f) & ' ,std(ES_prelim)); ES_mat_prelim(2,ii) = std(ES_prelim);
    fprintf(FID, '(%6.4f)   \\\\ [1ex] \n',std(ES_pmit));  ES_mat_pmit(2,ii) = std(ES_pmit);
    
end 
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular} \n');
fprintf(FID, '\\end{table} \n');
fclose(FID);

%% VaR results diff algo 
Print_PMitISEM_alg_comb(model, horizons, algos, M, N_sim, p_bar)

%% Pmit print
Print_pmit(model,model_tex,10,p_bar,N_sim)

%% Time-Precision Combined
save_on = true;
H = 250;
Plot_time_precision(model, save_on, H, p_bar, N_sim, M)


results_WN = Print_time_precision('WN',horizons,p_bar,N_sim,M);
results_arch = Print_time_precision('arch',horizons,p_bar,N_sim,M);
results_t_garch = Print_time_precision('t_garch2_noS',horizons,p_bar,N_sim,M);
results_t_gas = Print_time_precision('t_gas',horizons,p_bar,N_sim,M);


results_t_gas_ML = Print_time_precision('t_gas_ML',horizons,p_bar,N_sim,M);
results_WN_ML_true = Print_time_precision('WN_ML',horizons,p_bar,N_sim,M,'true');
results_WN_ML_mle = Print_time_precision('WN_ML',horizons,p_bar,N_sim,M,'mle');

%% Horizons
H = 10;
% model = 't_gas';
% model = 't_garch2_noS';
model = 'arch';

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);
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