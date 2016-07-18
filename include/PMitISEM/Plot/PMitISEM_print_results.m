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

switch model
    case 't_gas'
        model_tex = 'GAS(1,1)-$t$';        
    case 't_garch2_noS'
        model_tex = 'GARCH(1,1)-$t$';
    case 'arch'
        model_tex = 'ARCH(1,1)';
    case 'WN'
        model_tex = 'White Noise';
end
horizons = [10,20,40,100,250] ;
algos = {'Direct','Prelim','MitISEM','PMitISEM'};

VaR_mat_prelim = zeros(2,5);
VaR_mat_pmit = zeros(2,5);
ES_mat_prelim = zeros(2,5);
ES_mat_pmit = zeros(2,5);

%% Boxplot Combine
for h = horizons
    Boxplot_Combine(model, h, N_sim, p_bar, save_on)
end
 Boxplot_Combine(model, H, N_sim, p_bar, save_on)

for h = horizons
    name = ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_prelim','ES_prelim')
    name = ['results/PMitISEM/',model,'_PMitISEM_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_pmit','ES_pmit')
    Boxplot_PMitISEM(VaR_prelim,VaR_pmit,ES_prelim,ES_pmit,model,'PMitISEM',h,N_sim,true);
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

fname = ['results/PMitISEM/results_',model,'_alg_comp.tex'];
FID = fopen(fname, 'w+');
fprintf(FID, '\\begin{table}[h] \n');
fprintf(FID, '\\centering \n');

caption = ['\\caption{Results for the $99\\%%$ VaR and ES, in the ',...
    model_tex,' model, based on $N=',int2str(M),'$ candidate draws and $',...
    int2str(N_sim),'$ replications to obtain NSEs.} \n'];
fprintf(FID, caption);

label = ['\\label{tab:res_algos_',model,'} \n'];
fprintf(FID, label);
fprintf(FID, '\\begin{tabular}{ccccccccccc}  \n');
fprintf(FID, ' H & & $VaR_{naive}$ & $VaR_{adapt}$ & $VaR_{mit}$  & $VaR_{pmit}$ &  & $ES_{naive}$ & $ES_{adapt}$ & $ES_{mit}$ & $ES_{pmit}$ \\\\ \\hline \n');

for h = horizons
    VaR_direct = NaN;
    VaR_prelim = NaN;
    VaR_mit = NaN;
    VaR_pmit = NaN;
    ES_direct = NaN;
    ES_prelim = NaN;
    ES_mit = NaN;
    ES_pmit = NaN;
    
    fprintf(FID, '%s & & ',num2str(h));
    for algo = algos
        try
            name = ['results/PMitISEM/',model,'_',char(algo),'_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
            load(name)
        catch

        end
    end
    
    fprintf(FID, '%6.4f & ' , mean(VaR_direct));
    fprintf(FID, '%6.4f & ', mean(VaR_prelim));
    fprintf(FID, '%6.4f & ' , mean(VaR_mit));
    fprintf(FID, '%6.4f & & ', mean(VaR_pmit));
    fprintf(FID, '%6.4f & ' , mean(ES_direct));
    fprintf(FID, '%6.4f & ' , mean(ES_prelim));
    fprintf(FID, '%6.4f & ' , mean(ES_mit));
    fprintf(FID, '%6.4f  \\\\ \n', mean(ES_pmit));

    fprintf(FID, ' & & ');
    fprintf(FID, '(%6.4f) & ' , std(VaR_direct));
    fprintf(FID, '(%6.4f) & ' , std(VaR_prelim));
    fprintf(FID, '(%6.4f) & ' , std(VaR_mit));
    fprintf(FID, '(%6.4f) & & ', std(VaR_pmit));
    fprintf(FID, '(%6.4f) & ' , std(ES_direct));
    fprintf(FID, '(%6.4f) & ' , std(ES_prelim));
    fprintf(FID, '(%6.4f) & ' , std(ES_mit));
%     fprintf(FID, '(%6.4f)   \\\\ [1ex] \n', std(ES_pmit));
    fprintf(FID, '(%6.4f)   \\\\ \n', std(ES_pmit));
    
    fprintf(FID, ' & & ');
    fprintf(FID, '$[$%6.4f$]$ & ' , iqr(VaR_direct));
    fprintf(FID, '$[$%6.4f$]$ & ' , iqr(VaR_prelim));
    fprintf(FID, '$[$%6.4f$]$ & ' , iqr(VaR_mit));
    fprintf(FID, '$[$%6.4f$]$ & & ', iqr(VaR_pmit));
    fprintf(FID, '$[$%6.4f$]$ & ' , iqr(ES_direct));
    fprintf(FID, '$[$%6.4f$]$ & ' , iqr(ES_prelim));
    fprintf(FID, '$[$%6.4f$]$ & ' , iqr(ES_mit));
    fprintf(FID, '$[$%6.4f$]$  \\\\ [1ex] \n', iqr(ES_pmit));
end 
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular} \n');
fprintf(FID, '\\raggedright \n\n'); 
fprintf(FID, '\\vspace{5pt}\\footnotesize{NaN: it was not possible to generate the particular result with the corresponding algorithm.} \n');
fprintf(FID, '\\end{table} \n');
fclose(FID);

%% Pmit print
fname = ['results/PMitISEM/results_',model,'_pmit.tex'];
if strcmp(model,'WN')
    param = '\sigma^{2}';
elseif strcmp(model,'arch')
    param = '\alpha';
else %'t_garch','t_gas'
    param = '\theta';
end
h='10';
name = ['results/PMitISEM/',model,'_PMitISEM_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name)
FID = fopen(fname, 'w+');

fprintf(FID, '\\begin{table}[h] \n');
fprintf(FID, '\\centering \n');

caption = ['\\caption{Partial mixture properties for $H=10$ in the ', model_tex,' model.} \n'];
fprintf(FID, caption);

label = ['\\label{tab:pmits_',model,'} \n'];
fprintf(FID, label);
fprintf(FID, '\\begin{tabular}{cccc}  \n');
fprintf(FID, ' Subset & Parameters& No. of components $h_{s}$ & weighted $\\mu$ or $\\beta$  \\\\ \\hline \n');

for ii = 1:10
    fprintf(FID, '%d & ' ,ii);
    if (ii == 1)
        component = ['$\{(',param,',\varepsilon_{1})\}$'];
    else
        component = ['$\{\varepsilon_{',num2str(ii),'}\}$'];
    end
	fprintf(FID, '%s & ' ,component);
    fprintf(FID, '%d & ' ,length( pmit(ii).p));
    g = sprintf('%6.4f, ', pmit(ii).p*pmit(ii).mu);
    fprintf(FID, '[%s]   \\\\ [1ex] \n',g);
end

fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular} \n');
fprintf(FID, '\\end{table} \n');
fclose(FID);

%% Time-Precision
save_on = true;
H = 20;
Plot_time_precision(model, save_on, H, p_bar, N_sim, M)

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
        fn_vol = @(xx) randn(M,H); 
        fn_predict = @(xx,aa,zz) bsxfun(@times,vv,sqrt(zz)); 
        DD = 1;
end 

partition = [1,DD+2:H+DD];
Plot_hor(y, model, DD, H, p_bar, N_sim, M, BurnIn, save_on, kernel, fn_vol, fn_predict, fn_predict2, partition, fn_const_X, fn_input_X, GamMat)