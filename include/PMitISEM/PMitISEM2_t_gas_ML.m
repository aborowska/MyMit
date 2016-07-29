%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
 
x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 't_gas_ML';
algo = 'PMitISEM';

y = csvread('GSPC_ret_tgarch.csv');
y = 100*y;
T = size(y,1);
y_T = y(T,1);

p_bar = 0.01;
H = 250;

M = 10000;
N_sim = 20;
sim = 1;

plot_on = true;
save_on = true;

% Control parameters for PMitISEM
cont2 = MitISEM_Control;
cont2.mit.dfnc = 5;

VaR_pmit = zeros(N_sim,1);
ES_pmit = zeros(N_sim,1);
time_pmit = zeros(2,1);

%% PRELIM & BIG DRAW
name =  ['results/PMitISEM/',model,'_Direct_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name);

theta_mat = repmat(theta_mle,M,1);
f_mat = repmat(f_mle,M,1);    
nu_mat = theta_mat(:,5);

draw_hl = draw_hl(:,6:end);
[N,d] = size(draw_hl);
w_hl = ones(N,1);
lnk_hl = sum(duvt(draw_hl, theta_mle(1,5)*ones(N,1), H, true),2);

%% PMitISEM
partition = 1:H;
d = H;

% fn_const_X = @(xx) t_gas_hyper_const_X2(xx, y);
% fn_input_X = @(xx) t_gas_hyper_input_X(xx, y);
fn_const_X = @(xx) t_gas_ML_hyper_const_X2(xx, y_T, theta_mle, f_mle);
fn_input_X = @(xx) xx;

% kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
% kernel = @(a) posterior_t_gas_hl_hyper_mex(a, y, hyper, mean(VaR_prelim), GamMat);
kernel_init = @(xx) - MLtarget_t_gas_hl(xx, theta_mle, f_mle, y_T, mean(VaR_direct));
kernel = @(xx) MLtarget_t_gas_hl(xx, theta_mle, f_mle, y_T, mean(VaR_direct));


CV_old = cont2.mit.CV_old;
CV_tol = cont2.mit.CV_tol;

draw0 = draw_hl;
w0 = w_hl;
lnk0 = lnk_hl; %kernel(draw0);

% if (H > 10)
%     cont2.mit.iter_max = 1;
% else
    cont2.mit.iter_max = 1;%3;%6;%8;
% end
cont2.df.range = [1,10];
cont2.mit.Hmax = 1;
if (H >= 40)
    cont2.mit.dfnc = 10;
end

if (H == 250)
    cont2.mit.Hmax = 1;
end
cont = cont2;

tic
% [pmit, CV_mix, CV, iter, pmit_pre, pmit_pre2, pmit_adapt]  = PMitISEM(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont2, GamMat);
[pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM2(draw0, lnk0, w0, kernel, fn_const_X, fn_input_X, partition, d, cont2, GamMat);
time_pmit(1,1) = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'pmit','CV_mix','CV','iter')
end



%% VaR with PMit

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
pmit = pmit_step2;

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
pmit = pmit_step2_up;

tic
for sim = 1:N_sim  
    fprintf('\nVaR IS iter: %d\n',sim)

    draw_pmit = fn_p_rmvgt2(M, pmit, d, partition, [], fn_const_X, fn_input_X);             
    [lnk_opt, PL_mit] = MLtarget_t_gas_hl(draw_pmit, theta_mle, f_mle, y_T, Inf);
    lnd_opt = fn_dpmit2(draw_pmit, pmit, partition, fn_const_X, true, GamMat);        
    w_opt = exp(lnk_opt - lnd_opt)/M;
    [PL, ind] = sort(PL_mit);         
    w_opt = w_opt(ind,:);
    cum_w = cumsum(w_opt);
    ind_var = min(find(cum_w >= p_bar))-1; 
    VaR_pmit(sim,1) = PL(ind_var);
    ES = (w_opt(1:ind_var)/sum(w_opt(1:ind_var))).*PL(1:ind_var);
    ES_pmit(sim,1) = sum(ES(isfinite(ES)));    
    
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_pmit(sim,1), model, algo); 
end
time_pmit(2,1) = toc/N_sim;

VaR_step2 = VaR_pmit;
ES_step2 = ES_pmit;

% VaR_pmit = VaR_step2;
% ES_pmit = ES_step2;

VaR_step2_up = VaR_pmit;
ES_step2_up = ES_pmit;

% VaR_pmit = VaR_step2_up;
% ES_pmit = ES_step2_up;

% time_pmit(1,1) = time_pmit(1,1) + time_step2_up;


if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit')
end


y_pmit = predict_t_gas(theta_mat, y_T, f_mat, H, draw_pmit);
PL_pmit = fn_PL(y_pmit);
pmit_eff = sum(PL_pmit <= mean(VaR_direct))/M;


if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit','pmit_eff')
end

if plot_on
    labels_in = {'naive','pmit'};
    [VaR_outlier, ES_outlier] = Boxplot_PMitISEM(VaR_direct, VaR_pmit, ES_direct, ES_pmit, model, algo, H, N_sim, true, labels_in);

    [y_pmit, ~] = predict_t_gas(theta_mat, y_T, f_mat, H, draw_pmit);
    Plot_hor_pmit(y_pmit, y_T, mean(VaR_direct), model, algo, save_on)

    [Beta, Sigma, nu] = Plot_beta(pmit,model,H,save_on,2); % the last parmeter: version==2 ==> plot only the second beta coefficient
end