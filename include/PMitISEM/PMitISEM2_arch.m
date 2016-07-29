%% Initialisation
% clc
clear all
close all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));

x = (0:0.00001:50)' + 0.00001;
GamMat = gamma(x);

algo = 'PMitISEM';
model = 'arch';

% Data
data = csvread('GSPC_ret.csv');
data = 100*data;
ind_arch = find(data<=-5.5, 1, 'last' );
data = data(1:ind_arch,1);
data = data - mean(data);

T = length(data);
y_T = data(T);
S = var(data); % data variance for the variance targeting
         
mu_init = 0.03;
% Metropolis-Hastings for the preliminary
M = 10000; % number of draws for preliminary and IS computations
BurnIn = 1000;

N_sim = 20;
p_bar = 0.01;
H = 10; % forecast horizon
% d = H+1; % dimension of theta
% partition = [1,3:H+1];

plot_on = true;
save_on = false;

% Control parameters for PMitiISEM
cont2 = MitISEM_Control;
cont2.mit.dfnc = 5;

VaR_pmit = zeros(N_sim,1);
ES_pmit = zeros(N_sim,1);
time_pmit = zeros(2,1);

%% PRELIM & BIG DRAW
name =  ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name);

% WEIGHTS to initialise PMitISEM
% future disturbances are generated from the target thus have weights 1
% log kernel evaluation - only for the parameter draws, epsilons are drawn
% from the target so their weigths are 1
kernel = @(xx) posterior_arch(xx, data, S, true);
lnk_hl = kernel(draw_hl(:,1)); 
% log candidate evaluation
lnd_hl = dmvgt(draw_hl(:,1), mit1, true, GamMat);

% importance weights
w_hl = lnk_hl - lnd_hl;
w_hl = exp(w_hl - max(w_hl));

%% PMitISEM
partition = [1,3:H+1];
d = H+1;

fn_const_X = @(a) arch_const_X(a, y_T, S);
fn_input_X = @(xx) xx;
kernel = @(xx) posterior_arch_hl(xx, data, S, mean(VaR_prelim), true);

CV_old = cont2.mit.CV_old;
CV_tol = cont2.mit.CV_tol;

draw0 = draw_hl;
w0 = w_hl;
lnk0 = lnk_hl; %kernel(draw0);
% clear draw_hl w_hl lnk_hl lnd_hl

if (H == 250)
    cont2.mit.Hmax = 1; % <===== !!!
%     cont2.mit.Hmax1 = 5;  
%     cont2.mit.Hmax2 = 2;  
else
    cont2.mit.Hmax = 10; % <===== !!!
end
cont2.mit.iter_max = 2; % <===== !!!
cont2.df.range = [5,20];
cont2.mit.dfnc = 10;

cont = cont2;

% [pmit, CV_mix, CV, iter, pmit_pre, pmit_pre2, pmit_adapt]  = PMitISEM(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont2, GamMat);
% [pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM_debug(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont, GamMat)
tic
[pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM2(draw0, lnk0, w0, kernel, fn_const_X, fn_input_X, partition, d, cont2, GamMat);
time_pmit(1,1) =  toc;

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


s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
pmit = pmit_step3;

tic
for sim = 1:N_sim   
    fprintf('\nVaR IS iter: %d\n',sim)
     
    alpha1 = rmvgt2(M/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
    eps1 = randn(M/2,H); 
    draw1 = [alpha1, eps1];
    
    draw_pmit  = fn_p_rmvgt(M/2, pmit, d, partition, [], fn_const_X);  
    draw_opt = [draw1; draw_pmit];

    kernel = @(xx) posterior_arch(xx, data, S, true);
    lnk_opt = kernel(draw_opt(:,1));

    kernel = @(xx) - 0.5*(log(2*pi) + xx.^2);
    eps_pdf = sum(kernel(draw_opt(:,2:d)),2);
    lnk_opt = lnk_opt + eps_pdf;

    % optimal weights
    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,1), mit1, true, GamMat));
    exp_lnd2 = fn_dpmit(draw_opt, pmit, partition, fn_const_X, true, GamMat);
    exp_lnd2 = 0.5*exp(exp_lnd2);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % IS VaR estimation
    y_opt = predict_arch(draw_opt(:,1), y_T, S, H, draw_opt(:,2:H+1));  
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_pmit(sim,1) = IS_estim(1,1);
    ES_pmit(sim,1) = IS_estim(1,2);   
  
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_pmit(sim,1), model, algo);  
end      
time_pmit(2,1) = toc/N_sim;

VaR_step2 = VaR_pmit;
ES_step2 = ES_pmit;

% VaR_pmit = VaR_step2;
% ES_pmit = ES_step2;

VaR_step2_up = VaR_pmit;
ES_step2_up = ES_pmit;

VaR_step3 = VaR_pmit;
ES_step3 = ES_pmit;


kernel = @(xx) posterior_arch_hl(xx, data, S, mean(VaR_prelim), true);

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit')
end

y_pmit = predict_arch(draw_pmit(:,1), y_T, S, H, draw_pmit(:,2:H+1));  
PL_pmit = fn_PL(y_pmit);
pmit_eff = sum(PL_pmit <= mean(VaR_prelim))/(M/2);

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit','pmit_eff')
end

if plot_on
    [VaR_outlier, ES_outlier] = Boxplot_PMitISEM(VaR_prelim,VaR_pmit,ES_prelim,ES_pmit,model,algo,H,N_sim,true);
    
    y_pmit = predict_arch(draw_pmit(:,1), y_T, S, H, draw_pmit(:,2:H+1));  
    Plot_hor_pmit(y_pmit, y_T, mean(VaR_prelim),model,algo,save_on)
    Beta = Plot_beta(pmit,model,H,save_on,2); % the last parmeter: version==2 ==> plot only the second beta coefficient
end

%% Outliers detection
% V2 = VaR_outlier{2,1};
% VO2 = zeros(length(VaR_outlier{2,1}),1);
% for ii = 1:length(V2)
%     VO2(ii,1) = find(VaR_IS == V2(1,ii));
% end
% 
% E2 = ES_outlier{2,1};
% EO2 = zeros(length(ES_outlier{2,1}),1);
% for ii = 1:length(E2)
%     EO2(ii,1) = find(ES_IS == E2(1,ii));
% end
% 
% ind_redo = [VO2; EO2]; % redo the 'VaR with PMit' loop for these indicies
%  
% %%% REDO and resave
% 
% if save_on
%     name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
%     save(name,'VaR_prelim','ES_prelim','mit1','accept',...
%         'draw_hl','w_hl','lnk_hl','pmit','VaR_IS','ES_IS')
% end