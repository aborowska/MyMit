%% Initialisation
% clc
clear all
close all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));

x = (0:0.00001:50)' + 0.00001;
GamMat = gamma(x);

algo = 'MitISEM';
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
% mu_init = [1 0.03];

% Metropolis-Hastings for the preliminary
M = 10000; % number of draws for preliminary and IS computations
BurnIn = 1000;

N_sim = 20;
p_bar = 0.01;
H = 40; % forecast horizon
% d = H+1; % dimension of theta

plot_on = true;
save_on = false;

% Control parameters for MitISEM (cont) and PMitiISEM (cont2)
cont2 = MitISEM_Control;
cont2.mit.dfnc = 10;
cont2.mit.N = 10000;
cont2.resmpl_on = false;
cont2.df.range = [5,20];


VaR_mit = zeros(N_sim,1);
ES_mit = zeros(N_sim,1);
RNE_mit = zeros(N_sim,1);
time_mit = zeros(2,1);

%% PRELIM & BIG DRAW
name =  ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name);

% WEIGHTS to initialise MitISEM
% future disturbances are generated from the target thus have weights 1
% log kernel evaluation - only for the parameter draws, epsilons are drawn
% from the target so their weigths are 1
kernel = @(xx) posterior_arch(xx, data, S, true);
lnk_hl = kernel(draw_hl(:,1)); 
% log candidate evaluation
lnd_hl = dmvgt(draw_hl(:,1), mit1, true, GamMat);

% importance weights
w_hl = lnk_hl - lnd_hl;
max(w_hl)
w_hl = exp(w_hl - max(w_hl));


%% Standard MitISEM
[mu_hl, Sigma_hl] = fn_muSigma(draw_hl, w_hl);
d = size(draw_hl,2);

mit_hl.mu = mu_hl;
mit_hl.Sigma = Sigma_hl;
mit_hl.df = cont2.mit.dfnc;
mit_hl.p = 1;
% mu_init = mu_hl;
% mit_init = mit_hl;
 
kernel_init = @(a) - posterior_arch_hl(a, data, S, mean(VaR_prelim), true);
kernel = @(a) posterior_arch_hl(a, data, S, mean(VaR_prelim), true);
tic
if (H < 40)
    [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);
else
    [mit2, summary2] = MitISEM_new(mit_hl, kernel, mu_hl, cont2, GamMat);
end
time_mit(1,1) = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','mit2','summary2')
end


%% VaR with standard MitISEM
tic
for sim = 1:N_sim
    fprintf('\nVaR IS iter: %d\n',sim)

    alpha1 = rmvgt2(M/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
    eps1 = randn(M/2,H); 
    draw1 = [alpha1, eps1];

    draw2 = rmvgt2(M/2, mit2.mu, mit2.Sigma, mit2.df, mit2.p); 
    draw_opt = [draw1; draw2];

    %% IS weights
    kernel = @(xx) posterior_arch(xx, data, S, true);
    lnk_opt = kernel(draw_opt(:,1));

    kernel = @(xx) - 0.5*(log(2*pi) + xx.^2);
    eps_pdf = sum(kernel(draw_opt(:,2:d)),2);
    lnk_opt = lnk_opt + eps_pdf;

%     exp_lnd1 = 0.5*normpdf(draw_opt(:,2:H+1)).*dmvgt(draw_opt(:,1), mit1, false, GamMat);
    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,1), mit1, true, GamMat));
    exp_lnd2 = 0.5*dmvgt(draw_opt, mit2, false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % VaR and ES IS estimates
    y_opt = predict_arch(draw_opt(:,1), y_T, S, H, draw_opt(:,2:H+1));  
    ind_opt = (fn_PL(y_opt) <= mean(VaR_prelim));
    RNE_mit(sim,1) = fn_RNE(ind_opt, 'IS', w_opt);
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_mit(sim,1) = IS_estim(1,1);
    ES_mit(sim,1) = IS_estim(1,2);
    PL_opt = fn_PL(y_opt);

    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_mit(sim,1), model, algo);
    fprintf('IS 100*%4.2f%% ES estimate: %6.4f (%s, %s). \n', p_bar, ES_mit(sim,1), model, algo);  
end
time_mit(2,1) = toc/N_sim;

y2 = predict_arch(draw2(:,1), y_T, S, H, draw2(:,2:H+1));  
PL2 = fn_PL(y2);
mit_eff = sum(PL2 <= mean(VaR_prelim))/(M/2);

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','VaR_mit','ES_mit','mit2','summary2','mit_eff','time_mit','RNE_mit')
end


labels_in = {'prelim','mitisem'};
Boxplot_PMitISEM(VaR_prelim, VaR_mit, ES_prelim, ES_mit, model, algo, H, N_sim, true, labels_in);

% Plot_hor_pmit(y_pmit, y_T, mean(VaR_prelim),model,algo,save_on)