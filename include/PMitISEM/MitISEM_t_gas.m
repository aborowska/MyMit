%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 't_gas';
algo = 'MitISEM';

crisis = true;
recent = false;
old = false;
if crisis 
    y = csvread('GSPC_ret_updated.csv'); 
    results_path = 'results/PMitISEM/crisis/';
elseif recent
    y = csvread('GSPC_ret_updated_short.csv');
    results_path = 'results/PMitISEM/recent/';
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


% theta = [mu, omega, A, B, nu]
mu_init = [0, 0.01, 0.1, 0.89, 8];
d = size(mu_init,2);

plot_on = true;
save_on = true;

cont2 = MitISEM_Control;

p_bar = 0.01;
H = 20;     % prediction horizon 

VaR_mit = zeros(N_sim,1);
ES_mit = zeros(N_sim,1);
RNE_mit = zeros(N_sim,1);
time_mit = zeros(2,1);

%% PRELIM & BIG DRAW
name =  [results_path,model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name);

N = size(draw_hl,1);

% WEIGHTS to initialise MitISEM
% future disturbances are generated from the target thus have weights 1
% log kernel evaluation - only for the parameter draws, epsilons are drawn
% from the target so their weigths are 1

hyper = 0.01;
kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);

lnk_hl = kernel(draw_hl(:,1:d)); 
% log candidate evaluation
lnd_hl = dmvgt(draw_hl(:,1:d), mit1, true, GamMat);

% importance weights
w_hl = lnk_hl - lnd_hl;
w_hl = exp(w_hl - max(w_hl));

% % restrict the draws to nicely behaving
% [w_hl, ind] = sort(w_hl);
% draw_hl = draw_hl(ind,:);
% draw_hl = draw_hl(0.025*N:0.975*M,:);
% w_hl = w_hl(0.1*N:0.9*M,:);

[mu_hl, Sigma_hl] = fn_muSigma(draw_hl, w_hl);
Sigma_hl = reshape(Sigma_hl,d+H,d+H);
Sigma_hl(d+1:d+H,d+1:d+H) = eye(H);
Sigma_hl = reshape(Sigma_hl,1,(d+H)^2);

mit_hl.mu = mu_hl;
mit_hl.Sigma = Sigma_hl;
mit_hl.df = cont2.mit.dfnc;
mit_hl.p = 1;
% mu_init = mu_hl;
% mit_init = mit_hl;

% cont2.mit.Hmax = 10;
if (H == 10)
    cont2.mit.dfnc = 5;
    cont2.mit.Hmax = 1;  % <<<<<< !!
    cont2.mit.iter_max = 10;  % <<<<<< !! 
%     cont2.mit.CV_tol = 0.01;
elseif (H == 20)
    cont2.mit.dfnc = 10;
    cont2.mit.Hmax = 2;  % <<<<<< !!
    cont2.mit.iter_max = 10;  % <<<<<< !!    
end

cont = cont2;
kernel_init = @(a) - posterior_t_gas_hl_hyper_mex(a, y, hyper, mean(VaR_prelim), GamMat);
kernel = @(a) posterior_t_gas_hl_hyper_mex(a, y, hyper, mean(VaR_prelim), GamMat);

tic
% if (H < 20)
%     [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);
% else
    [mit2, summary2] = MitISEM_new(mit_hl, kernel, mu_hl, cont2, GamMat);
% end
time_mit(1,1) = toc;

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','mit2','summary2')
end


%% QERMit 2:  MONTE CARLO VaR_mit and ES_mit (and their NSEs) ESTIMATION 
% use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
% to estiamte VaR and ES for theta and y (or alpha in eps)

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

tic    
for sim = 1:N_sim
    resampl_on = false;
    fprintf('VaR IS sim = %i.\n', sim);
    draw1 = rmvgt2(M/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
    eps1 = zeros(M/2, H);
    for hh = 1:H
       eps1(:,hh) = trnd(draw1(:,d)); % ERRORS ARE iid T!!
    end
    draw1_eps1 = [draw1(1:M/2,:), eps1];
    draw2 = rmvgt2(M/2, mit2.mu, mit2.Sigma, mit2.df, mit2.p); 
    draw_opt = [draw1_eps1; draw2];

    % IS weights
%     kernel = @(a) posterior_t_garch_noS_mex(a, y ,S, GamMat);
    kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
    lnk = kernel(draw_opt(:,1:d));
    eps_pdf = duvt(draw_opt(:,d+1:H+d), draw_opt(:,d), H, true);
    lnk = lnk + eps_pdf;

    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,1:d), mit1, true, GamMat));
    exp_lnd2 = 0.5*dmvgt(draw_opt, mit2, false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd = log(exp_lnd);

    w_opt = fn_ISwgts(lnk, lnd, false);

    %VaR and ES IS estimates 
    f_T = volatility_t_gas_mex(draw_opt(:,1:d), y);
    y_opt = predict_t_gas(draw_opt(:,1:d), y_T, f_T, H, draw_opt(:,d+1:H+d));
    ind_opt = (fn_PL(y_opt) <= mean(VaR_prelim));
    RNE_mit(sim,1) = fn_RNE(ind_opt, 'IS', w_opt);    
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_mit(sim,1) = IS_estim(1,1);
    ES_mit(sim,1) = IS_estim(1,2);
    
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_mit(sim,1), model, algo);  
end
time_mit(2,1) = toc/N_sim;

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','mit2','summary2','VaR_mit','ES_mit','time_mit','RNE_mit')
end

f2 = volatility_t_gas_mex(draw2(:,1:d), y);
y2 = predict_t_gas(draw2(:,1:d), y_T, f2, H, draw2(:,d+1:H+d));
PL2 = fn_PL(y2);
mit_eff = sum(PL2 <= mean(VaR_prelim))/(M/2);

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','mit2','summary2','VaR_mit','ES_mit','time_mit','mit_eff','RNE_mit')
end

labels_in = {'prelim','mitisem'};
Boxplot_PMitISEM(VaR_prelim, VaR_mit, ES_prelim, ES_mit, model, algo, H, N_sim, true, labels_in);