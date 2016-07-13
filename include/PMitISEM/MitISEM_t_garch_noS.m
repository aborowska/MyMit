%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 't_garch2_noS';
algo = 'MitISEM';

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

cont2 = MitISEM_Control;
cont2.mit.dfnc = 5;
cont2.mit.Hmax = 10;
cont2.df.range = [1, 10];

p_bar = 0.01;
H = 10;     % prediction horizon 

VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
accept = zeros(N_sim,1);

VaR_mit = zeros(N_sim,1);
ES_mit = zeros(N_sim,1);
time_mit = zeros(2,1);

%% PRELIM & BIG DRAW
name =  ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name);

% WEIGHTS to initialise MitISEM
% future disturbances are generated from the target thus have weights 1
% log kernel evaluation - only for the parameter draws, epsilons are drawn
% from the target so their weigths are 1
% % kernel = @(xx) posterior_t_garch_noS_mex(xx, data, S, GamMat);
hyper = 0.1; 
kernel = @(a) posterior_t_garch_noS_hyper_mex(a, data, S, GamMat, hyper);
lnk_hl = kernel(draw_hl(:,1:d)); 
% log candidate evaluation
lnd_hl = dmvgt(draw_hl(:,1:d), mit1, true, GamMat);

% importance weights
w_hl = lnk_hl - lnd_hl;
w_hl = exp(w_hl - max(w_hl));

[mu_hl, Sigma_hl] = fn_muSigma(draw_hl, w_hl);
Sigma_hl = reshape(Sigma_hl,d+H,d+H);
Sigma_hl(d+1:d+H,d+1:d+H) = eye(H);
Sigma_hl = reshape(Sigma_hl,1,(d+H)^2);
% cont2.mit.N = 10000;
% cont2.mit.Hmax = 1;

mit_hl.mu = mu_hl;
mit_hl.Sigma = Sigma_hl;
mit_hl.df = cont2.mit.dfnc;
mit_hl.p = 1;
% mu_init = mu_hl;
% mit_init = mit_hl;
cont = cont2;
% cont2.mit.CV_tol = 0.12;

% kernel = @(a) posterior_t_garch_hl_noS_mex(a, data, S, mean(VaR_prelim), GamMat);
kernel_init = @(a) - posterior_t_garch_hl_noS_hyper_mex(a, data, S, mean(VaR_prelim), GamMat, hyper);
kernel = @(a) posterior_t_garch_hl_noS_hyper_mex(a, data, S, mean(VaR_prelim), GamMat, hyper);
tic
% if (H < 40)
%     [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);
% else
    [mit2, summary2] = MitISEM_new(mit_hl, kernel, mu_hl, cont2, GamMat);
% end
time_mit(1,1) = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'mit2','summary2')
end


%% QERMit 2:  MONTE CARLO VaR_mit and ES_mit (and their NSEs) ESTIMATION 
% use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
% to estiamte VaR and ES for theta and y (or alpha in eps)
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
%     kernel = @(a) posterior_t_garch_noS_mex(a, data ,S, GamMat);
    kernel = @(a) posterior_t_garch_noS_hyper_mex(a, data, S, GamMat, hyper);
    lnk = kernel(draw_opt(:,1:d));
    eps_pdf = duvt(draw_opt(:,d+1:H+d), draw_opt(:,d), H, true);
    lnk = lnk + eps_pdf;

    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,1:d), mit1, true, GamMat));
    exp_lnd2 = 0.5*dmvgt(draw_opt, mit2, false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd = log(exp_lnd);

    w_opt = fn_ISwgts(lnk, lnd, false);

    %VaR and ES IS estimates 
    h_T = volatility_t_garch_noS_mex(draw_opt(:,1:d), data, S);
    [y_opt, ~] = predict_t_garch_noS(draw_opt(:,1:d), y_T, h_T, H, draw_opt(:,d+1:d+H));
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_mit(sim,1) = IS_estim(1,1);
    ES_mit(sim,1) = IS_estim(1,2);
    
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_mit(sim,1), model, algo);  
end
time_mit(2,1) = toc/N_sim;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'mit2','summary2','VaR_mit','ES_mit','time_mit')
end

h2 = volatility_t_garch_noS_mex(draw2(:,1:d), data, S);
y2 = predict_t_garch_noS(draw2(:,1:d), y_T, h2, H, draw2(:,d+1:H+d));
PL2 = fn_PL(y2);
mit_eff = sum(PL2 <= mean(VaR_prelim))/(M/2);

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'mit2','summary2','VaR_mit','ES_mit','time_mit','mit_eff')
end

labels_in = {'prelim','mitisem'};
Boxplot_PMitISEM(VaR_prelim, VaR_mit, ES_prelim, ES_mit, model, algo, H, N_sim, true, labels_in);