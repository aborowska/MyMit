%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

y = csvread('GSPC_ret_tgarch.csv');
y = 100*y;
% y = csvread('GSPC_ret.csv');
% y = y - mean(y);
data = y;

T = size(data,1);
y_T = data(T);
S = var(data);

p_bar = 0.01;
H = 2; % forecast horizon

M = 10000;
BurnIn = 1000;
N_sim = 20;

% L = true;
% hyper = 1;
% theta = [alpha, beta, mu, nu]
% mu_init = [0.03, 0.9, 0.03, 6];
mu_init = [0.065 0.93 0.048 8.4];

algo = 'PMitISEM';
model = 't_garch';

plot_on = true;
save_on = true;

% Control parameters for MitISEM (cont) and PMitiISEM (cont2)
MitISEM_Control
N = cont.mit.N;
cont.mit.dfnc = 5;
cont.resmpl_on = false;

cont2 = cont;
cont2.mit.iter_max = 1;
cont2.df.range = [1,20];

VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
accept = zeros(N_sim,1);

VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);

 
%% QERMit 1a.:
kernel_init = @(a) - posterior_t_garch_mex(a, data , S, GamMat);
kernel = @(a) posterior_t_garch_mex(a, data, S, GamMat);

[mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);

for sim = 1:N_sim  
    fprintf('\nPrelim sim = %i.\n', sim);

    %% QERMit 1b.:
    % generate set opf draws of theta using independence MH with
    % candiate from MitISEM; then simulate returns based on the draw of theta 
    [theta1, accept(sim,1)] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept(sim,1), model, algo);
    theta1 = theta1(BurnIn+1:M+BurnIn,:);

    %% High loss, 10 days horizon
    % approximate the high loss distribution of (theta,eps*) where eps*={eps_T+1,...,eps_T+hp}
    h_T = volatility_t_garch_mex(theta1, data, S);
    [y_H, eps1] = predict_t_garch(theta1, y_T, S, h_T, H);

    ind_real = (imag(sum(y_H,2))==0);
    M_real = sum(ind_real); 
    y_H = y_H(ind_real,:);
    theta1 = theta1(ind_real,:);  
    eps1 = eps1(ind_real,:);

    % get the preliminary 10-day-ahead 99% VaR estimate as the 100th of the ascendingly sorted percentage loss values
    [PL, ind] = sort(fn_PL(y_H));

    VaR_prelim(sim,1) = PL(round(p_bar*M_real));
    ES_prelim(sim,1) = mean(PL(round(1:p_bar*M)));   
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
end

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept')
end

if plot_on
    Plot_hor_direct(y_H, y_T, VaR_prelim(sim,1), model, save_on);
end

% If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
kernel = @(xx) posterior_t_garch_mex(xx, data, S, GamMat);
y_predict = @(draw) predict_t_garch_new(draw(:,1:4), data, S, H, draw(:,5:end));
tic
[draw_hl, VaR_est, ~, ~] = BigDraw(cont2.mit.N, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat, 4);
toc

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept','draw_hl','VaR_est')
end

% log kernel evaluation - only for the parameter draws, epsilons are drawn
% from the target so their weigths are 1
kernel = @(xx) posterior_t_garch_mex(xx, data, S, GamMat);
lnk_hl = kernel(draw_hl(:,1:4)); 
% log candidate evaluation
lnd_hl = dmvgt(draw_hl(:,1:4), mit1, true, GamMat);

% importance weights
w_hl = lnk_hl - lnd_hl;
w_hl = exp(w_hl - max(w_hl));

%% PMitISEM
partition = [1,6:H+4];
d = H+4;
% S = var(data);

fn_const_X = @(xx) t_garch_const_X(xx, data, S);
kernel = @(a) posterior_t_garch_hl_mex(a, data, S, mean(VaR_prelim), GamMat);
CV_old = cont.mit.CV_old;
CV_tol = cont.mit.CV_tol;

draw0 = draw_hl;
w0 = w_hl;
lnk0 = lnk_hl; %kernel(draw0);
% clear draw_hl w_hl lnk_hl lnd_hl

cont2.mit.iter_max = 5;
cont2.df.range = [1,10];
cont2.mit.Hmax = 4;
cont = cont2;
% [pmit, CV_mix, CV, iter, pmit_pre, pmit_pre2, pmit_adapt] = PMitISEM(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont2, GamMat);
[pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM_debug(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont, GamMat)

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept','draw_hl','VaR_est','pmit','CV_mix','CV')
end

%% VaR with PMit
for sim = 1:N_sim   
    fprintf('\nVaR IS iter: %d\n',sim)
     
    theta1 = rmvgt2(M/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
    eps1 = zeros(M/2, H);
    for hh = 1:H
        eps1(:,hh) = trnd(theta1(:,4)); % ERRORS ARE iid T!!
    end
    draw1 = [theta1, eps1]; 
    clear theta1 eps1
    
    draw_pmit  = fn_p_rmvgt(M/2, pmit, d, partition, [], fn_const_X);  
    draw_opt = [draw1; draw_pmit];
    clear draw1 %draw_pmit
 
    kernel = @(xx) posterior_t_garch_mex(xx, data ,S, GamMat);
    lnk_opt = kernel(draw_opt(:,1:4)); 
    
%     eps_pdf = zeros(M, 1);
%     for hh = 1:H
%         eps_pdf = eps_pdf + log(tpdf(draw_opt(:,4+hh),draw_opt(:,4)));
%     end
    eps_pdf = duvt(draw_opt(:,4+1:H+4), draw_opt(:,4), H, true);
    lnk_opt = lnk_opt + eps_pdf;
    
    % optimal weights
%     [s1, s2] = fn_partition_ends(partition, d, 1);
    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,1:4), mit1, true, GamMat));
    exp_lnd2 = fn_dpmit(draw_opt, pmit, partition, fn_const_X, true, GamMat);

    exp_lnd2 = 0.5*exp(exp_lnd2);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % IS VaR estimation
    h_T = volatility_t_garch_mex(draw_opt(:,1:4), data, S);
    [y_opt, ~] = predict_t_garch(draw_opt(:,1:4), y_T, S, h_T, H, draw_opt(:,5:d));
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);   
  
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_IS(sim,1), model, algo);  
end


if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept','draw_hl','VaR_est','pmit','CV_mix','CV','VaR_IS','ES_IS')
end

if plot_on
    Boxplot_PMitISEM(VaR_prelim,VaR_IS,ES_prelim,ES_IS,model,algo,H,N_sim,save_on);
    
    h_T = volatility_t_garch_mex(draw_pmit(:,1:4), data, S);
    [y_pmit, ~] = predict_t_garch(draw_pmit(:,1:4), y_T, S, h_T, H, draw_pmit(:,5:d));
    Plot_hor_pmit(y_pmit, y_T, mean(VaR_prelim),model,algo,save_on)
end