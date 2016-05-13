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
data = y;

T = size(data,1);
y_T = data(T);
S = var(data);
 
M = 10000;
BurnIn = 1000;
N_sim = 20;

% hyper = 1;
% theta = [alpha, beta, mu, nu]
% mu_init = [0.03, 0.9, 0.03, 6];
mu_init = [0.065 0.93 0.048 8.4];
d = size(mu_init,2);

model = 't_garch';

plot_on = true;
save_on = true;

MitISEM_Control
N = cont.mit.N;
cont.mit.dfnc = 5;
cont.resmpl_on = false;

cont2 = cont;  
cont2.mit.Hmax = 2;
cont2.df.range = [1, 10];

p_bar = 0.01;
H = 10;     % prediction horizon 

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
accept_direct = zeros(N_sim,1);

VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
accept = zeros(N_sim,1);

VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);


%% Direct
algo = 'Direct';
kernel_init = @(a) - posterior_t_garch_mex(a, data, S, GamMat);

[mu, Sigma] = fn_initopt(kernel_init, mu_init);
cont_direct = cont;
mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',cont_direct.mit.dfnc);

for sim = 1:N_sim
    fprintf('\nPrelim sim = %i.\n', sim);
    kernel = @(a) posterior_t_garch_mex(a, data, S, GamMat);
    [theta_direct, accept_direct(sim,1)] = Mit_MH(M+BurnIn, kernel, mit_direct, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept_direct(sim,1), model, algo);
    theta_direct = theta_direct(BurnIn+1:M+BurnIn,:);

    h_direct = volatility_t_garch_mex(theta_direct, data, S);
    [y_direct, eps_direct] = predict_t_garch(theta_direct, y_T, S, h_direct, H);

    ind_real = (sum(imag(y_direct),2)==0);
    M_real = sum(ind_real); 
    y_direct = y_direct(ind_real,:);
    theta_direct= theta_direct(ind_real,:);  
    eps_direct = eps_direct(ind_real,:);

    [PL_direct, ind] = sort(fn_PL(y_direct));
    VaR_direct(sim,1) = PL_direct(round(p_bar*M_real));
    ES_direct(sim,1) = mean(PL_direct(round(1:p_bar*M)));   
     
    fprintf('Direct 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_direct(sim,1), model, algo);
end

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','mit_direct','accept_direct')
end

%% QERMit 1a.: 
algo = 'MitISEM';
kernel_init = @(a) - posterior_t_garch_mex(a, data , S, GamMat);
kernel = @(a) posterior_t_garch_mex(a, data, S, GamMat);

[mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
if save_on
    save(name,'mit1','summary1')
end

%% QERMit 1b.: 
for sim = 1:N_sim  
    fprintf('\nPrelim sim = %i.\n', sim);
    kernel = @(a) posterior_t_garch_mex(a, data, S, GamMat);
    [theta1, accept(sim,1)] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept(sim,1), model, algo);
    theta1 = theta1(BurnIn+1:M+BurnIn,:);

    h_T = volatility_t_garch_mex(theta1, data, S);
    [y_H, eps_H] = predict_t_garch(theta1, y_T, S, h_T, H);

    ind_real = (sum(imag(y_H),2)==0);
    M_real = sum(ind_real); 
    y_H = y_H(ind_real,:);
    theta1 = theta1(ind_real,:);  
    eps_H = eps_H(ind_real,:);

    [PL_H, ind] = sort(fn_PL(y_H));
    VaR_prelim(sim,1) = PL_H(round(p_bar*M_real));
    ES_prelim(sim,1) = mean(PL_H(round(1:p_bar*M)));   
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
end         

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept','summary1')
end
% If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
kernel = @(xx) posterior_t_garch_mex(xx, data, S, GamMat);
y_predict = @(draw) predict_t_garch_new(draw(:,1:d), data, S, H, draw(:,d+1:d+H));
tic
[draw_hl, VaR_est, ~, ~] = BigDraw(cont2.mit.N, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat, d);
toc

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept','summary1','draw_hl','VaR_est')
end

% log kernel evaluation - only for the parameter draws, epsilons are drawn
% from the target so their weigths are 1
kernel = @(xx) posterior_t_garch_mex(xx, data, S, GamMat);
lnk_hl = kernel(draw_hl(:,1:d)); 
% log candidate evaluation
lnd_hl = dmvgt(draw_hl(:,1:d), mit1, true, GamMat);

% importance weights
w_hl = lnk_hl - lnd_hl;
w_hl = exp(w_hl - max(w_hl));

[mu_hl, Sigma_hl] = fn_muSigma(draw_hl, w_hl);

% cont2.mit.N = 10000;
cont2.mit.Hmax = 10;

mit_hl.mu = mu_hl;
mit_hl.Sigma = Sigma_hl;
mit_hl.df = cont2.mit.dfnc;
mit_hl.p = 1;
% mu_init = mu_hl;
% mit_init = mit_hl;

kernel = @(a) posterior_t_garch_hl_mex(a, data, S, mean(VaR_prelim), GamMat);
[mit2, summary2] = MitISEM_new(mit_hl, kernel, mu_hl, cont2, GamMat);

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept','summary1','draw_hl','VaR_est','mit2','summary2')
end
%% QERMit 2:  MONTE CARLO VaR_IS and ES_IS (and their NSEs) ESTIMATION 
% use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
% to estiamte VaR and ES for theta and y (or alpha in eps)
    
for sim = 1:N_sim
    resampl_on = false;
    fprintf('NSE sim = %i.\n', sim);
    draw1 = rmvgt2(M/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
    eps1 = zeros(M/2, H);
    for hh = 1:H
       eps1(:,hh) = trnd(draw1(:,d)); % ERRORS ARE iid T!!
    end
    draw1_eps1 = [draw1(1:M/2,:), eps1];
    draw2 = rmvgt2(M/2, mit2.mu, mit2.Sigma, mit2.df, mit2.p); 
    draw_opt = [draw1_eps1; draw2];

    % IS weights
    kernel = @(a) posterior_t_garch_mex(a, data ,S, GamMat);
    lnk = kernel(draw_opt(:,1:d));
    eps_pdf = duvt(draw_opt(:,d+1:H+d), draw_opt(:,d), 10, true);
    lnk = lnk + eps_pdf;

    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,1:d), mit1, true, GamMat));
    exp_lnd2 = 0.5*dmvgt(draw_opt, mit2, false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd = log(exp_lnd);

    w_opt = fn_ISwgts(lnk, lnd, false);

    %VaR and ES IS estimates 
    h_T = volatility_t_garch_mex(draw_opt(:,1:d), data, S);
    [y_opt, ~] = predict_t_garch(draw_opt(:,1:d), y_T, S, h_T, H, draw_opt(:,d+1:d+H));
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);
    
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_IS(sim,1), model, algo);  
end


if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept','summary1',...
        'draw_hl','VaR_est','mit2','summary2','VaR_IS','ES_IS')
end


labels_in = {'prelim','mitisem'};
Boxplot_PMitISEM(VaR_prelim, VaR_IS, ES_prelim, ES_IS, model, algo, H, N_sim, true, labels_in);
    