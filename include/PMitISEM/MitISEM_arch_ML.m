%% NON BAYESIAN VERSION 
% PARAMETERS SET TO THEIR ML ESTIMATES

%% Initialisation
% clc
clear all
close all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));

x = (0:0.00001:50)' + 0.00001;
GamMat = gamma(x);

% algo = 'MitISEM';
model = 'arch_ml';
data = csvread('GSPC_ret.csv');
data = 100*data;

% QERMit ARCH 
ind_arch = find(data<=-5.5, 1, 'last' );
data = data(1:ind_arch,1);
data = data - mean(data);

T = length(data);
y_T = data(T);
S = var(data); % data variance for the variance targeting
         
mu_init = 0.03;
% mu_init = [1 0.03];

M = 10000;
N_sim = 100;
p_bar = 0.01;
H = 100; % forecast horizon
% d = H+1; % dimension of theta
% partition = [1,3:H+1];

plot_on = true;
save_on = true;

% Control parameters for MitISEM (cont) and PMitiISEM (cont2)
MitISEM_Control
cont2 = cont;
cont2.mit.iter_max = 1;
cont2.df.range = [1,10];

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);

VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);

%% Direct
algo = 'Direct';
kernel_init = @(a) - posterior_arch(a, data, S, true);
% ML estimate 
[mu, Sigma] = fn_initopt(kernel_init, mu_init);
alpha = repmat(mu,M,1);
 
for sim = 1:N_sim
    fprintf('\nDirect sim = %i.\n', sim);
    eps_direct = randn(M,H);    
    y_direct = predict_arch(alpha, y_T, S, H, eps_direct);
    PL_direct = sort(fn_PL(y_direct));
    VaR_direct(sim,1) = PL_direct(p_bar*M);
    ES_direct(sim,1) = mean(PL_direct(1:p_bar*M));
    fprintf('Direct 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_direct(sim,1), model, algo);
end

%% Choose the starting point (mu_hl) for the constuction of the approximaton
kernel = @(xx) posterior_arch_hl(xx, data, S, true);
y_predict = @(draw) predict_arch(draw(:,1), y_T, S, H, draw(:,2:end));  
[draw_hl, VaR_est, ~, ~] = BigDraw(cont2.mit.N, H, 0, p_bar, mu, kernel, y_predict, GamMat);
draw_hl = draw_hl(:,2:H+1);  
% importance weights --> ARE EQUAL 
w_hl = ones(M,1)/M;

%% Standard MitISEM --> ONLY FOR THE EPSILONS!
algo = 'MitISEM';
[mu_hl, Sigma_hl] = fn_muSigma(draw_hl, w_hl);
d = size(draw_hl,2);
Sigma_hl_mat = reshape(Sigma_hl,d,d);

mit_hl.mu = mu_hl;
mit_hl.Sigma = Sigma_hl;
mit_hl.df = cont2.mit.dfnc;
mit_hl.p = 1;
% mu_init = mu_hl; 
% mit_init = mit_hl;
cont2.mit.Hmax = 10;
cont2.mit.dfnc = 5;
% cont = cont2;

kernel_init = @(xx) - MLtarget_arch_hl(mu, xx, y_T, S, mean(VaR_direct), true);
kernel = @(xx) MLtarget_arch_hl(mu, xx, y_T, S, mean(VaR_direct), true);
if (H < 2)
    [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);
else
    [mit2, summary2] = MitISEM_new(mit_hl, kernel, mu_hl, cont2, GamMat);
end
    
if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','mit2','draw_hl','summary2')
end

%% VaR with standard MitISEM
for sim = 1:N_sim
    fprintf('\nVaR IS iter: %d\n',sim)

    eps1 = randn(M/2,H); 
    eps2 = rmvgt2(M/2, mit2.mu, mit2.Sigma, mit2.df, mit2.p); 
    draw_opt = [eps1; eps2];

    %% IS weights
    kernel = @(xx) - 0.5*(log(2*pi) + xx.^2);
    lnk_opt = sum(kernel(draw_opt),2);

    exp_lnd1 = 0.5*exp(lnk_opt);
    exp_lnd2 = 0.5*dmvgt(draw_opt, mit2, false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % VaR and ES IS estimates
    y_opt = predict_arch(alpha, y_T, S, H, draw_opt);  
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);
    PL_opt = fn_PL(y_opt);

    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_IS(sim,1), model, algo);
    fprintf('IS 100*%4.2f%% ES estimate: %6.4f (%s, %s). \n', p_bar, ES_IS(sim,1), model, algo);  
end
 
if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','mit2','draw_hl','summary2','VaR_IS','ES_IS')
end

labels_in = {'direct','mitisem'};
Boxplot_PMitISEM(VaR_direct, VaR_IS, ES_direct, ES_IS, model, H, N_sim, true, labels_in);
