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
model = 'arch';
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

% Metropolis-Hastings for the preliminary
M = 10000; % number of draws for preliminary and IS computations
BurnIn = 1000;

N_sim = 20;
p_bar = 0.01;
H = 10; % forecast horizon
% d = H+1; % dimension of theta
% partition = [1,3:H+1];

plot_on = true;
save_on = true;

% Control parameters for MitISEM (cont) and PMitiISEM (cont2)
MitISEM_Control
cont.mit.dfnc = 5;
cont.mit.N = 10000;
cont.resmpl_on = false;

cont2 = cont;
cont2.mit.iter_max = 1;
cont2.df.range = [1,10];


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
kernel_init = @(a) - posterior_arch(a, data, S, true);
% kernel_init = @(a) - posterior_arch_noS(a, data, true);

[mu, Sigma] = fn_initopt(kernel_init, mu_init);
cont_direct = cont;
mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',cont_direct.mit.dfnc);


for sim = 1:N_sim
    fprintf('\nDirect sim = %i.\n', sim);
    kernel = @(a) posterior_arch(a, data, S, true);
    % kernel = @(a) posterior_arch_noS(a, data, true);
    [alpha_direct, accept_direct(sim,1)] = Mit_MH(M+BurnIn, kernel, mit_direct, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept_direct(sim,1), model, algo);
    alpha_direct = alpha_direct(BurnIn+1:M+BurnIn,:);

    eps_direct = randn(M,H);
    
    y_direct = predict_arch(alpha_direct, y_T, S, H, eps_direct);
%     y_direct = predict_arch_noS(alpha_direct, y_T, S, H, eps_direct);

    PL_direct = sort(fn_PL(y_direct));
    VaR_direct(sim,1) = PL_direct(p_bar*M);
    ES_direct(sim,1) = mean(PL_direct(1:p_bar*M));
    
    fprintf('Direct 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_direct(sim,1), model, algo);
end

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','mit_direct','accept_direct')
end

%% QERMit 1a.: 
algo = 'MitISEM';
kernel_init = @(a) - posterior_arch(a, data, S, true);
kernel = @(a) posterior_arch(a, data, S, true);
[mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);

for sim = 1:N_sim
    fprintf('\nPrelim sim = %i.\n', sim);
    kernel = @(a) posterior_arch(a, data, S, true);

    [alpha1, accept(sim,1)] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept(sim,1), model, algo);
    alpha1 = alpha1(BurnIn+1:M+1000);

    eps1 = randn(M,H);
    y_H = predict_arch(alpha1, y_T, S, H, eps1);
    % get the preliminary VaR estimate as the 100th of the ascendingly sorted percentage loss values
    [PL_H, ind] = sort(fn_PL(y_H));
    VaR_prelim(sim,1) = PL_H(p_bar*M);
    ES_prelim(sim,1) = mean(PL_H(1:p_bar*M));    
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
end

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept')
end

%% Choose the starting point (mu_hl) for the constuction of the approximaton
 % If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
% arch model
        kernel = @(xx) posterior_arch(xx, data, S, true);
        % y_H = y_predict(draw_mm); 
        y_predict = @(draw) predict_arch(draw(:,1), y_T, S, H, draw(:,2:end));  
% cont2.mit.N =10000; % the desired number of high-loss draws         

        [draw_hl, VaR_est, ~, ~] = BigDraw(cont2.mit.N, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat);

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept','draw_hl','VaR_est')
end

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

% cont2.mit.N = 10000;
cont2.mit.Hmax = 10;

mit_hl.mu = mu_hl;
mit_hl.Sigma = Sigma_hl;
mit_hl.df = cont2.mit.dfnc;
mit_hl.p = 1;
% mu_init = mu_hl;
% mit_init = mit_hl;
 
kernel_init = @(a) - posterior_arch_hl(a, data, S, mean(VaR_prelim), true);
kernel = @(a) posterior_arch_hl(a, data, S, mean(VaR_prelim), true);
if (H < 3)
    [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);
else
    [mit2, summary2] = MitISEM_new(mit_hl, kernel, mu_hl, cont2, GamMat);
end
    
if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept','draw_hl','w_hl','lnk_hl','mit2','summary2')
end


%% VaR with standard MitISEM
for sim = 1:N_sim
    resampl_on = false;
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
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);
    PL_opt = fn_PL(y_opt);

    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_IS(sim,1), model, algo);
    fprintf('IS 100*%4.2f%% ES estimate: %6.4f (%s, %s). \n', p_bar, ES_IS(sim,1), model, algo);  
end
kernel = @(xx) posterior_arch_hl(xx, data, S, mean(VaR_prelim), true);


if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','VaR_IS','ES_prelim','ES_IS','mit2','mit1','accept','draw_hl','VaR_est')
end


labels_in = {'prelim','mitisem'};
Boxplot_PMitISEM(VaR_prelim, VaR_IS, ES_prelim, ES_IS, model, algo, H, N_sim, true, labels_in);

