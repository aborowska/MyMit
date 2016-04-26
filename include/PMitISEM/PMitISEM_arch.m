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
% Metropolis-Hastings for the preliminary
M = 10000;
BurnIn = 1000;

N_sim = 20;
p_bar = 0.01;
H = 10; % forecast horizon
% d = H+1; % dimension of theta
% partition = [1,3:H+1];

plot_on = false;
print_on  = false;
plot_on2 = false;
save_on = false;

MitISEM_Control
cont.mit.dfnc = 5;
cont.mit.N = 10000;
cont.resmpl_on = false;

cont2 = cont;
% cont2.mit.Hmax = 10;
% cont2.mit.dfnc = 5;
% cont2.mit.CV_tol = 0.1;
cont2.df.range = [1, 40];

VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
accept = zeros(N_sim,1);

VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);

%% QERMit 1a.: 
kernel_init = @(a) - posterior_arch(a, data, S, true);
kernel = @(a) posterior_arch(a, data, S, true);
[mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);


for sim = 1:N_sim
    fprintf('\nPrelim sim = %i.\n', sim);
    kernel = @(a) posterior_arch(a, data, S, true);

    [alpha1, accept(sim,1)] = Mit_MH(M+1000, kernel, mit1, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept(sim,1), model, algo);
    alpha1 = alpha1(1001:M+1000);

    eps1 = randn(M,H);

    y_H = predict_arch(alpha1, y_T, S, H, eps1);

    % get the preliminary VaR estimate as the 100th of the ascendingly sorted percentage loss values
    [PL_H, ind] = sort(fn_PL(y_H));
    VaR_prelim(sim,1) = PL_H(p_bar*M); % VaR_prelim = 0; VaR_prelim = Inf;
    ES_prelim(sim,1) = mean(PL_H(1:p_bar*M));    
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
end

if save_on
    name = ['results/PMitISEM/',model,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept')
end


if plot_on
    Plot_hor_direct(y_H,y_T,VaR_prelim(sim,1),model);
end

% Choose the starting point (mu_hl) for the constuction of the approximaton
%         alpha1_hl = alpha1(ind); % ind sorts the draws ascendingly in the corresponding profit/losses
%         eps_hl = eps1(ind,:); 
%         draw_hl = [alpha1_hl, eps_hl]; 
%         draw_hl =  draw_hl((find(PL_H < mean(VaR_prelim))),:);  

% If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
% arch model
        kernel = @(xx) posterior_arch(xx, data, S, true);
        % y_H = y_predict(draw_mm); 
        y_predict = @(draw) predict_arch(draw(:,1), y_T, S, H, draw(:,2:end));  

        [draw_hl, VaR_est, ~, ~] = BigDraw(M, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat);



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
kernel = @(xx) posterior_arch_hl(xx, data, S, mean(VaR_prelim), true);
CV_old = cont.mit.CV_old;
CV_tol = cont.mit.CV_tol;

draw0 = draw_hl;
w0 = w_hl;
lnk0 = kernel(draw0);
clear draw_hl w_hl lnk_hl lnd_hl
[pmit, CV_mix, CV, iter] = PMitISEM(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont, GamMat);

if save_on
    name = ['results/PMitISEM/',model,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','pmit')
end


%% VaR with PMit
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
    [s1, s2] = fn_partition_ends(partition, d, 1);
    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,s1:s2), mit1, true, GamMat));
    exp_lnd2 = fn_dpmit(draw_opt, pmit, partition, fn_const_X, true, GamMat);

    exp_lnd2 = 0.5*exp(exp_lnd2);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % IS VaR estimation
    y_opt = predict_arch(draw_opt(:,1), y_T, S, H, draw_opt(:,2:H+1));  
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);   
  
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_IS(sim,1), model, algo);  
end


if save_on
    name = ['results/PMitISEM/',model,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','VaR_IS','ES_prelim','ES_IS','pmit')
end


if plot_on
    Boxplot_PMitISEM(VaR_prelim,VaR_IS,ES_prelim,ES_IS,model,H,N_sim);
    
    y_pmit = predict_arch(draw_pmit(:,1), y_T, S, H, draw_pmit(:,2:H+1));  
    Plot_hor_pmit(y_pmit, y_T, mean(VaR_prelim),model)
end
 