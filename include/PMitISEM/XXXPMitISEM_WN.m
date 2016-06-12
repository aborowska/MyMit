clear all
close all

%% Initialization
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s); 

addpath(genpath('include/'));
plot_on = true;
save_on = true;

x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);

model = 'WN';
algo = 'PMitISEM';

% Artificial, white noise data 
T = 10000;
y = randn(T,1); 
y = y - mean(y);

% sigma is the VARIANCE of the error term, i.e. y_t ~ NID(0, sigma)
sigma_init = 0.9;

% Control parameters for MitISEM (cont) and PMitiISEM (cont2)
MitISEM_Control
cont.mit.dfnc = 5;
cont.mit.N = 10000;
cont.resmpl_on = false;

cont2 = cont;
cont2.mit.iter_max = 5;
cont2.df.range = [1,10];
% hyper parameters for the prior for sigma2(inv. gamma)
a = 1; % if a == 0, then the flat prior is used; if a == 1, then the conjugate prior (inv. gamma)
b = 1; 
% posterior parameters:
a_post = a + T/2;
b_post = b + sum(y)/2;

sim = 1;
N_sim = 20;

% Metropolis-Hastings for the parameters
M = 10000; % number of draws for preliminary and IS computations
BurnIn = 1000;

H = 10; % forecast horizon
p_bar = 0.01;
% d = H+1; % dimension of theta
% partition = [1,3:H+1];
% partition = [1,3,5*(1:49)+2];
% partition = [1:2:251]; 


VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
accept = zeros(N_sim,1);

VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);
time_pmit = zeros(2,1);

% Construct the approximation to the parameter posterior
kernel_init = @(x) - posterior_debug(x, y, a, b, true);
kernel = @(x) posterior_debug(x, y, a, b, true);
tic
[mit1, summary1] = MitISEM_new(kernel_init, kernel, sigma_init, cont, GamMat);
time_pmit(1,1) = toc;



%% Preliminary VaR
tic
for sim = 1:N_sim 
    fprintf('\nVaR prelim iter: %d\n',sim)
    kernel = @(x) posterior_debug(x, y, a, b, true);
    [sigma1, accept(sim,1)] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    fprintf('(%s) MH acceptance rate: %4.2f. \n', model, accept(sim,1));
    sigma1 = sigma1((BurnIn+1):(M+BurnIn));

    %% Future disturbances
    eps_H = randn(M,H); % --> if future disturbances drawn from the target then their weights are 1

    % the returns coresponding to disturbances: y = sqrt(sigma)*eps
    y_H = bsxfun(@times,eps_H,sqrt(sigma1)); 

    % preliminary VaR
    [PL, ind] = sort(fn_PL(y_H));
%     clear y_H
    VaR_prelim(sim,1) = PL(p_bar*M);  
    ES_prelim(sim,1) = mean(PL(1:p_bar*M));    
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
    % theoretical VaR (with the FLAT prior) (or, the likelihood)
%     VaR_true = norminv(p_bar,0,sqrt(H));
%     fprintf('Theoretical 100*%4.2f%% VaR (under flat prior): %6.4f.\n', p_bar,VaR_true);
end
time_pmit(1,1) = time_pmit(1,1) + toc/N_sim;

if plot_on
    Plot_hor_direct(y_H,y(end), VaR_prelim(sim,1),model, save_on);
end

%% Constuct pmit using the high loss draws   
% If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
% arch model
% kernel = @(xx) posterior_arch(xx, data, S, true);
% y_H = y_predict(draw_mm); 
kernel = @(xx) posterior_debug(xx, y, a, b, true);
y_predict = @(draw) bsxfun(@times,draw(:,2:end),sqrt(draw(:,1))); 
% cont2.mit.N = 10000; % the desired number of high-loss draws         
[draw_hl, VaR_est, ~, ~] = BigDraw(cont2.mit.N, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat);

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept','draw_hl','VaR_est')
end

% WEIGHTS to initialise PMitISEM
% future disturbances are generated from the target thus have weights 1
% log kernel evaluation - only for the parameter draws
kernel = @(xx) posterior_debug(xx, y, a, b, true);
lnk_hl = kernel(draw_hl(:,1)); 
% log candidate evaluation
lnd_hl = dmvgt(draw_hl(:,1), mit1, true, GamMat);
% weights
w_hl = lnk_hl - lnd_hl;
w_hl = exp(w_hl - max(w_hl));


%% PMIT ALL
% patial mixture - based on multidimensional stuctures with 4 mit fields
% mit_struct = struct('mu',[],'Sigma',[],'df',[],'p',[]);
% pmit_struct = struct('mu',cell(1,S),'Sigma',cell(1,S),'df',cell(1,S),'p',cell(1,S));

partition = [1,3:H+1];
% partition = [1,4:2:11];
% partition = [1,3,5*(1:49)+2];
% partition = [1:2:251];
S = length(partition);
d = H+1;

fn_const_X = @(a) WN_const_X(a);

kernel = @(x) posterior_debug_hl(x, y, a, b, mean(VaR_prelim), true); 
CV_old = cont.mit.CV_old;
CV_tol = cont.mit.CV_tol;

draw0 = draw_hl;
w0 = w_hl;
lnk0 = lnk_hl; %kernel(draw0);
% clear draw_hl w_hl lnk_hl lnd_hl

cont2.mit.iter_max = 5;
cont.df.range = [1,20];
cont = cont2;
% profile on
% [pmit, CV_mix, CV, iter, pmit_pre, pmit_pre2, pmit_adapt] = PMitISEM(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont2, GamMat);
tic
[pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM_debug(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont, GamMat);
time_pmit(1,1) = time_pmit(1,1) + toc;

% profile off
% profile viewer

 if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','accept','draw_hl','VaR_est','pmit','CV_mix','CV')
 end

% kernel_init = @(x) - posterior_debug(x, y, a, b, true);
% [mu, Sigma] = fn_initopt(kernel_init, sigma_init);
% mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',5);

%% VaR with PMit
tic
for sim = 1:N_sim 
    fprintf('\nVaR IS iter: %d\n',sim)
     
    sigma1 = rmvgt2(M/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
%     sigma1 = rmvgt2(M/2, mit_direct.mu, mit_direct.Sigma, mit_direct.df, mit_direct.p); 
    eps1 = randn(M/2,H); 
    draw1 = [sigma1, eps1];
    
    draw_pmit  = fn_p_rmvgt(M/2, pmit, d, partition, [], fn_const_X);  
    draw_opt = [draw1; draw_pmit];
    clear draw1 %draw_pmit
    
    kernel = @(x) posterior_debug(x, y, a, b, true);
    lnk_opt = kernel(draw_opt(:,1));

    kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
    eps_pdf = sum(kernel(draw_opt(:,2:d)),2);
    lnk_opt = lnk_opt + eps_pdf;

    % optimal weights 
    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,1), mit1, true, GamMat));
%     exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,1), mit_direct, true, GamMat));
    exp_lnd2 = fn_dpmit(draw_opt, pmit, partition, fn_const_X, true, GamMat);

    exp_lnd2 = 0.5*exp(exp_lnd2);
    exp_lnd = exp_lnd1 + exp_lnd2;
    clear exp_lnd1 exp_lnd2
    lnd_opt = log(exp_lnd);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % IS VaR estimation
    y_opt = bsxfun(@times,draw_opt(:,2:d),sqrt(draw_opt(:,1)));  
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);   
  
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_IS(sim,1), model, algo);  
end
time_pmit(2,1) = toc/N_sim;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','summary1','accept',...
        'draw_hl','VaR_est','pmit','CV_mix','CV','VaR_IS','ES_IS','time_pmit')
end

if plot_on
    Boxplot_PMitISEM(VaR_prelim,VaR_IS,ES_prelim,ES_IS,model,algo,H,N_sim,save_on);

    y_pmit = bsxfun(@times,draw_pmit(:,2:d),sqrt(draw_pmit(:,1)));  
    Plot_hor_pmit(y_pmit, y(end), mean(VaR_prelim),model,algo,save_on)
end