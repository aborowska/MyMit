% https://e76d6ebf22ef8d7e079810f3d1f82ba1e5f145d5.googledrive.com/host/0B2GQktu-wcTiWDB2R2t2a2tMUG8/
% https://www.chrisstucchio.com/blog/2013/bayesian_bandit.html
% http://tdunning.blogspot.nl/2012/02/bayesian-bandits.html
% https://dataorigami.net/blogs/napkin-folding/79031811-multi-armed-bandits

clear all
close all
cd ../
addpath(genpath('include/'));
cd ./Bandit

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
 
x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

% true theta
theta_true = [0.5, 0.49, 0.4];%, 0.3];
K = size(theta_true,2);
G = 10000;
BurnIn = 1000;

% simulated data
% [S, F, theta] = bernoulli_sampling(G, K, theta_true); %Thomposon sampling
% [~, arm_selected] = max(theta,[],2);

N = [100, 120];%, 130]; % no. of trials
Y = [sum(rand(N(1,1),1)<theta_true(1,1)),...
     sum(rand(N(1,2),1)<theta_true(1,2))];%,...
    sum(rand(N(1,3),1)<theta_true(1,3))]; % no. of successes

cont = MitISEM_Control;
cont.mit.dfnc = 5;
cont.mit.N = 10000;
cont.resmpl_on = false;

mu_init = [0.5, 0.5];%, 0.5];

kernel_init = @(xx) - posterior_binomial(N, Y, xx);
kernel = @(xx) posterior_binomial(N, Y, xx);
kernel(mu_init);
[mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
  
xx = 0:0.01:1;
Mit1 = MitISEM_plot(mit1, 3, xx, xx, GamMat);
     

% MH
kernel = @(xx) posterior_binomial(N, Y, xx);
[theta1, accept] = Mit_MH(G+BurnIn, kernel, mit1, GamMat);
fprintf('MH acceptance rate: %4.2f. \n', accept);
theta1 = theta1(BurnIn+1:G+BurnIn,:);

regret1 = regret_binomial(theta1);
p_bar = 0.95;
[regret1_sort, ind] = sort(regret1);
PVR_prelim = regret1_sort(round(p_bar*G));
PVR_cum_prelim = mean(regret1_sort(round(p_bar*G:G)));   


theta1_sort = theta1(ind,:);
theta_hl = theta1_sort(round(p_bar*G)+1:G,:);
for ii = 1:K
    figure(ii)
    subplot(1,2,1)
    hist(theta1(:,ii))
   
    subplot(1,2,2)
    hist(theta_hl(:,ii))
end

kernel = @(xx) posterior_binomial(N, Y, xx);
% log kernel
lnk_hl = kernel(theta_hl); 
% log candidate evaluation
lnd_hl = dmvgt(theta_hl, mit1, true, GamMat);
% importance weights
w_hl = lnk_hl - lnd_hl;
w_hl = exp(w_hl - max(w_hl));

[mu_hl, Sigma_hl] = fn_muSigma(theta_hl, w_hl);
cont2 = cont;
mit_hl.mu = mu_hl;
mit_hl.Sigma = Sigma_hl;
mit_hl.df = cont2.mit.dfnc;
mit_hl.p = 1;
% mu_init_hl = mean(theta_hl,1); 

% kernel_init = @(xx) - posterior_binomial_hl(N, Y, xx, PVR_prelim);
kernel = @(xx) posterior_binomial_hl(N, Y, xx, PVR_prelim);
cont2.mit.Hmax = 2;
% [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_init_hl, cont2, GamMat);
[mit2, summary2] = MitISEM_new(mit_hl, kernel, mu_hl, cont2, GamMat);

Mit2 = MitISEM_plot(mit2, 3, xx, xx, GamMat);



% IS

draw1 = rmvgt2(G/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
draw2 = rmvgt2(G/2, mit2.mu, mit2.Sigma, mit2.df, mit2.p); 
draw_opt = [draw1; draw2];

% IS weights
kernel = @(xx) posterior_binomial(N, Y, xx);
lnk = kernel(draw_opt);

exp_lnd1 = 0.5*dmvgt(draw_opt, mit1, false, GamMat);
exp_lnd2 = 0.5*dmvgt(draw_opt, mit2, false, GamMat);
exp_lnd = exp_lnd1 + exp_lnd2;
lnd = log(exp_lnd);

w_opt = fn_ISwgts(lnk, lnd, false);

regret_opt = regret_binomial(draw_opt);
regret_opt_sort = sort(regret_opt);


dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    PVR_IS = IS_estim(1,1);
    PVR_cum_IS = IS_estim(1,2);
