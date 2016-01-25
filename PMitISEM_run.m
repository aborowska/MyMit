clear all
close all

addpath(genpath('include/'));
% d = 4;
% 
% mit_adapt.mu = [1 2 3 4];
% mit_adapt.Sigma = reshape(eye(4), 1, d^2);
% mit_adapt.df = 5;
% mit_adapt.p = 1;
% 
% 
% partition = [1,4];
% 
% S = length(partition);
% 
% 
%    
% fn_const_X = @(a) WN_const_X(a)
% 
% [pmit, summary] = PMitISEM(kernel_init, kernel, fn_const_X, mu_init, partition, cont, GamMat)

%% WN

model = 'WN';
algo = 'MitISEM';
p_bar = 0.01;


x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);
%% Artificial, white noise data 
T = 10000;
y = randn(T,1); 
y = y - mean(y);

%% Parameter draws 
% sigma is the VARIANCE of the error term, i.e. y_t ~ NID(0, sigma)
% (should be sigma2)
sigma_init = 0.9;

% Control parameters for  MitISEM 
MitISEM_Control
cont.mit.dfnc = 5;
cont.mit.N = 10000;
cont.resmpl_on = false;
% hyper parameters for the prior for sigma2(inv. gamma)
a = 1; % if a == 0, then the flat prior is used; if a == 1, then the conjugate prior (inv. gamma)
b = 1; 

kernel_init = @(x) - posterior_debug(x, y, a, b, true);
kernel = @(x) posterior_debug(x, y, a, b, true);
[mit1, summary1] = MitISEM_new(kernel_init, kernel, sigma_init, cont, GamMat);

% Metropolis-Hastings for the parameters
M = 1000000;
[sigma1, accept ] = Mit_MH(M+1000, kernel, mit1, GamMat);
fprintf('(%s) MH acceptance rate: %4.2f. \n', model, accept);
sigma1 = sigma1(1001:M+1000);
% lnk = kernel(sigma1);


%% Future disturbances
H = 10; % forecast horizon
% y_H = randn(M,H); % --> if future disturbances drawn from the prior then
% their weights are 1
  
% to have some stating mixture of t's --> construct a mit for errors and
% draw from that

% mit approximation to the normal distribution
mit_init.mu = 0;
mit_init.Sigma = 1;
mit_init.df = 1;
mit_init.p = 1;
kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
mu_init = 1;
[mit2, summary1] = MitISEM_new(mit_init, kernel, mu_init, cont, GamMat);

% draw from mit2 H times
[eps_H, ~ ] = fn_rmvgt_robust(M*H, mit2, kernel, false);
eps_H = reshape(eps_H,M,H);
% lnk = reshape(lnk,M,H);
y_H = bsxfun(@times,eps_H,sqrt(sigma1)); % the coresponding returns 

[PL, ind] = sort(fn_PL(y_H));
VaR_prelim = PL(p_bar*M);  
ES_prelim = mean(PL(1:p_bar*M));    
fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim, model, algo);
clear PL y_H 


draw_hl = [sigma1, eps_H];
clear eps_H
draw_hl = draw_hl(ind,:);
draw_hl = draw_hl(PL<=VaR_prelim,:);  

% log kernel evaluation
kernel = @(x) posterior_debug_hl(x, y, a, b, Inf, true); % to skip the condition that the losses are below VaR_prelim
lnk = kernel(draw_hl); 
ind = find(lnk~=-Inf);
lnk = lnk(ind,:);     
draw_hl = draw_hl(ind,:);

% log candidate evaluation
lnd = dmvgt(draw_hl(:,1), mit1, true, GamMat);
for h = 1:H
    lnd = lnd + dmvgt(draw_hl(:,h+1), mit2, true, GamMat);
end
w_hl =  fn_ISwgts(lnk, lnd, false);


%% Use draw_hl and w_hl to start the Partial MitISEM
partition = 1:H+1;
d = 11;
S = length(partition);
 
% 
% [mu_hl, Sigma_hl] = fn_muSigma(draw_hl(:,1), w_hl);
% mit_hl.mu = mu_hl;
% mit_hl.Sigma = Sigma_hl;
% mit_hl.df = 1;
% mit_hl.p = 1;
% 
% pmit(S) = mit_hl; % preallocate
% 
% for s = 1:S
%     [mu_hl, Sigma_hl] = fn_muSigma(draw_hl(:,s), w_hl);
%     pmit(s).mu = mu_hl;
%     pmit(s).Sigma = Sigma_hl;
%     pmit(s).df = 1;
%     pmit(s).p = 1;
% end



% patial mixture
pmit(1) = mit1;
for s = 2:S
	pmit(s) = mit2;
end

kernel = @(x) posterior_debug_hl(x, y, a, b, Inf, true); % to skip the condition that the losses are below VaR_prelim

lnk = kernel(draw_hl);
% lnd = 0;
% for s = 1:S
%     [s1, s2] = fn_partition_ends(partition, d, s);
%     lnd = lnd + dmvgt(draw_hl(:,s1:s2), pmit(s), true, GamMat);
% end
% w =  fn_ISwgts(lnk, lnd, false);
pmit(1) = mit1;
fn_const_X = @(a) WN_const_X(a);
theta = draw_hl;
norm = cont.mit.norm;

[pmit, summary] = PMitISEM(kernel_init, kernel, fn_const_X, mu_init, partition, cont, GamMat);

