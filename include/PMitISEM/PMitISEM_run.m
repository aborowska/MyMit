clear all
close all

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s); 

addpath(genpath('include/'));
plot_on = false;

%% Constnts 

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

wn_plot0

% Metropolis-Hastings for the parameters
M = 100000;
[sigma1, accept ] = Mit_MH(M+1000, kernel, mit1, GamMat);
fprintf('(%s) MH acceptance rate: %4.2f. \n', model, accept);
sigma1 = sigma1(1001:M+1000);

%% Future disturbances
H = 10; % forecast horizon
eps_H = randn(M,H); % --> if future disturbances drawn from the target then their weights are 1
  
%% to have some stating mixture of t's --> construct a mit for errors and
% % draw from that
% % mit approximation to the normal distribution
% mit_init.mu = 0;
% mit_init.Sigma = 1;
% mit_init.df = 1;
% mit_init.p = 1;
% kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
% mu_init = 1;
% [mit2, summary1] = MitISEM_new(mit_init, kernel, mu_init, cont, GamMat);
% 
% % draw from mit2 H times
% [eps_H, ~ ] = fn_rmvgt_robust(M*H, mit2, kernel, false);
% eps_H = reshape(eps_H,M,H);
% % lnk = reshape(lnk,M,H);
%%
% the returns coresponding to disturbances
y_H = bsxfun(@times,eps_H,sqrt(sigma1)); 
% preliminary VaR
[PL, ind] = sort(fn_PL(y_H));
VaR_prelim = PL(p_bar*M);  
ES_prelim = mean(PL(1:p_bar*M));    
fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim, model, algo);
% theoretical VaR (with the FLAT prior) (or, the likelihood)
fprintf('Theoretical 100*%4.2f%% VaR (under flat prior): %6.4f.\n', p_bar, norminv(p_bar,0,sqrt(H)));
 
if plot_on
    hold on
    plot(PL)
    plot(sort(sqrt(10)*randn(M,1)),'k')
    hold off
end
% clear PL y_H 

% High loss draws = the target of the truncated H-days-ahead return
% distibution
draw = [sigma1, eps_H];
clear eps_H
draw_hl = draw(ind,:);
draw_hl = draw_hl(PL<=VaR_prelim,:);  

% Is there any difference in the emipircal distribution of high loss
% distubances in the subsequitive horizon periods?
if plot_on
    figure(1)
    for ii = 2:H+1
        subplot(2,5,ii-1) 
    %     hist(draw_hl(:,ii),20)
        [f,xi] = ksdensity(draw_hl(:,ii));
        hold on 
        histnorm(draw_hl(:,ii),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
        Me = mean(draw_hl(:,ii));
        plot([Me Me], [0 min(2*max(f),1)],'k','LineWidth',2);
        hold off
        str = sprintf('h = %i, mean = %6.4f',ii-1,Me);
        title(str)
    end
end

% WEIGHTS to initialise PMitISEM
% future disturbances are generated from the target thus have weights 1
% log kernel evaluation
kernel = @(x) posterior_debug(x, y, a, b, true);
lnk = kernel(draw(:,1));
lnk_hl = kernel(draw_hl(:,1)); 

% log candidate evaluation
lnd = dmvgt(draw(:,1), mit1, true, GamMat);
lnd_hl = dmvgt(draw_hl(:,1), mit1, true, GamMat);

% log weights
w = lnk - lnd;
w = exp(w - max(w));

w_hl = lnk_hl - lnd_hl;
w_hl = exp(w_hl - max(w_hl));



%% 
% log kernel evaluation
% kernel = @(x) posterior_debug_hl(x, y, a, b, Inf, true); % to skip the condition that the losses are below VaR_prelim
% lnk = kernel(draw_hl); 
% ind = find(lnk~=-Inf);
% lnk = lnk(ind,:);     
% draw_hl = draw_hl(ind,:);
% 
% % log candidate evaluation
% lnd = dmvgt(draw_hl(:,1), mit1, true, GamMat);
% for h = 1:H
%     lnd = lnd + dmvgt(draw_hl(:,h+1), mit2, true, GamMat);
% end
% w_hl =  fn_ISwgts(lnk, lnd, false);

%% Use draw_hl and w_hl to start the Partial MitISEM
partition = 1:H+1;
d = 11;
S = length(partition);
 
%% [mu_hl, Sigma_hl] = fn_muSigma(draw_hl(:,1), w_hl);
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

%% patial mixture
pmit(1) = mit1;

for s = 2:S    
    [mu_adapt, Sigma_adapt] = fn_muSigma(draw_hl(:,s), w_hl);
    mit_adapt.mu = mu_adapt;
    mit_adapt.Sigma = Sigma_adapt;
    mit_adapt.df = cont.mit.dfnc;
    mit_adapt.p = 1;

    pmit(s) = fn_optimt(draw_hl(:,s), mit_adapt, w_hl, cont, GamMat);
end

% draw from pmit
draw_pmit = zeros(M,S);

kernel = @(x) posterior_debug(x, y, a, b, true);
[draw_pmit(:,1), lnk] = fn_rmvgt_robust(M, pmit(1), kernel, false);
lnd = dmvgt(draw_pmit(:,1), pmit(1), true, GamMat);
kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
for ss = 2:S
    [draw_pmit(:,ss), lnk_s] = fn_rmvgt_robust(M, pmit(ss), kernel, false);
    lnk = lnk + lnk_s;
    lnd = lnd + dmvgt(draw_pmit(:,ss), pmit(ss), true, GamMat);
end   
w_pmit = lnk - lnd;

y_pmit = bsxfun(@times,draw_pmit(:,2:S),sqrt(draw_pmit(:,1))); 
% y_pmit = sort(fn_PL(return_pmit));
% the returns coresponding to disturbances
y_H = bsxfun(@times,eps_H,sqrt(sigma1)); 
% preliminary VaR
[PL, ind] = sort(fn_PL(y_H));

y_opt = [y_H; y_pmit];
% w_opt = [w; w_pmit];
dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
IS_estim = fn_PL(dens, 1);
   
if plot_on
    figure(2)
    hold on
    plot(PL)
    plot(sort(return_pmit),'r')
    plot(ones(M,1)*VaR_prelim,'m')
    hold off
    legend('Direct returns','MitISEM returns','VaR prelim')
end
% IS VaR estimation











% APPLY ISEM TO MARGINAL
        pmit(1) = fn_optimt(draw_hl(:,1), mit1, w_hl, cont, GamMat);
        lnd = 
        w_current = lnk - 
        
        ind_w = fn_select(w_hl, cont.mit.ISpc);
        theta_nc = draw_hl(ind_w,:);
        w_nc = w_hl(ind_w);

        % NEW COMPONENT
        % compute new component's mode and scale from IS weights
        mit_nc.p = cont.mit.pnc;
        mit_nc.df = cont.mit.dfnc;
        [mit_nc.mu, mit_nc.Sigma] = fn_muSigma(theta_nc(:,1), w_nc);

        % COBINE OLD AND NC
        % combine the old mixture mit_new and the new component mit_nc      
        mit_new = fn_updateMit(pmit(1), mit_nc); 

        % UPDATE COMBINED
        % update mode, scale and df  of all mixture components
        mit_new = fn_optimt(draw_hl(:,1), mit_new, w_hl, cont, GamMat);
   








fn_const_X = @(a) WN_const_X(a);
theta = draw_hl;
norm = cont.mit.norm;

[pmit, summary] = PMitISEM(kernel_init, kernel, fn_const_X, mu_init, partition, cont, GamMat);

