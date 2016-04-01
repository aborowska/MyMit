%%
clear all
addpath('include/');

data = csvread('GSPC_ret.csv');
data = 100*data;
% QERMit ARCH 
I_arch = find(data<=-5.5);
data = data(1:576,1);
S = var(data);

theta = 0.03;
M = 10000;


% Control
cont.mit.N = 1000; %1e5;
cont.mit.Hmax = 10;

cont.mit.CV_tol = 0.1; %0.1
cont.mit.AR_tol = 0.1; %0.1
cont.mit.CV_old = 100;
cont.mit.AR_old = 100;
cont.mit.meth = 'CV';

cont.mit.ISpc = 0.1;
cont.mit.pnc = 0.1; % probability of a new component
cont.mit.tol_pr = 0;

cont.EM.maxit = 1000;
cont.EM.tol = 0.001;

cont.df.maxit = 1000;
cont.df.opt = true;
cont.df.range = [0.01,40];
cont.df.tol = eps^(0.25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
L = true;
kernel_init = @(a) - posterior_arch(a, data, S, L);
kernel = @(a) posterior_arch(a, data, S, L);
mu0 = theta;

[mit, summary] = MitISEM(kernel_init, kernel, mu0, cont);
mit1 = mit;
% for QERMit 2.
[draw1, lnk1, ind_red1, time1] = fn_rmvgt_robust(M/2, mit1, kernel);
draw1 = [draw1,randn(M/2,1)];

% Fig. 5
figure(1)
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');
Mit1 = MitISEM_plot(mit1, 2);
title('Approximation to the posterior density')
plotTickLatex2D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% QERMit 1b.: generate set opf draws of alpha using independence MH with
% candiate from MitISEM; then simulate returns from normal with variance
% based on the draw of alpha 

tic
[alpha, accept] = MitISEM_MH(M, kernel, mit);
time_MH = toc;
y_T = data(length(data));
stdev =  sqrt(S+(y_T^2-S)*alpha);
y_T1 = stdev.*randn(M,1);

% get the preliminary VaR estimate as the 100th of the ascendingly sorted percentage loss values
y_T1 = sort(y_T1);
VaR_prelim = y_T1(100);
%-5.3749

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% QERMit 1c.: get mit approximation of the conditional joint density of
% parameters and future returns given the returns are below VaR_prelim
% approximation of the joiunt high loss distribution
% here: not future returns but furure disturbances  (varepsilons)
c = 100*log(1 + VaR_prelim/100);
eps_bord = c./stdev;
high_loss_subspace = [alpha, eps_bord];
% high loss error - from truncated normal
% eps_hl = norminv(rand(M,1).*normcdf(eps_bord));
% d = posterior_arch_hl(theta, data, S, eps_bord, L)
 

% Fig. 6.1
figure(2)
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');
hold on
scatter(draw1(:,1),draw1(:,2),10,'k.')
xlabel('$$\alpha_1$$')
ylabel('$$\varepsilon_{T+1}$$')
scatter(high_loss_subspace(:,1),high_loss_subspace(:,2),10,'r.')
hold off
title('Joint density $$p(\alpha_{1},\varepsilon_{T+1}|y)$$ and hight loss subspace')
plotTickLatex2D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
L = true;
kernel_init = @(a) - posterior_arch_hl(a, data, S, VaR_prelim, L);
kernel = @(a) posterior_arch_hl(a, data, S, VaR_prelim, L);
% theta = [alpha, eps]
theta = [0.15, -3];
mu0 = theta;
[mit, summary2] = MitISEM(kernel_init, kernel, mu0, cont);

mit2 = mit;
% Mit = MitISEM_plot(mit2,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% QERMit 2: use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
% to estiamte VaR and ES for theta and y (or alpha in eps)
[draw2, lnk2, ind_red2, time2] = fn_rmvgt_robust(M/2, mit2, kernel);

draw_opt = [draw1; draw2];
