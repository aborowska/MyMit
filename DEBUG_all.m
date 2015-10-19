%% Initialisation
% clc
clear all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath('include/');
addpath('.git/');

x_gam = (0:0.00001:50)' + 0.00001;
GamMat = gamma(x_gam);

T = 10000;
y = randn(T,1);


sigma_init = 0.9;
M = 10000;
N_sim = 2;

p_bar = 0.05;
p_bar2 = 0.05;

MitISEM_Control

% hyper params:
a = 1;
b = 1; 
% logkernel
kernel_init = @(x) - posterior_debug(x, y, a, b);
kernel = @(x) posterior_debug(x, y, a, b);
% kernel([-10:1:10]')
% posterior_debug([-1,0,1,1.5,0.2,-1]', [0,1,2,3,4]', 1, 1)
% mu_init = sigma_init;
cont.mit.dfnc = 1;
[mit1_df1, summary1_df1] = MitISEM(kernel_init, kernel, sigma_init, cont, GamMat);

% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s); 
% cont.mit.dfnc = 5;
% [mit1_df5, summary1_df5] = MitISEM(kernel_init, kernel, sigma_init, cont, GamMat);


% draw from posterior
[draw1_df1, lnk1_df1, ~] = fn_rmvgt_robust(M, mit1_df1, kernel);
% [draw1_df5, lnk1_df5, ~] = fn_rmvgt_robust(M, mit1_df5, kernel);

% Posterior, theoretical moments
a_post = a + T/2;
b_post = b + sum(y.^2)/2;
mean_post = b_post/(a_post-1);
var_post = (b_post.^2)/(((a_post-1)^2)*(a_post-2));
std_post = sqrt(var_post);

% MH
[sigma1_df1, accept_df1] = Mit_MH(M+1000, kernel, mit1_df1, GamMat);
fprintf('(df1) MH acceptance rate: %4.2f. \n',accept_df1);
% [sigma1_df5, accept_df5] = Mit_MH(M+1000, kernel, mit1_df5, GamMat);
% fprintf('(df5) MH acceptance rate: %4.2f. \n',accept_df5);

sigma1_df1 = sigma1_df1(1001:M+1000);
% sigma1_df5 = sigma1_df5(1001:M+1000);

% % [h1, ~, bounds1] = autocorr(sigma1_df1,50);
% % % [h5, ~, bounds5] = autocorr(sigma1_df5,50);
% % autocorr(sigma1_df1,50);
% % % autocorr(sigma1_df5,50);
% % L1 = min(find((h1<bounds1(1,1)) & (h1>bounds1(2,1))));
% % % L5 = min(find((h5<bounds5(1,1)) & (h5>bounds5(2,1))));
% % IF1 = 1 + 2*sum(h1(1:L1));
% % % IF5 = 1 + 2*sum(h5(1:L5));


% Future logreturns
eps1_df1 = randn(M,1);
y_T1_df1 = sqrt(sigma1_df1).*eps1_df1; 
% eps1_df5 = randn(M,1);
% y_T1_df5 = sqrt(sigma1_df5).*eps1_df5; 

 
[PL_T1_df1, ind_df1] = sort(fn_PL(y_T1_df1));
VaR_prelim_df1 = PL_T1_df1(p_bar*M);
fprintf('(df1) Preliminary VAR estimate: %6.4f. \n',VaR_prelim_df1);
VaR_const_df1 = PL_T1_df1(p_bar2*M);
fprintf('(df1) VAR estimate to construct mit2: %6.4f. \n',VaR_const_df1);

% [PL_T1_df5, ind_df5]= sort(fn_PL(y_T1_df5));
% VaR_prelim_df5 = PL_T1_df5(p_bar*M);
% fprintf('(df5) Preliminary VAR estimate: %6.4f. \n',VaR_prelim_df5);

draw_df1 = [sigma1_df1, eps1_df1];
% draw_df5 = [sigma1_df5, eps1_df5];

sigma1_hl_df1 = sigma1_df1(ind_df1); sigma1_hl_df1 = sigma1_hl_df1(1:p_bar*M); % mean: 1.0075
% sigma1_hl_df5 = sigma1_df5(ind_df5); sigma1_hl_df5 = sigma1_hl_df5(1:p_bar*M); % mean: 1.0075

eps_hl_df1 = eps1_df1(ind_df1); eps_hl_df1 = eps_hl_df1(1:p_bar*M); % mean: -2.65
% eps_hl_df5 = eps1_df5(ind_df5); eps_hl_df5 = eps_hl_df5(1:p_bar*M); % mean: -2.67

%% High loss df1
kernel_init = @(x) - posterior_debug_hl(x, y, a, b, VaR_const_df1);
kernel = @(x) posterior_debug_hl(x, y, a, b, VaR_const_df1);
mu_hl = [1, -3];

cont2 = cont;
% cont2.mit.Hmax = 2;
cont2.mit.dfnc = 1;
[mit2_df1, summary2] = MitISEM(kernel_init, kernel, mu_hl, cont2, GamMat);
[draw2_df1, lnk2_df1, ~] = fn_rmvgt_robust(M, mit2_df1, kernel);
draw_opt_df1 = [draw_df1; draw2_df1];
lnk_opt_df1 = lnk1_df1 - 0.5*(log(2*pi) + eps1_df1.^2);
lnk_opt_df1 = [lnk_opt_df1; lnk2_df1];

% profile on
exp_lnd1_df1 = 0.5*normpdf(draw_opt_df1(:,2)).*dmvgt(draw_opt_df1(:,1),mit1_df1,false, GamMat);
% profile off
% profile viewer
% 
% profile on
exp_lnd2_df1 = 0.5*dmvgt(draw_opt_df1,mit2_df1,false, GamMat);
% profile off
% profile viewer

exp_lnd_df1 = exp_lnd1_df1 + exp_lnd2_df1;
lnd_opt_df1 = log(exp_lnd_df1);
w_opt_df1 =  fn_ISwgts(lnk_opt_df1, lnd_opt_df1, false);
y_opt_df1 = sqrt(draw_opt_df1(:,1)).*draw_opt_df1(:,2);
dens = struct('y',y_opt_df1,'w',w_opt_df1,'p_bar',p_bar);
IS_estim_df1 = fn_PL(dens, 1);
VaR_IS_df1 = IS_estim_df1(1,1);
ES_IS_df1 = IS_estim_df1(1,2);

PL_opt_df1 = fn_PL(y_opt_df1);
[PL_opt_h1_df1, ind] = sort(PL_opt_df1);
hold on
plot(PL_opt_h1_df1)
pos =  max(find(PL_opt_h1_df1<=VaR_IS_df1));
scatter(pos, VaR_IS_df1,'MarkerFaceColor','red')
hold off


%% High loss df5
kernel_init = @(x) - posterior_debug_hl(x, y, a, b, VaR_prelim_df5);
kernel = @(x) posterior_debug_hl(x, y, a, b, VaR_prelim_df5);
mu_hl = [1, -3];

cont2 = cont;
% cont2.mit.Hmax = 2;
cont2.mit.dfnc = 5;
[mit2_df5, summary2_df5] = MitISEM(kernel_init, kernel, mu_hl, cont2, GamMat);
[draw2_df5, lnk2_df5, ~] = fn_rmvgt_robust(M, mit2_df5, kernel);
draw_opt_df5 = [draw_df5; draw2_df5];
lnk_opt_df5 = lnk1_df5 - 0.5*(log(2*pi) + eps1_df5.^2);
lnk_opt_df5 = [lnk_opt_df5; lnk2_df5];

% profile on
exp_lnd1_df5 = 0.5*normpdf(draw_opt_df5(:,2)).*dmvgt(draw_opt_df5(:,1),mit1_df5,false, GamMat);
% profile off
% profile viewer

% profile on
exp_lnd2_df5 = 0.5*dmvgt(draw_opt_df5,mit2_df5,false, GamMat);
% profile off
% profile viewer

exp_lnd_df5 = exp_lnd1_df5 + exp_lnd2_df5;
lnd_opt_df5 = log(exp_lnd_df5);
w_opt_df5 =  fn_ISwgts(lnk_opt_df5, lnd_opt_df5, false);
y_opt_df5 = sqrt(draw_opt_df5(:,1)).*draw_opt_df5(:,2);
dens = struct('y',y_opt_df5,'w',w_opt_df5,'p_bar',p_bar);
IS_estim_df5 = fn_PL(dens, 1);
VaR_IS_df5 = IS_estim_df5(1,1);
ES_IS_df5 = IS_estim_df5(1,2);

PL_opt_df5= fn_PL(y_opt_df5);
[PL_opt_h1_df5, ind] = sort(PL_opt_df5);
hold on
plot(PL_opt_h1)
pos =  max(find(PL_opt_h1<=VaR_IS));
scatter(pos, VaR_IS,'MarkerFaceColor','red')
hold off


%%  
if plot_on
    figure(3)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0.85:0.01:1.15;
    yy = -5:0.01:5;
    Mit2 = MitISEM_plot(mit2_df5, 3, xx, yy, GamMat);
    title('(White noise logreturns) Approximation to the high loss density $$q_{2,Mit}(\sigma^2,\varepsilon_{T+1})$$.')
    xlabel('$$\sigma^2$$')
    ylabel('$$\varepsilon_{T+1}$$')
%    plotTickLatex2D;
    set(gca,'TickLabelInterpreter','latex')
    campos([-3,-30,150]);
    if print_on
        name = 'figures/q2_mit_mitisem.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
end