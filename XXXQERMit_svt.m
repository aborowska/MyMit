 %% Initialisation
 
clear all
addpath(genpath('include/'));
addpath('results/');

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);
prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -2.5*log(0.025), -log(gamma(2.5))];

logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;
    
model = 'svt';
y = csvread('IBM_ret.csv');
y = 100*y;

% fprintf('IBM logreturns kurtosis: %6.4f, and skewness: %6.4f.\n', kurtosis(data), skewness(data));

% http://faculty.chicagobooth.edu/nicholas.polson/research/papers/jpr2.pdf
T = size(y,1);

plot_on = true;
print_on = true;

if plot_on
	figure(1)
	set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
    xx = linspace(2007,2012,T);
    plot(xx, y)
    set(gca,'XTickLabel',{'2007',' ','2008',' ','2009',' ','2010',' ','2011',' ','2012'})
	title('IBM log-returns $$y$$')
%    plotTickLatex2D;
	set(gca,'TickLabelInterpreter','latex')
	if print_on
		name = 'figures/sv_data.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end

EMitISEM_Control
cont.mit.dfnc = 5;
N = cont.mit.N;
 
par_NAIS_init.b = zeros(T,1);
par_NAIS_init.C = ones(T,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1a.:
% theta = [c, phi, sigma2, nu]
% mu_init = [0.5, 0.98, 0.15^2, 7];
mu_init = [0.5, 0.98, 0.15^2, 7];

d = size(mu_init,2);
% SML ... <<<<<<<<<<<<<<<<<<<<<<
load('SMLt_ibm.mat', 'par_SV_opt', 'hess_SV_opt') 

Sigma = inv(hess_SV_opt);
Sigma(1,1) = 1;
% Sigma(1:3,1:3) = Sigma(1:3,1:3)/10;
Sigma(4,:) = Sigma(4,:)/20;
Sigma(:,4) = Sigma(:,4)/20;
Sigma = reshape(Sigma,1,d^2);


mit_init.mu = par_SV_opt;
mit_init.Sigma = Sigma;
mit_init.p = 1;
mit_init.df = 5;

mit = mit_init;
kernel = @(a) prior_svt(a,true,prior_const);
% load('sv_mitisem_dom_all.mat', 'mit_adapt', 'mit_new');
% [theta_adapt, prior_adapt] = fn_rmvgt_robust(N, mit_adapt, kernel);
% [theta_new, prior_new] = fn_rmvgt_robust(N, mit_new, kernel);
% [mit1, theta1, x1, w1, par_NAIS1, summary1] = EMitISEM(mit_init, kernel,  cont, GamMat, prior_const);
% clear GamMat
% save('results/sv_mitisem.mat');
% x_gam = (0:0.00001:100)'+0.00001;
% GamMat = gamma(x_gam);
% 
% load('sv_mitisem_dom_all.mat', 'mit_new', 'theta', 'x', 'w', 'b', 'C'); 
% par_NAIS.b = b';
% par_NAIS.C = C';
% clear b C
[theta1, prior1] = fn_rmvgt_robust(N, mit, kernel);
[x_1, lng_y1, lnw_x1, par_NAIS_1] = fn_nais_wgts(y, theta1, par_NAIS_init, cont.nais);
lnk1 = lng_y1 + lnw_x1 + prior1;
lnd1 = dmvgt(theta1, mit, true, GamMat);
lnw1 = lnk1 - lnd1;

ind_nan = find(isnan(lng_y1)~=1);
N = length(ind_nan);
x_1 = x_1(ind_nan,:);
theta1 = theta1(ind_nan,:);
prior1 = prior1(ind_nan,:);

lng_y1 = lng_y1(ind_nan,1);
lnw_x1 = lnw_x1(ind_nan,1);
lnk1 = lnk1(ind_nan,1);
lnd1 = lnd1(ind_nan,1);
lnw1 = lnw1(ind_nan,1);
par_NAIS_1.b = par_NAIS_1.b(:,ind_nan);
par_NAIS_1.C = par_NAIS_1.C(:,ind_nan);

w_nais = fn_ISwgts(prior1, lnd1, true);

b_pmean_1 = sum(repmat(prior1',T,1).*par_NAIS_1.b,2)/sum(prior1);
C_pmean_1 = sum(repmat(prior1',T,1).*par_NAIS_1.C,2)/sum(prior1);
b_mean_1 = mean(par_NAIS_1.b,2);
C_mean_1 = mean(par_NAIS_1.C,2);
b_wmean_1 = sum(repmat(w_nais',T,1).*par_NAIS_1.b,2);
C_wmean_1 = sum(repmat(w_nais',T,1).*par_NAIS_1.C,2);


if plot_on
    figure(2)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
	set(gcf,'defaulttextinterpreter','latex');
    
    subplot(2,2,1)  
    [f,xi] = ksdensity(theta1(:,1));
    hold on
    histnorm(theta1(:,1),20)
    h = findobj(gca, 'Type','patch');
    set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
    plot(xi,f,'r');
    hold off
    title('$$c$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(2,2,2)
    [f,xi] = ksdensity(theta1(:,2));
    hold on
    histnorm(theta1(:,2),20)
    h = findobj(gca, 'Type','patch');
    set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
    plot(xi,f,'r');
    hold off
    title('$$\phi$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(2,2,3)
    [f,xi] = ksdensity(theta1(:,3));
    hold on
    histnorm(theta1(:,3),20)
    h = findobj(gca, 'Type','patch');
    set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
    plot(xi,f,'r');
    hold off
    title('$$\sigma^{2}_{\eta}$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(2,2,4)
    [f,xi] = ksdensity(theta1(:,4));
    hold on
    histnorm(theta1(:,4),20)
    h = findobj(gca, 'Type','patch');
    set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
    plot(xi,f,'r');
    hold off
    title('$$\nu$$')
	set(gca,'TickLabelInterpreter','latex')
    
%     suptitle('Empirical approximate posterior parameter density $$q_{1,\zeta}(\theta|y)y$$.' )
    if print_on
		name = 'figures/svt_draws_param_init.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(3)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
	set(gcf,'defaulttextinterpreter','latex');
    
    subplot(2,1,1)
    hold on
    plot(b_wmean_1)    
%     plot(b_mean_1) 
    hold off 
    title('$$b$$')
     set(gca,'TickLabelInterpreter','latex')

    subplot(2,1,2)
    hold on
    plot(C_wmean_1)    
%     plot(C_mean_1) 
    hold off 
    title('$$C$$')
	set(gca,'TickLabelInterpreter','latex')
    
%     suptitle('Average NAIS parameters given the draws from $$q_{1,\zeta}(\theta|y)y$$')

	if print_on
		name = 'figures/svt_NAIS_param_init.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end



%% QERMit 1b.:
% generate set opf draws of theta using independence MH with
% candiate from MitISEM; then simulate returns based on the draw of theta 
            % [theta, x, accept] = EMit_MH(M+1000, y, kernel, mit, par_NAIS, cont, GamMat);
            % fprintf('(MitISEM) MH acceptance rate: %6.4f. \n',accept);
            % save('results/sv_mitisem_theta.mat','theta', 'accept');
            % theta = theta(1001:M+1000,:);


%% High loss, 1 days horizon
% approximate the high loss distribution of (theta,eps*,eta*) where 
% eps*={eps_T+1,...,eps_T+hp}
% eta*={eta_T+1,...,eta_T+hp}
p_bar = 0.01;
hp = 1;

c1 = theta1(:,1);
phi1 = theta1(:,2);
sigma21 = theta1(:,3);
nu1 = theta1(:,4);
rho1 = (nu1-2)./nu1;
eta_h1 = randn(N,1);
eps_h1 = trnd(repmat(nu1,1,hp));

x_h1 = c1 + phi1.*(x_1(:,end) - c1) + sqrt(sigma21).*eta_h1;
y_h1 = sqrt(rho1).*exp(0.5*x_h1).*eps_h1;

[PL_h1, ind] = sort(fn_PL(y_h1));
eps_hl_init = eps_h1(ind,:); 
eps_hl_init = eps_hl_init(1:floor(p_bar*N),:);
eta_hl_init = eta_h1(ind,:); 
eta_hl_init = eta_hl_init(1:floor(p_bar*N),:);
theta_hl_init = theta1(ind,:);
theta_hl_init = theta_hl_init(1:floor(p_bar*N),:);
VaR_prelim = PL_h1(floor(p_bar*N));
fprintf('hp = %i, VaR_prelim = %4.5f. \n', hp, VaR_prelim)

mu_init_hl_init = [mean(theta_hl_init,1), mean(eta_hl_init), mean(eps_hl_init)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1c.:
% get mit approximation of the conditional joint density of
% parameters and future returns given the returns are below VaR_prelim
% approximation of the joint high loss distribution
% here: not future returns but future disturbances  (varepsilons)

draw_hl_init = [theta_hl_init, eta_hl_init, eps_hl_init];
lnk_hl_init = kernel(theta_hl_init) + prior_const(1,1) - 0.5*(eta_hl_init).^2 + duvt(eps_hl_init, theta_hl_init(:,4), hp, true);

w_hl_init = lnk_hl_init/sum(lnk_hl_init);
[mu_hl_init, Sigma_hl_init] = fn_muSigma(draw_hl_init, w_hl_init);

mit_hl_init.mu = mu_hl_init;
mit_hl_init.Sigma = Sigma_hl_init;
mit_hl_init.p = 1;
mit_hl_init.df = 5;

cont.DUPA = true;

% kernel = @(a) prior_sv(a,true,prior_const);
mit = mit_hl_init;
kernel = @(a) posterior_svt_hl(y, a, VaR_prelim, par_NAIS_init, prior_const, cont.nais); 
% [d, x, lng_y, lnw_x, par_NAIS] = posterior_sv_hl(y, theta, VaR, par_NAIS_init, prior_const, cont) 

% clear GamMat
% save('results/sv_mitisem_hl_n.mat');
% x_gam = (0:0.00001:100)'+0.00001;
% GamMat = gamma(x_gam);
% load('sv_mitisem_hl.mat', 'theta', 'lnk', 'x', 'lng_y', 'lnw_x', 'par_NAIS');

[theta_hl, lnk_hl, ~, x_hl, lng_y_hl, lnw_x_hl, par_NAIS_hl] = fn_rmvgt_robust(N, mit, kernel, cont.DUPA);

lnd_hl = dmvgt(theta_hl, mit_hl_init, true, GamMat);
prior_hl = lnk_hl - lng_y_hl - lnw_x_hl;

w_nais = fn_ISwgts(prior_hl, lnd_hl, true);
b_wmean_hl = sum(repmat(w_nais',T,1).*par_NAIS_1.b,2);
C_wmean_hl = sum(repmat(w_nais',T,1).*par_NAIS_1.C,2);
b_mean_hl = mean(par_NAIS_hl.b,2);
C_mean_hl = mean(par_NAIS_hl.C,2);

load('SMLt_ibm_smooth.mat', 'theta_smooth', 'V_smooth');
x_mean_hl = mean(x_hl,1);
x_mean_1 = mean(x_1,1);
x_median_hl = median(x_hl,1);
x_median_1 = median(x_1,1);

if plot_on
    figure(4)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
	set(gcf,'defaulttextinterpreter','latex');
    
    subplot(2,2,1)
    [f,xi] = ksdensity(theta_hl(:,1));
    hold on    
    histnorm(theta_hl(:,1),20)   
    h = findobj(gca, 'Type','patch');
    set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
    plot(xi,f,'r');
    hold off
    title('$$c$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(2,2,2)
    [f,xi] = ksdensity(theta_hl(:,2));
    hold on       
    histnorm(theta_hl(:,2),20)   
    h = findobj(gca, 'Type','patch');
    set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
    plot(xi,f,'r');
    hold off
    title('$$\phi$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(2,2,3)
    [f,xi] = ksdensity(theta_hl(:,3));
    hold on       
    histnorm(theta_hl(:,3),20)
    h = findobj(gca, 'Type','patch');
    set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
    plot(xi,f,'r');
    hold off
    title('$$\sigma^{2}_{\eta}$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(2,2,4)
    [f,xi] = ksdensity(theta_hl(:,4));
    hold on       
    histnorm(theta_hl(:,4),20)
    h = findobj(gca, 'Type','patch');
    set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')    
    plot(xi,f,'r');
    hold off
    title('$$\nu$$')
	set(gca,'TickLabelInterpreter','latex')
    
%     suptitle('Empirical approximate joint posterior density $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)y$$.' )
    if print_on
		name = 'figures/svt_draws_param_hl.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(5)
    set(gcf,'units','normalized','outerposition',[0 0 1 0.33]);
	set(gcf,'defaulttextinterpreter','latex');
    
    subplot(1,2,1)
    [f,xi] = ksdensity(theta_hl(:,5));
    hold on    
    histnorm(theta_hl(:,5),20)
    h = findobj(gca, 'Type','patch');
    set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
    plot(xi,f,'r');
    hold off
    title('$$\eta_{T+1}$$')
	set(gca,'TickLabelInterpreter','latex')
    
    subplot(1,2,2)
    [f,xi] = ksdensity(theta_hl(:,6));
    hold on       
    histnorm(theta_hl(:,6),20)
    h = findobj(gca, 'Type','patch');
    set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
    plot(xi,f,'r');
    hold off
    title('$$\varepsilon_{T+1}$$')
	set(gca,'TickLabelInterpreter','latex')
    
%     suptitle('Empirical approximate joint posterior density $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)y$$.' )
    if print_on
		name = 'figures/svt_draws_err_hl.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
    end   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(6)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
	set(gcf,'defaulttextinterpreter','latex');
    
    subplot(2,1,1)
    hold on
    plot(b_wmean_hl)    
%     plot(b_mean_hl) 
    hold off  
    title('$$b$$')
    set(gca,'TickLabelInterpreter','latex')

    subplot(2,1,2)
    hold on
    plot(C_wmean_hl)    
%     plot(C_mean_hl) 
    hold off  
    title('$$C$$')
	set(gca,'TickLabelInterpreter','latex')
    
%     suptitle('Average NAIS parameters given the draws from $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)y$$.')

	if print_on
		name = 'figures/svt_NAIS_param_hl.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(7)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
	set(gcf,'defaulttextinterpreter','latex');
    
    subplot(2,1,1)

    hold on
    plot(b_mean_1) 
    plot(b_mean_hl,'r')
%     plot(par_NAIS_1.b(:,1))
%     plot(par_NAIS_hl.b(:,1),'r')
    hold off
    title('$$b$$')
    set(gca,'TickLabelInterpreter','latex')

    subplot(2,1,2)
    hold on
    plot(C_mean_1)
    plot(C_mean_hl,'r')

%     plot(par_NAIS_1.C(:,1))
%     plot(par_NAIS_hl.C(:,1),'r')
    hold off
    title('$$C$$')
	set(gca,'TickLabelInterpreter','latex')
    
%     suptitle('Average NAIS parameters given the draws from $$q_{1,\zeta}(\theta|y)$$ and $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)$$.')

	if print_on
		name = 'figures/svt_NAIS_param_both.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
    end
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(8)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
    hold on
    plot(x_median_1) 
    plot(x_median_hl,'r')
    plot(theta_smooth,'k','LineWidth',1)
    hold off
    title('Median signal path given the draws from $$q_{1,\zeta}(\theta|y)$$ and $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)$$.')
    set(gca,'TickLabelInterpreter','latex')

	if print_on
		name = 'figures/svt_median_x_both.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
    end
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(9)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
    hold on
    plot(x_mean_1) 
    plot(x_mean_hl,'r')
    plot(theta_smooth,'k','LineWidth',1)
    hold off
    title('Average signal path given the draws from $$q_{1,\zeta}(\theta|y)$$ and $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)$$.')
    set(gca,'TickLabelInterpreter','latex')

	if print_on
		name = 'figures/svt_mean_x_both.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 2:
% use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
% to estiamte VaR and ES for theta and y (or alpha in eps)

theta_opt = [theta1,eta_h1,eps_h1;theta_hl];
x_opt = [x_1;x_hl];
par_NAIS_opt.b = [par_NAIS_1.b , par_NAIS_hl.b];
par_NAIS_opt.C = [par_NAIS_1.C , par_NAIS_hl.C];

kernel = @(aa,bb,cc) posterior_svt_whole(y, aa, bb, prior_const, cont.nais, cc);
% [d, x, lng_y, lnw_x, par_NAIS] = posterior_sv_whole(y, theta, par_NAIS_init, prior_const, cont, x_in) 
[lnk_opt, x_opt, lng_y_opt, lnw_x_opt, par_NAIS_opt] = kernel(theta_opt,par_NAIS_opt,x_opt);

%% IS weights
exp_lnd1 = 0.5*normpdf(theta_opt(:,5)).*duvt(theta_opt(:,6), theta_opt(:,4), 1, false).*dmvgt(theta_opt(:,1:4),mit_init,false, GamMat);
exp_lnd2 = 0.5*dmvgt(theta_opt,mit_hl_init,false, GamMat);
exp_lnd = exp_lnd1 + exp_lnd2;
lnd_opt = log(exp_lnd);
w_opt =  fn_ISwgts(lnk_opt-lnw_x_opt, lnd_opt, false);


c_opt = theta_opt(:,1);
phi_opt = theta_opt(:,2);
sigma2_opt = theta_opt(:,3);
nu_opt = theta_opt(:,4);
rho_opt = (nu_opt-2)./nu_opt;
eta_opt = theta_opt(:,5);
eps_opt = theta_opt(:,6);

x_opt_h1 = c_opt + phi_opt.*(x_opt(:,end) - c_opt) + sqrt(sigma2_opt).*eta_opt;
y_opt_h1 = sqrt(rho_opt).*exp(0.5*x_opt_h1).*eps_opt;
[PL_opt_h1, ind] = sort(fn_PL(y_opt_h1));
ES_estim = sum(PL_opt_h1(PL_opt_h1<VaR_prelim))/(2*N*p_bar);
VaR_estim = PL_opt_h1(floor(p_bar*2*N));

if plot_on    
    figure(10)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
    hold on
%     plot(PL_opt_h1)
    plot(PL_opt_h1((PL_opt_h1>-20) & (PL_opt_h1<15)))

    hold off
%     subplot(3,1,1)
%     subplot(3,1,2)
%     subplot(3,1,3)
    title('Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$.')
    set(gca,'TickLabelInterpreter','latex')

	if print_on
		name = 'figures/svt_predict.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end

clear b_mean_1 b_mean_hl b_pmean_1 b_wmean_1 b_wmean_hl 
clear C_mean_1 C_mean_hl C_pmean_1 C_wmean_1 C_wmean_hl 
clear c1 c_opt phi1 phi_opt sigma21 sigma2_opt nu1 rho1 nu_opt rho_opt
clear eta_opt eps_opt eta_h1 eps_h1
clear lnd1 lnd_hl lnd_opt
clear PL_h1 PL_opt PL_opt_h1
clear theta_smooth V_smooth
clear GamMat x_gam fig h f xi xx

save(['results/svt_mitisem.mat']);