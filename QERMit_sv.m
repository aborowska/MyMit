%% Initialisation
% clc
clear all
close all
% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s); 
addpath(genpath('include/'));

v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014b)')
    v_new = 1;
else
    v_new = 0;
end

plot_on = false;
print_on = false;
save_on = false;
plot_on2 = true;

% User prompts
usr_prompt = 0;

if usr_prompt
    model = 'Select model [sv/svt]: ';
    model = input(model,'s');

    N_sim = 'Select no. of MC replications: [5/10/20] ';
    N_sim = input(N_sim);
   
    hp = 'Select no. of days-ahead for VaR estimation: [1] ';
    hp = input(hp);

    p_bar = 'Select the quantile for VaR estimation: [0.01/0.02/0.05] ';
    p_bar = input(p_bar);
else
    model = 'sv';
    N_sim = 10;
    hp = 1;
    p_bar = 0.01;
end

%% Constants
x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

% prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -2.5*log(0.025), -log(gamma(2.5))];
prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  2.5*log(0.025), -log(gamma(2.5))];
% prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -0.5*log(2), -log(gamma(0.5))];

logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
% % logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;
logpdf_invgamma = @(x) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(x) - 0.025./x;
% % logpdf_chi2 = @(x) prior_const(1,3) + prior_const(1,4) -0.5*log(x) - 0.5*x;

%% Data
% y = csvread('IBM_ret.csv');
% y = 100*y;
% % fprintf('IBM logreturns kurtosis: %6.4f, and skewneN_sim: %6.4f.\n', kurtosis(data), skewneN_sim(data));
% % http://faculty.chicagobooth.edu/nicholas.polson/research/papers/jpr2.pdf

% QERMit ARCH
y = csvread('GSPC_ret.csv');
y = 100*y;
ind_arch = find(y<=-5.5, 1, 'last' );
y = y(1:ind_arch,1);
y = y - mean(y);


T = size(y,1);

% SV_plot0;

%% Initialisation
EMitISEM_Control
cont.resmpl_on = false;

if usr_prompt
    N = 'Select number of draws: [1000/2000/5000/10000] ';
    N = input(N);
    cont.mit.N = N;
else
    cont.mit.N = 10000;
    N = cont.mit.N;
end


par_NAIS_init.b = zeros(T,1);
par_NAIS_init.C = ones(T,1); 

%% QERMit 1a.:
% SML ... <<<<<<<<<<<<<<<<<<<<<<
% load('SML_ibm.mat', 'par_SV_opt', 'heN_sim_SV_opt') 
% theta = [c, phi, sigma2, nu]
if strcmp(model,'sv')
    mu_init = [0.5, 0.98, 0.15^2];
%     load('SML_ibm.mat', 'par_SV_opt', 'V_SV_corr_opt') 
    load('SML_arch.mat', 'par_SV_opt', 'V_SV_corr_opt') 
else
    mu_init = [0.5, 0.98, 0.15^2, 7];
    load('SMLt_ibm.mat', 'par_SV_opt', 'V_SV_corr_opt') 
end

d = size(mu_init,2);

% Sigma = inv(hess_sim_SV_opt);
Sigma = V_SV_corr_opt;
Sigma = reshape(Sigma,1,d^2);
mit_init.mu = par_SV_opt;
mit_init.Sigma = Sigma;
mit_init.p = cont.mit.pnc;
mit_init.df = cont.mit.dfnc;

    
if strcmp(model,'sv')
    kernel_prior = @(a) prior_sv(a, prior_const); 
    kernel = @(a) posterior_sv(y, a, par_NAIS_init, prior_const, cont.nais);
else
    kernel_prior = @(a) prior_svt(a, prior_const); 
    kernel = @(a) posterior_svt(y, a, par_NAIS_init, prior_const, cont.nais);
end
[mit1, theta1, x1, w1, lnk1, lng_y1, lnw_x1, CV1] = EMitISEM(y, mit_init, kernel_prior, kernel, cont, GamMat);

        
c1 = theta1(:,1);
phi1 = theta1(:,2);
sigma21 = theta1(:,3);

eta_h1_1 = randn(N,1);

if strcmp(model,'sv')
    eps_h1_1 = randn(N,1);
    x_h1_1 = c1 + phi1.*(x1(:,end) - c1) + sqrt(sigma21).*eta_h1_1;
    y_h1_1 = exp(0.5*x_h1_1).*eps_h1_1;
else
    nu1 = theta1(:,4);
    rho1 = (nu1-2)./nu1;
    eps_h1_1 = trnd(repmat(nu1,1,hp));

    x_h1_1 = c1 + phi1.*(x1(:,end) - c1) + sqrt(sigma21).*eta_h1_1;
    y_h1_1 = sqrt(rho1).*exp(0.5*x_h1_1).*eps_h1_1;
end
PL_h1_1 = fn_PL(y_h1_1);
    
SV_plot1;

if save_on
    if strcmp(model,'sv')
        save(['results/sv_mit1_',int2str(N),'.mat'],'mit1','theta1', 'x1', 'w1', 'lnk1','CV1','cont','p_bar','N');
    else
        save(['results/svt_mit1_',int2str(N),'.mat'],'theta1', 'x1', 'w1', 'lnk1','CV1','cont','p_bar','N');
    end
end

% M = 10*N;
M = 4000;

%% QERMit 1b.:
% generate set of draws of theta using independence MH with candiate from MitISEM; 
% simulate returns based on the draw of theta and future return paths     
% compute VaR_prelim
N_sim = 20;
VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
accept = zeros(N_sim,1);
for sim = 1:N_sim    
    [theta, x, lnw, lnk, ~, ~, ~, accept(sim,1)] = EMit_MH(y, M+1000, kernel_prior, kernel, mit1, GamMat, true);
    fprintf('(MitISEM) MH acceptance rate: %6.4f. \n',accept(sim,1));

    theta = theta(1001:M+1000,:);
    x = x(1001:M+1000,:);
    lnw = lnw(1001:M+1000,:);
    lnk = lnk(1001:M+1000,:);

    if (sim == N_sim)
        SV_autocorr;
    end
    % High loss, 1 days horizon
    % approximate the high loss distribution of (theta,eps*,eta*) where 
    % eps*={eps_T+1,...,eps_T+hp}
    % eta*={eta_T+1,...,eta_T+hp}

    c = theta(:,1);
    phi = theta(:,2);
    sigma2 = theta(:,3);
    
    eta_h1 = randn(M,1);
    
    if strcmp(model,'sv')
        eps_h1 = randn(M,1);
        x_h1 = c + phi.*(x(:,end) - c) + sqrt(sigma2).*eta_h1;
        y_h1 = exp(0.5*x_h1).*eps_h1;
    else
        nu = theta(:,4);
        rho = (nu-2)./nu;
        eps_h1 = trnd(repmat(nu,1,hp));

        x_h1 = c + phi.*(x(:,end) - c) + sqrt(sigma2).*eta_h1;
        y_h1 = sqrt(rho).*exp(0.5*x_h1).*eps_h1;
    end
    
    ind_real = (imag(y_h1)==0);
    M_real = sum(ind_real); 
    y_h1 = y_h1(ind_real,:);
    eps_h1 = eps_h1(ind_real,:);
    eta_h1 = eta_h1(ind_real,:);
    theta = theta(ind_real,:);
    lnw = lnw(ind_real,:);
    lnk = lnk (ind_real,:);
        
    [PL_h1, ind] = sort(fn_PL(y_h1));
   
    VaR_prelim(sim,1) = PL_h1(round(p_bar*M_real));
    ES_prelim(sim,1) = mean(PL_h1(round(1:p_bar*M)));   
    fprintf('p_bar = %4.2f, VaR_prelim = %4.5f. \n', p_bar, VaR_prelim(sim,1))
end

VaR_prelim_MC =  VaR_prelim;
VaR_prelim = mean(VaR_prelim_MC);
M_real = max(find(PL_h1 < VaR_prelim)); 

%% QERMit 1c.:
% get mit approximation of the conditional joint density of
% parameters and future returns given the returns are below VaR_prelim
% approximation of the joint high loss distribution
% here: not future returns but future disturbances  (varepsilons)
eps_hl_init = eps_h1(ind,:); 
eps_hl_init = eps_hl_init(1:M_real,:);

eta_hl_init = eta_h1(ind,:); 
eta_hl_init = eta_hl_init(1:M_real,:);

theta_hl_init = theta(ind,:);
theta_hl_init = theta_hl_init(1:M_real,:);

lnw_hl_init = lnw(ind,:);
lnw_hl_init = lnw_hl_init(1:M_real,:);

lnk_hl_init = lnk(ind,:);
lnk_hl_init = lnk_hl_init(1:M_real,:);

draw_hl_init = [theta_hl_init, eta_hl_init, eps_hl_init];


% >>>>>>>>>>> ?????
if strcmp(model,'sv')
    lnk_hl_init = lnk_hl_init + 2*prior_const(1,1) - 0.5*(eta_hl_init).^2 - 0.5*(eps_hl_init).^2;
else
    lnk_hl_init = lnk_hl_init + prior_const(1,1) - 0.5*(eta_hl_init).^2 + duvt(eps_hl_init, theta_hl_init(:,4), hp, true);
end
% >>>>>>>>>>> ?????
w_hl_init = lnk_hl_init/sum(lnk_hl_init);
[mu_hl_init, Sigma_hl_init] = fn_muSigma(draw_hl_init, w_hl_init);

mit_hl_init.mu = mu_hl_init;
mit_hl_init.Sigma = Sigma_hl_init;
mit_hl_init.p = 1;
mit_hl_init.df = 1;

if strcmp(model,'sv')
    kernel_prior = @(a) prior_sv_hl(a, prior_const); 
    kernel = @(a) posterior_sv_hl(y, a, VaR_prelim, par_NAIS_init, prior_const, cont.nais); 
else
    kernel_prior = @(a) prior_svt_hl(a, prior_const); 
    kernel = @(a) posterior_svt_hl(y, a, VaR_prelim, par_NAIS_init, prior_const, cont.nais); 
end

cont2 = cont;
cont2.mit.Hmax = 2;
% cont.mit.CV_tol = 0.1;
[mit2, theta2, x2, w2, lnk2, lng_y2, lnw_x2, CV2] = EMitISEM(y, mit_hl_init, kernel_prior, kernel, cont2, GamMat);

if save_on
    if strcmp(model,'sv')
        save(['results/sv_mit2_',int2str(N),'.mat'],'mit2','theta2','x2','w2','lnk2','lng_y2','lnw_x2','CV2','cont','p_bar','N','VaR_prelim');
    else
        save(['results/svt_mit2_',int2str(N),'.mat'],'mit2','theta2','x2','w2','lnk2','lng_y2','lnw_x2','CV2','cont','p_bar','N','VaR_prelim');
    end
end

c2 = theta2(:,1);
phi2 = theta2(:,2);
sigma22 = theta2(:,3);

if strcmp(model,'sv')
    eta_h1_2 = theta2(:,4);
    eps_h1_2 = theta2(:,5);
    x_h1_2 = c2 + phi2.*(x2(:,end) - c2) + sqrt(sigma22).*eta_h1_2;
    y_h1_2 = exp(0.5*x_h1_2).*eps_h1_2;
else
    nu2 = theta2(:,4);
	rho2 = (nu2-2)./nu2;
    eta_h1_2 = theta2(:,5);
    eps_h1_2 = theta2(:,6);

    x_h1_2 = c2 + phi2.*(x2(:,end) - c2) + sqrt(sigma22).*eta_h1_2;
    y_h1_2 = sqrt(rho2).*exp(0.5*x_h1_2).*eps_h1_2;   
end

PL_h1_2 = fn_PL(y_h1_2);
 
SV_plot2;

%% QERMit 2:
% use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
% to estiamte VaR and ES for theta and y (or alpha in eps)


VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONTE CARLO VaR_IS and ES_IS (and their NSEs) ESTIMATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 2000;
M_real = zeros(N_sim,1);

for sim = 1:N_sim

    fprintf('\n')
    fprintf('NSE sim = %i.\n', sim);
    fprintf('\n')

    %%% MIT1 %%%
    if strcmp(model,'sv')
        kernel_prior = @(a) prior_sv(a, prior_const); 
        kernel = @(a) posterior_sv(y, a, par_NAIS_init, prior_const, cont.nais);
    else
        kernel_prior = @(a) prior_svt(a, prior_const); 
        kernel = @(a) posterior_svt(y, a, par_NAIS_init, prior_const, cont.nais);
    end
    [theta1, x1, lnw1, lnk1, lng_y1, lnw_x1, ~, ~] = EMit_MH(y, M, kernel_prior, kernel, mit1, GamMat, false);

    eta_h1_1 = randn(M,1);

    if strcmp(model,'sv')
        eps_h1_1 = randn(M,1);
    else
        nu1 = theta1(:,4);
        rho1 = (nu1-2)./nu1;
        eps_h1_1 = trnd(repmat(nu1,1,hp));
    end    

    %%% MIT2 %%% 
    if strcmp(model,'sv')
        kernel_prior = @(a) prior_sv_hl(a, prior_const); 
        kernel = @(a) posterior_sv_hl(y, a, VaR_prelim, par_NAIS_init, prior_const, cont.nais); 
    else
        kernel_prior = @(a) prior_svt_hl(a, prior_const); 
        kernel = @(a) posterior_svt_hl(y, a, VaR_prelim, par_NAIS_init, prior_const, cont.nais); 
    end

    [theta2, x2, lnw2, lnk2, lng_y2, lnw_x2, ~, ~] = EMit_MH(y, M, kernel_prior, kernel, mit2, GamMat, false);


    %%% OPT %%%
    theta_opt = [theta1, eta_h1_1, eps_h1_1; theta2];
    % x_opt = [x1;x2];
    x_opt_end = [x1(:,end); x2(:,end)];
    % clear x1 x2
    % >>>>>>>>>>> ?????
    if strcmp(model,'sv')
        lnk_opt = lnk1 + 2*prior_const(1,1) - 0.5*(eta_h1_1).^2 - 0.5*(eps_h1_1).^2;
    else
        lnk_opt = lnk1 + prior_const(1,1) - 0.5*(eta_h1_1).^2 + duvt(eps_h1_1, theta1(:,4), hp, true);
    end

    lnk_opt = [lnk_opt; lnk2];
%     lng_y_opt = [lng_y1; lng_y2];
%     lnw_x_opt = [lnw_x1; lnw_x2];

    %% IS weights
    if strcmp(model,'sv')
        exp_lnd1 = 0.5*normpdf(theta_opt(:,4)).*normpdf(theta_opt(:,5)).*dmvgt(theta_opt(:,1:3), mit1, false, GamMat);
    else
        exp_lnd1 = 0.5*normpdf(theta_opt(:,5)).*duvt(theta_opt(:,6), theta_opt(:,4), 1, false).*dmvgt(theta_opt(:,1:4), mit1, false, GamMat);
    end
    exp_lnd2 = 0.5*dmvgt(theta_opt, mit2, false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);

    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    c_opt = theta_opt(:,1);
    phi_opt = theta_opt(:,2);
    sigma2_opt = theta_opt(:,3);
 
    if strcmp(model,'sv')
        eta_opt = theta_opt(:,4);
        eps_opt = theta_opt(:,5);  
        x_opt_h1 = c_opt + phi_opt.*(x_opt_end - c_opt) + sqrt(sigma2_opt).*eta_opt;
        y_opt_h1 = exp(0.5*x_opt_h1).*eps_opt;
    else
        nu_opt = theta_opt(:,4);
        rho_opt = (nu_opt-2)./nu_opt;
        eta_opt = theta_opt(:,5);
        eps_opt =  theta_opt(:,6);
        x_opt_h1 = c_opt + phi_opt.*(x_opt_end - c_opt) + sqrt(sigma2_opt).*eta_opt;
        y_opt_h1 = sqrt(rho_opt).*exp(0.5*x_opt_h1).*eps_opt;
    end

%     ind_real = (imag(y_opt_h1)==0);
%     M_real(sim,1) = sum(ind_real); 

    dens = struct('y',y_opt_h1,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);
    
    fprintf('(%s) IS 100*%4.2f%% VAR estimate: %6.4f. \n', model, p_bar, VaR_IS(sim,1));
    fprintf('(%s) IS 100*%4.2f%% ES estimate: %6.4f. \n', model, p_bar, ES_IS(sim,1));  

end

SV_plot3;


mean_VaR_prelim = mean(VaR_prelim_MC);
mean_ES_prelim = mean(ES_prelim);

NSE_VaR_prelim = std(VaR_prelim_MC);
NSE_ES_prelim = std(ES_prelim);

mean_VaR_IS = mean(VaR_IS);
mean_ES_IS = mean(ES_IS);

NSE_VaR_IS = std(VaR_IS);
NSE_ES_IS = std(ES_IS);

fprintf('(%s) 100*%4.2f%% VaR prelim (mean) estimate: %6.4f. \n', model, p_bar, mean_VaR_prelim);
fprintf('(%s) NSE VaR prelim: %6.4f. \n', model, NSE_VaR_prelim);
fprintf('(%s) VaR prelim: [%6.4f, %6.4f]. \n \n', model, mean_VaR_prelim - NSE_VaR_prelim, mean_VaR_prelim + NSE_VaR_prelim);

fprintf('(%s) 100*%4.2f%% VaR IS (mean) estimate: %6.4f. \n', model, p_bar, mean_VaR_IS);
fprintf('(%s) NSE VaR IS estimate: %6.4f. \n', model, NSE_VaR_IS);
fprintf('(%s) VaR: [%6.4f, %6.4f]. \n \n', model, mean_VaR_IS - NSE_VaR_IS, mean_VaR_IS + NSE_VaR_IS);

fprintf('(%s) 100*%4.2f%% ES prelim (mean) estimate: %6.4f. \n', model, p_bar, mean_ES_prelim);
fprintf('(%s) NSE ES prelim: %6.4f. \n', model, NSE_ES_prelim);
fprintf('(%s) ES prelim: [%6.4f, %6.4f]. \n \n', model, mean_ES_prelim - NSE_ES_prelim, mean_ES_prelim + NSE_ES_prelim);

fprintf('(%s) 100*%4.2f%% ES IS (mean) estimate: %6.4f. \n', model, p_bar, mean_ES_IS);
fprintf('(%s) NSE ES IS estimate: %6.4f. \n', model, NSE_ES_IS);
fprintf('(%s) ES: [%6.4f, %6.4f]. \n', model, mean_ES_IS - NSE_ES_IS, mean_ES_IS + NSE_ES_IS);


if plot_on2
    figure(590+100*p_bar)
%         set(gcf, 'visible', 'off');
%         set(gcf,'defaulttextinterpreter','latex');
    boxplot([VaR_prelim_MC, VaR_IS],'labels',{'VaR_prelim MC','VaR_IS'})        
    title(['(',model,' M = ',num2str(M),') ','100*', num2str(p_bar),'% VaR estimates: prelim and IS.'])
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    if print_on
        name = ['figures/(',model,')', num2str(p_bar),'_VaR_box_',num2str(M),'.png'];
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%

    figure(5900+100*p_bar)
%         set(gcf, 'visible', 'off');
%         set(gcf,'defaulttextinterpreter','latex');
    hold on; 
    bar(VaR_IS,'FaceColor',[0 0.4470 0.7410], 'EdgeColor','w'); 
    plot(0:(N_sim+1), (mean_VaR_prelim - NSE_VaR_prelim)*ones(N_sim+2,1),'r--'); 
    plot(0:(N_sim+1), (mean_VaR_prelim + NSE_VaR_prelim)*ones(N_sim+2,1),'r--'); 
    plot(0:(N_sim+1), mean_VaR_prelim*ones(N_sim+2,1),'r'); 
    hold off;
    title(['(',model,' M = ',num2str(M),') ','100*', num2str(p_bar),'% VaR IS estimates and the mean VaR prelim (+/- NSE VaR prelim).'])

    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    if print_on
        name = ['figures/(',model,')', num2str(p_bar),'_VaR_bar_',num2str(M),'.png'];
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
end
if save_on
    gen_out2;
end


%% Finish
% clear c1 c2 c_opt phi1 phi2 phi_opt sigma21 sigma22 sigma2_opt
% clear eta_opt eps_opt eta_hl eps_hl
% clear lnd1 lnd_hl lnd_opt
% clear PL_h1 PL_h1_hl PL_opt PL_opt_h1
% clear theta_smooth V_smooth
clear GamMat x_gam fig h f xi xx
% clear theta_opt x_opt
% clear exp_lnd1 exp_lnd2 exp_lnd lnd_opt w_opt x_opt_h1 y_opt_h1

% if save_on
%     if strcmp(model,'sv')
%         save(['results/sv_mitisem_',int2str(N),'.mat']);
%     else
%         save(['results/svt_mitisem_',int2str(N),'.mat']);
%     end
% end
