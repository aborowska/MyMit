%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 'sv'; % 'svt'
algo = 'Prelim';

old = false;
crisis = false;
recent = true;
if old 
    y = csvread('IBM_ret.csv');
    results_path = 'results/EQERMit/old/';
elseif crisis 
    y = csvread('GSPC_ret_updated.csv'); 
    results_path = 'results/EQERMit/crisis/';
elseif recent
    y = csvread('GSPC_ret_updated_short_sv.csv');
    results_path = 'results/EQERMit/recent/';
 end

y = 100*y;
T = size(y,1);
y_T = y(T);

M = 10000;
BurnIn = 1000;
N_sim = 20;
p_bar = 0.01;
H = 1;     % prediction horizon 

plot_on = false;
save_on = true;


%% Constants
if strcmp(model, 'sv')
    prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  2.5*log(0.025), ...
        -log(gamma(2.5))];
else
    hyper = 0.01; % 0.01 - VERY UNINFORMATIVE PRIOR
    prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  2.5*log(0.025), ...
        -log(gamma(2.5)), hyper]; % the last one is lambda for nu
end

% % Prior specification: (in prior_sv or prior_svt)
% logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
% logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
% % % logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;
% logpdf_invgamma = @(x) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(x) - 0.025./x;
% % % logpdf_chi2 = @(x) prior_const(1,3) + prior_const(1,4) -0.5*log(x) - 0.5*x;
% if strcmp(model, 'svt')
%     logpdf_exp = @(x) log(prior_const(1,5)) - prior_const(1,5)*(x - 2);  
% end

par_NAIS_init.b = zeros(T,1);
par_NAIS_init.C = ones(T,1); 

%% Initialisation
cont_prelim = EMitISEM_Control(model);
cont_prelim.mit.dfnc = 3;

VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
RNE_prelim = zeros(N_sim,1);
RNE_ES_prelim =  zeros(N_sim,1);
accept_prelim = zeros(N_sim,1);
time_prelim = zeros(2,1);
 
  
%% SML 
% load('SML_ibm.mat', 'par_SV_opt', 'heN_sim_SV_opt') 
% theta = [c, phi, sigma2, nu]

if old 
    if strcmp(model,'sv')
        load('SML_ibm.mat', 'par_SV_opt', 'V_SV_corr_opt') 
    else
        load('SMLt_ibm.mat', 'par_SV_opt', 'V_SV_corr_opt') 
    end   
elseif crisis 
    if strcmp(model,'sv')
        load('SML_gspc_updated_sv.mat', 'par_SV_opt', 'V_SV_corr_opt') 
     else
        load('SMLt_gspc_updated_sv.mat', 'par_SV_opt', 'V_SV_corr_opt') 
    end
elseif recent
    if strcmp(model,'sv')
        load('SML_gspc_updated_short_sv.mat', 'par_SV_opt', 'V_SV_corr_opt') 
     else
        load('SMLt_gspc_updated_short_sv.mat', 'par_SV_opt', 'V_SV_corr_opt') 
    end    
end

d = size(par_SV_opt,2);

% Sigma = inv(hess_sim_SV_opt);
Sigma = V_SV_corr_opt;
Sigma = reshape(Sigma,1,d^2);
mit_init.mu = par_SV_opt;
mit_init.Sigma = Sigma;
mit_init.p = cont_prelim.mit.pnc;
mit_init.df = cont_prelim.mit.dfnc;
    
if strcmp(model,'sv')
    kernel = @(aa) posterior_sv(y, aa, par_NAIS_init, prior_const, cont_prelim.nais);
else
    kernel = @(aa) posterior_svt(y, aa, par_NAIS_init, prior_const, cont_prelim.nais);
end

[mit1, theta1, x1, w1, lnk1, lng_y1, lnw_x1, CV1] = EMitISEM(mit_init, kernel, cont, GamMat);

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name, 'cont_prelim', 'mit1', 'theta1', 'x1', 'w1', 'lnk1', 'lng_y1', 'lnw_x1', 'CV1')
end
 

%% Generate set of draws of theta using independence MH with naive candiate
% (based on the SML estimates) 
% simulate returns based on the draw of theta and future return paths     
% compute VaR_prelim

tic
for sim = 1:N_sim    
    fprintf('\nPrelim sim = %i.\n', sim);
   
    [theta, x, lnw, lnk, lng_y, lnw_x, ~, accept_prelim(sim,1)] = EMit_MH(M+BurnIn, d, kernel, mit1, GamMat, true);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept_prelim(sim,1), model, algo);

    theta = theta(BurnIn+1:M+BurnIn,:);
    x = x(BurnIn+1:M+BurnIn,:);
    
    if strcmp(model,'sv')
        [y_prelim, eta_hp, eps_hp] = predict_sv(theta, x(:,end), H);
    else
        [y_prelim, eta_hp, eps_hp] = predict_svt(theta, x(:,end), H);       
    end
    
    ind_real = find(sum(imag(y_prelim),2)==0);
    M_real = length(ind_real); 
    fprintf('M_real = %i.\n',M_real)
    y_prelim = y_prelim(ind_real,:);
        
    PL_prelim_ind = fn_PL(y_prelim);
    [PL_prelim, ind] = sort(PL_prelim_ind);
    VaR_prelim(sim,1) = PL_prelim(round(p_bar*M_real));
    ES_prelim(sim,1) = mean(PL_prelim(round(1:p_bar*M_real)));   
     
    ind_prelim = double((PL_prelim_ind <= VaR_prelim(sim,1)));
    RNE_prelim(sim,1) = fn_RNE(ind_prelim, 'MH',[],'Q');
    ind_prelim = PL_prelim_ind((PL_prelim_ind <= VaR_prelim(sim,1)));
    RNE_ES_prelim(sim,1) = fn_RNE(ind_prelim, 'MH',[],'Q');
    
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
end
time_prelim(2,1) = toc/N_sim;
 
if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','accept_prelim','time_prelim','RNE_prelim','RNE_ES_prelim','-append')
end


%% QERMit 1c.:
% get mit approximation of the conditional joint density of
% parameters and future returns given the returns are below VaR_prelim

% take all the variables corresponding to the returns below VaR_prelim
% --> to constuct mit_init_hl, i.e. a starting mixture for mit2 
% where mu and Sigma are IS based

M_low = find(PL_prelim < mean(VaR_prelim), 1, 'last' );

eps_hp = eps_hp(ind_real,:);
eta_hp = eta_hp(ind_real,:);
theta = theta(ind_real,:);
lnw = lnw(ind_real,:);
lnk = lnk(ind_real,:);
    
eps_hl_init = eps_hp(ind,:); 
eps_hl_init = eps_hl_init(1:M_low,:);

eta_hl_init = eta_hp(ind,:); 
eta_hl_init = eta_hl_init(1:M_low,:);

theta_hl_init = theta(ind,:);
theta_hl_init = theta_hl_init(1:M_low,:);

lnw_hl_init = lnw(ind,:);
lnw_hl_init = lnw_hl_init(1:M_low,:);

lnk_hl_init = lnk(ind,:);
lnk_hl_init = lnk_hl_init(1:M_low,:);

draw_hl_init = [theta_hl_init, eta_hl_init, eps_hl_init];

% if strcmp(model,'sv')
%     lnk_hl_init = lnk_hl_init + 2*prior_const(1,1) - 0.5*(eta_hl_init).^2 - 0.5*(eps_hl_init).^2;
% else
%     lnk_hl_init = lnk_hl_init + prior_const(1,1) - 0.5*(eta_hl_init).^2 + duvt(eps_hl_init, theta_hl_init(:,4), hp, true);
% end

ind_fin = isfinite(lnk_hl_init);
lnk_hl_init = lnk_hl_init(ind_fin,:); 
draw_hl_init = draw_hl_init(ind_fin,:);
if isempty(strfind(model, 't'))
    lnd_hl_init = dmvgt(draw_hl_init(:,1:3), mit1, true, GamMat);
else
    lnd_hl_init = dmvgt(draw_hl_init(:,1:4), mit1, true, GamMat);
end

w_hl_init =  fn_ISwgts(lnk_hl_init, lnd_hl_init, false);
[mu_hl_init, Sigma_hl_init] = fn_muSigma(draw_hl_init, w_hl_init);

mit_hl_init.mu = mu_hl_init;
mit_hl_init.Sigma = Sigma_hl_init;
mit_hl_init.p = 1;
mit_hl_init.df = cont_prelim.mit.dfnc; %1
 
if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'mit_hl_init','-append')
end
