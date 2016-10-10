%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 'svt'; % 'svt'
algo = 'MitISEM';

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

cont2 = EMitISEM_Control(model);
cont2.mit.dfnc = 15;
cont2.df.range = [5,20];

mit_hl_init.df = cont2.mit.dfnc;

VaR_mit = zeros(N_sim,1);
ES_mit = zeros(N_sim,1);
RNE_mit = zeros(N_sim,1);
time_mit = zeros(2,1);

%% PRELIM
name =  [results_path,model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name,'VaR_prelim','mit1','mit_hl_init');
    
if strcmp(model,'sv')
    kernel = @(aa) posterior_sv_hl(y, aa, mean(VaR_prelim), par_NAIS_init, prior_const, cont2.nais); 
else
    kernel = @(aa) posterior_svt_hl(y, aa, mean(VaR_prelim), par_NAIS_init, prior_const, cont2.nais); 
end
cont = cont2
mit_init = mit_hl_init
[mit2, theta2, x2, w2, lnk2, lng_y2, lnw_x2, CV2] = EMitISEM(mit_hl_init, kernel, cont2, GamMat);
mit2 = mit_new 

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name, 'cont2', 'mit2', 'theta2', 'x2', 'w2', 'lnk2', 'lng_y2', 'lnw_x2', 'CV2')
end
%  
% mit2 = mit_new;
% theta2 = theta;
% x2 = x;
% w2 = w_norm;
% lnk2 = lnk;
% lng_y2 = lng_y;
% lnw_x2 = lnw_x;
% CV2 = CV;

%% Generate set of draws of theta using independence MH with naive candiate
% (based on the SML estimates) 
% simulate returns based on the draw of theta and future return paths     
% compute VaR_prelim


tic
for sim = 1:N_sim    
    fprintf('VaR IS sim = %i.\n', sim);
        
    theta1 = rmvgt2(M/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
    eta_h1_1 = randn(M/2,H);
    if strcmp(model,'sv')
        eps_h1_1 = randn(M/2,H);
    else
        nu1 = theta1(:,4);
        rho1 = (nu1-2)./nu1;
        eps_h1_1 = trnd(repmat(nu1,1,H));
    end      
    theta1 = [theta1, eta_h1_1, eps_h1_1];
    theta2 = rmvgt2(M/2, mit2.mu, mit2.Sigma, mit2.df, mit2.p); 
    theta_mit = [theta1; theta2];
    
    if strcmp(model,'sv')
        kernel = @(aa) posterior_sv(y, aa, par_NAIS_init, prior_const, cont2.nais);
    else
        kernel = @(aa) posterior_svt(y, aa, par_NAIS_init, prior_const, cont2.nais);
    end
    
    lnk = zeros(M,1);
    x = zeros(M,cont2.nais.HP+1);
    lng_y = zeros(M,1);
    lnw_x = zeros(M,1);
    eps_bar = zeros(M,1);
    eps_sim = zeros(M,1);
    C_sim = zeros(M,1);
    lnp_T = zeros(M,cont2.nais.HP+1);
    RND = zeros(M,1);
     
    for ii = 1:(M/1000)
        fprintf('ii = %i\n',ii)
        ind = (1:1000) + (ii-1)*1000; 
        [lnk(ind,:), x(ind,:), lng_y(ind,:), lnw_x(ind,:), eps_bar(ind,:), eps_sim(ind,:), C_sim(ind,:), lnp_T(ind,:), RND(ind,:)] = kernel(theta_mit(ind,:));      
    end
    
    if strcmp(model,'sv')
        lnk = lnk + 2*prior_const(1,1) - 0.5*(theta_mit(:,4)).^2 - 0.5*(theta_mit(:,5)).^2;
    else
        lnk = lnk + prior_const(1,1) - 0.5*(theta_mit(:,5)).^2 + duvt(theta_mit(:,6), theta_mit(:,4), 1, true); 
    end
  
    if strcmp(model,'sv')
        exp_lnd1 = 0.5*normpdf(theta_mit(:,4)).*normpdf(theta_mit(:,5)).*dmvgt(theta_mit(:,1:3), mit1, false, GamMat);
    else
        exp_lnd1 = 0.5*normpdf(theta_mit(:,5)).*duvt(theta_mit(:,6), theta_mit(:,4), 1, false).*dmvgt(theta_mit(:,1:4), mit1, false, GamMat);
    end
    exp_lnd2 = 0.5*dmvgt(theta_mit, mit2, false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);
       
    w_opt = fn_ISwgts(lnk, lnd_opt, false);    
    
    if strcmp(model,'sv')
        y_mit = predict_sv(theta_mit, x(:,end), H);
    else
        y_mit = predict_svt(theta_mit, x(:,end), H);       
    end
    
    ind_opt = (fn_PL(y_mit) <= mean(VaR_prelim));
    RNE_mit(sim,1) = fn_RNE(ind_opt, 'IS', w_opt);
    dens = struct('y',y_mit,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_mit(sim,1) = IS_estim(1,1);
    ES_mit(sim,1) = IS_estim(1,2);

    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_mit(sim,1), model, algo);  
end
time_mit(2,1) = toc/N_sim;
 
if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_mit','ES_mit','time_mit','RNE_mit','-append')
end