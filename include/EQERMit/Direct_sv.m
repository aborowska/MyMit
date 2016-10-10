%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 'svt'; % 'svt'
algo = 'Direct';

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
N_sim = 1;
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
cont_direct = EMitISEM_Control(model);
cont_direct.mit.dfnc = 3;

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
RNE_direct = zeros(N_sim,1);
RNE_ES_direct =  zeros(N_sim,1);
accept_direct = zeros(N_sim,1);
time_direct = zeros(2,1);
 
  
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
        load('SML_gspc_updated.mat', 'par_SV_opt', 'V_SV_corr_opt') 
     else
        load('SMLt_gspc_updated.mat', 'par_SV_opt', 'V_SV_corr_opt') 
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
mit_direct.mu = par_SV_opt;
mit_direct.Sigma = Sigma;
mit_direct.p = cont_direct.mit.pnc;
mit_direct.df = cont_direct.mit.dfnc;
    
if strcmp(model,'sv')
    kernel = @(aa) posterior_sv(y, aa, par_NAIS_init, prior_const, cont_direct.nais);
else
    kernel = @(aa) posterior_svt(y, aa, par_NAIS_init, prior_const, cont_direct.nais);
end

%% Generate set of draws of theta using independence MH with naive candiate
% (based on the SML estimates) 
% simulate returns based on the draw of theta and future return paths     
% compute VaR_direct

tic
for sim = 1:N_sim    
    fprintf('\nDirect sim = %i.\n', sim);
        
    [theta, x, ~, ~, ~, ~, ~, accept_direct(sim,1)] = EMit_MH(M+BurnIn, d, kernel, mit_direct, GamMat, true);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept_direct(sim,1), model, algo);

    theta = theta(BurnIn+1:M+BurnIn,:);
    x = x(BurnIn+1:M+BurnIn,:);
    
    if strcmp(model,'sv')
        y_direct = predict_sv(theta, x(:,end), H);
    else
        y_direct = predict_svt(theta, x(:,end), H);       
    end
    
    ind_real = find(sum(imag(y_direct),2)==0);
    M_real = length(ind_real); 
    fprintf('M_real = %i.\n',M_real)
    y_direct = y_direct(ind_real,:);
        
    PL_direct_ind = fn_PL(y_direct);
    PL_direct = sort(PL_direct_ind);
    VaR_direct(sim,1) = PL_direct(round(p_bar*M_real));
    ES_direct(sim,1) = mean(PL_direct(round(1:p_bar*M_real)));   
     
    ind_direct = double((PL_direct_ind <= VaR_direct(sim,1)));
    RNE_direct(sim,1) = fn_RNE(ind_direct, 'MH',[],'Q');
    ind_direct = PL_direct_ind((PL_direct_ind <= VaR_direct(sim,1)));
    RNE_ES_direct(sim,1) = fn_RNE(ind_direct, 'MH',[],'Q');
    
    fprintf('Direct 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_direct(sim,1), model, algo);
end
time_direct(2,1) = toc/N_sim;
 
if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','mit_direct','accept_direct','time_direct','RNE_direct','RNE_ES_direct')
end