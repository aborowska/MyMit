%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

algo = 'PMitISEM';
model = 't_garch2_noS';

y = csvread('GSPC_ret_tgarch.csv');
y = 100*y;
% y = csvread('GSPC_ret.csv');
% y = y - mean(y);
data = y;

T = size(data,1);
y_T = data(T);
S = var(data);

p_bar = 0.01;
H = 100; % forecast horizon

M = 10000;
BurnIn = 1000;
N_sim = 20;
sim = 1;
% L = true;
% hyper = 1;
% theta = [omega, alpha, beta, mu, nu]
% mu_init = [0.008, 0.07, 0.9, 0.01, 10];
% % mu_init = [0.065 0.93 0.048 8.4];
mu_init = [0.009, 0.07, 0.9, 0.05, 11];
DD = size(mu_init,2);

plot_on = true;
save_on = false;

% Control parameters for MitISEM (cont) and PMitiISEM (cont2)
cont2 = MitISEM_Control;
cont2.mit.dfnc = 5;

VaR_pmit = zeros(N_sim,1);
ES_pmit = zeros(N_sim,1);
RNE_pmit = zeros(N_sim,1);
time_pmit = zeros(2,1);

%% PRELIM & BIG DRAW
name =  ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name);

% WEIGHTS to initialise PMitISEM
% future disturbances are generated from the target thus have weights 1
% log kernel evaluation - only for the parameter draws, epsilons are drawn
% from the target so their weigths are 1
% % kernel = @(xx) posterior_t_garch_noS_mex(xx, data, S, GamMat);
hyper = 0.1; 
kernel = @(a) posterior_t_garch_noS_hyper_mex(a, data, S, GamMat, hyper);
lnk_hl = kernel(draw_hl(:,1:DD)); 
% log candidate evaluation
lnd_hl = dmvgt(draw_hl(:,1:DD), mit1, true, GamMat);

% importance weights
w_hl = lnk_hl - lnd_hl;
w_hl = exp(w_hl - max(w_hl));

%% PMitISEM
partition = [1,DD+2:H+DD];
d = H+DD;
% S = var(data);

fn_const_X = @(xx) t_garch_noS_const_X2(xx, data, S);
fn_input_X = @(xx) t_garch_noS_input_X(xx, data, S);
% kernel = @(xx) posterior_t_garch_hl_noS_mex(xx, data, S, mean(VaR_prelim), GamMat);
kernel = @(a) posterior_t_garch_hl_noS_hyper_mex(a, data, S, mean(VaR_prelim), GamMat, hyper);

CV_old = cont2.mit.CV_old;
CV_tol = cont2.mit.CV_tol;

draw0 = draw_hl;
w0 = w_hl;
lnk0 = lnk_hl; %kernel(draw0);

cont2.mit.iter_max = 1;
if (H==250)
    cont2.mit.Hmax = 1;
end
cont2.df.range = [1,10]; %<<<==== was: [1;10]
cont = cont2;

tic
% [pmit, CV_mix, CV, iter, pmit_pre, pmit_pre2, pmit_adapt]  = PMitISEM(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont2, GamMat);
[pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM2(draw0, lnk0, w0, kernel, fn_const_X, fn_input_X, partition, d, cont2, GamMat);
time_pmit(1,1) = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','pmit','CV_mix','CV','iter')
end

%% VaR with PMit
tic
for sim = 1:N_sim   
    fprintf('\nVaR IS iter: %d\n',sim)

    theta1 = rmvgt2(M/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
    eps1 = zeros(M/2, H);
    for hh = 1:H
        eps1(:,hh) = trnd(theta1(:,DD)); % ERRORS ARE iid T!!
    end
    draw1 = [theta1, eps1];
    clear theta1 eps1
    input_X_1 = fn_input_X(draw1);

%         draw_pmit  = fn_p_rmvgt(M/2, pmit, d, partition, [], fn_const_X);  
    [draw_pmit, ~, input_X_pmit] = fn_p_rmvgt2(M/2, pmit, d, partition, [], fn_const_X, fn_input_X);         

    draw_opt = [draw1; draw_pmit];

    input_X.theta = draw_opt;
    input_X.h_T = [input_X_1.h_T; input_X_pmit.h_T];
    clear draw1 %draw_pmit
    clear input_X_1 input_X_pmit

%     kernel = @(xx) posterior_t_garch_noS_mex(xx, data ,S, GamMat);
    kernel = @(a) posterior_t_garch_noS_hyper_mex(a, data, S, GamMat, hyper);
    lnk_opt = kernel(draw_opt(:,1:5)); 

%     eps_pdf = zeros(M, 1);
%     for hh = 1:H
%         eps_pdf = eps_pdf + log(tpdf(draw_opt(:,DD+hh),draw_opt(:,DD)));
%     end
    eps_pdf = duvt(draw_opt(:,DD+1:H+DD), draw_opt(:,DD), H, true);
    lnk_opt = lnk_opt + eps_pdf;

    % optimal weights
    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,1:DD), mit1, true, GamMat));
%         exp_lnd2 = fn_dpmit(draw_opt, pmit, partition, fn_const_X, true, GamMat);
    exp_lnd2 = fn_dpmit2(input_X, pmit, partition, fn_const_X, true, GamMat);        

    exp_lnd2 = 0.5*exp(exp_lnd2);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % IS VaR estimation
%         h_T = volatility_t_garch_noS_mex(draw_opt(:,1:DD), data, S);
    h_T = input_X.h_T;     
    y_opt = predict_t_garch_noS(draw_opt(:,1:DD), y_T, h_T, H, draw_opt(:,DD+1:H+DD));   
%     ind_opt = (fn_PL(y_opt) <= mean(VaR_prelim));
%     RNE_pmit(sim,1) = fn_RNE(ind_opt, 'IS', w_opt);     
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_pmit(sim,1) = IS_estim(1,1);
    ES_pmit(sim,1) = IS_estim(1,2);   

    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_pmit(sim,1), model, algo);  
end
time_pmit(2,1) = toc/N_sim;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit','RNE_pmit')
end

h_pmit = volatility_t_garch_noS_mex(draw_pmit(:,1:DD), data, S);
y_pmit = predict_t_garch_noS(draw_pmit(:,1:DD), y_T, h_pmit, H, draw_pmit(:,DD+1:H+DD));
PL_pmit = fn_PL(y_pmit);
pmit_eff = sum(PL_pmit <= mean(VaR_prelim))/(M/2);

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'_step2_up.mat'];
    save(name,'cont2','pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit','pmit_eff','RNE_pmit')
end

if plot_on
    [VaR_outlier, ES_outlier] = Boxplot_PMitISEM(VaR_prelim,VaR_pmit,ES_prelim,ES_pmit,model,algo,H,N_sim,save_on);

    h_T = volatility_t_garch_noS_mex(draw_pmit(:,1:DD), data, S);
    [y_pmit, ~] = predict_t_garch_noS(draw_pmit(:,1:DD), y_T, h_T, H, draw_pmit(:,DD+1:H+DD));
    Plot_hor_pmit(y_pmit, y_T, mean(VaR_prelim),model,algo, save_on)

    Beta = Plot_beta(pmit,model,H,save_on,2); % the last parmeter: version==2 ==> plot only the second beta coefficient
end

%% Outliers detection
% % preliminary
% V1 = VaR_outlier{1,1};
% VO1 = zeros(length(VaR_outlier{1,1}),1);
% for ii = 1:length(V1)
%     VO1(ii,1) = find(VaR_prelim == V1(1,ii));
% end
% 
% E1 = ES_outlier{1,1};
% EO1 = zeros(length(ES_outlier{1,1}),1);
% for ii = 1:length(E1)
%     EO1(ii,1) = find(ES_prelim == E1(1,ii));
% end
% 
% % IS
% V2 = VaR_outlier{2,1};
% VO2 = zeros(length(VaR_outlier{2,1}),1);
% for ii = 1:length(V2)
%     VO2(ii,1) = find(VaR_IS == V2(1,ii));
% end
% 
% E2 = ES_outlier{2,1};
% EO2 = zeros(length(ES_outlier{2,1}),1);
% for ii = 1:length(E2)
%     EO2(ii,1) = find(ES_IS == E2(1,ii));
% end
% 
% ind_redo = [VO2; EO2]; % redo the 'VaR with PMit' loop for these indicies
%  
% %%% REDO and resave
% 
% if save_on
%     name = ['results/PMitISEM/',model,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
%     save(name,'VaR_prelim','ES_prelim','mit1','accept',...
%         'draw_hl','w_hl','lnk_hl','pmit','VaR_IS','ES_IS')
% end