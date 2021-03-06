%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
 
x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 't_gas';
algo = 'PMitISEM';

y = csvread('GSPC_ret_tgarch.csv');
y = 100*y;
T = size(y,1);
y_T = y(T,1);

p_bar = 0.01;
H = 100;

M = 10000;
BurnIn = 1000;
N_sim = 20;
sim = 1;

% theta = [mu, omega, A, B, nu]
mu_init = [0, 0.01, 0.1, 0.89, 8];
DD = size(mu_init,2);

plot_on = true;
save_on = false;

% Control parameters for PMitISEM
cont2 = MitISEM_Control;
cont2.mit.dfnc = 5;

VaR_pmit = zeros(N_sim,1);
ES_pmit = zeros(N_sim,1);
RNE_pmit = zeros(N_sim,1);
time_pmit = zeros(2,1);

%% PRELIM & BIG DRAW
name =  ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name);

hyper = 0.01;
kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);

lnk_hl = kernel(draw_hl(:,1:DD)); 
% log candidate evaluation
lnd_hl = dmvgt(draw_hl(:,1:DD), mit1, true, GamMat);

% importance weights
w_hl = lnk_hl - lnd_hl;
w_hl = exp(w_hl - max(w_hl));

%% PMitISEM
partition = [1,DD+2:H+DD];
d = H+DD;

fn_const_X = @(xx) t_gas_hyper_const_X2(xx, y);
fn_input_X = @(xx) t_gas_hyper_input_X(xx, y);
% kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
kernel = @(a) posterior_t_gas_hl_hyper_mex(a, y, hyper, mean(VaR_prelim), GamMat);


CV_old = cont2.mit.CV_old;
CV_tol = cont2.mit.CV_tol;

draw0 = draw_hl;
w0 = w_hl;
lnk0 = lnk_hl; %kernel(draw0);

% if (H > 10)
%     cont2.mit.iter_max = 1;
% else
    cont2.mit.iter_max = 1;%6;%8;
% end
cont2.df.range = [5,15];
cont2.mit.Hmax = 1;
% if (H >= 40)
    cont2.mit.dfnc = 10;
% end
if (H==250)
    cont2.mit.Hmax = 1;
end
cont = cont2;

tic
% [pmit, CV_mix, CV, iter, pmit_pre, pmit_pre2, pmit_adapt]  = PMitISEM(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont2, GamMat);
[pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM2(draw0, lnk0, w0, kernel, fn_const_X, fn_input_X, partition, d, cont2, GamMat);
time_pmit(1,1) = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'pmit','CV_mix','CV','iter')
end



%% VaR with PMit

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
pmit = pmit_step2;

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
pmit = pmit_step2_up;

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

%     draw_pmit  = fn_p_rmvgt(M/2, pmit, d, partition, [], fn_const_X);  
    [draw_pmit, ~, input_X_pmit] = fn_p_rmvgt2(M/2, pmit, d, partition, [], fn_const_X, fn_input_X);         

    draw_opt = [draw1; draw_pmit];   
    
    input_X.theta = draw_opt;
    input_X.f_T = [input_X_1.f_T; input_X_pmit.f_T];
    clear draw1 %draw_pmit
    clear input_X_1 input_X_pmit
    
    kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
    lnk_opt = kernel(draw_opt(:,1:5)); 

    eps_pdf = duvt(draw_opt(:,DD+1:H+DD), draw_opt(:,DD), H, true);
    lnk_opt = lnk_opt + eps_pdf;

    % optimal weights
    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,1:DD), mit1, true, GamMat));
%     exp_lnd2 = fn_dpmit(draw_opt, pmit, partition, fn_const_X, true, GamMat);
    exp_lnd2 = fn_dpmit2(input_X, pmit, partition, fn_const_X, true, GamMat);        
    
    exp_lnd2 = 0.5*exp(exp_lnd2);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % IS VaR estimation
%     f_T = volatility_t_gas_mex(draw_opt(:,1:DD), y);
    f_T = input_X.f_T;
    y_opt = predict_t_gas(draw_opt(:,1:DD), y_T, f_T, H, draw_opt(:,DD+1:H+DD));
    ind_opt = (fn_PL(y_opt) <= mean(VaR_prelim));
    RNE_pmit(sim,1) = fn_RNE(ind_opt, 'IS', w_opt); 
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_pmit(sim,1) = IS_estim(1,1);
    ES_pmit(sim,1) = IS_estim(1,2);   

    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_pmit(sim,1), model, algo);  
end
time_pmit(2,1) = toc/N_sim;

VaR_step2 = VaR_pmit;
ES_step2 = ES_pmit;

VaR_step2_up = VaR_pmit;
ES_step2_up = ES_pmit;

% time_pmit(1,1) = time_pmit(1,1) + time_step2_up;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit','RNE_pmit')
end


f_pmit = volatility_t_gas_mex(draw_pmit(:,1:DD), y);
y_pmit = predict_t_gas(draw_pmit(:,1:DD), y_T, f_pmit, H, draw_pmit(:,DD+1:H+DD));
PL_pmit = fn_PL(y_pmit);
pmit_eff = sum(PL_pmit <= mean(VaR_prelim))/(M/2);

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit','pmit_eff','RNE_pmit')
end

if plot_on
    [VaR_outlier, ES_outlier] = Boxplot_PMitISEM(VaR_prelim,VaR_pmit,ES_prelim,ES_pmit,model,algo,H,N_sim,true);

    f_T = volatility_t_gas_mex(draw_pmit(:,1:DD), y);
    [y_pmit, ~] = predict_t_gas(draw_pmit(:,1:DD), y_T, f_T, H, draw_pmit(:,DD+1:H+DD));
    Plot_hor_pmit(y_pmit, y_T, mean(VaR_prelim),model,algo, save_on)

    [Beta, Sigma, nu] = Plot_beta(pmit,model,H,save_on,2); % the last parmeter: version==2 ==> plot only the second beta coefficient
end