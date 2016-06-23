%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
 
x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

y = csvread('GSPC_ret_tgarch.csv');
y = 100*y;
T = size(y,1);
y_T = y(T,1);


model = 't_gas';
algo = 'Prelim';


plot_on = true;
save_on = true;

% Control parameters for MitISEM
cont1 = MitISEM_Control;
cont1.mit.dfnc = 5;
cont1.mit.N = 10000;
cont1.resmpl_on = false;
 
cont2 = MitISEM_Control;
cont2.mit.dfnc = 5;

    
M = 10000;
BurnIn = 1000;
N_sim = 20;
p_bar = 0.01;
H = 10;


VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
accept = zeros(N_sim,1);
time_prelim = zeros(2,1);


% theta = [mu, omega, A, B, nu]
mu_init = [0, 0.01, 0.1, 0.89, 8];
DD = size(mu_init,2);



hyper = 0.01;
% kernel_init = @(xx) - posterior_t_gas(xx, y, hyper, true, GamMat);
% tic
% kernel = @(xx) posterior_t_gas(xx, y, hyper, true, GamMat);
% kernel(mu_init)
% toc

kernel_init = @(xx) - posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
% tic
kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
% kernel(mu_init)
% toc
% [mu, Sigma] = fn_initopt(kernel_init, mu_init);

[mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont1, GamMat);

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'mit1','summary1')
end

% load(name)

tic
for sim = 1:N_sim  
    fprintf('\nPrelim sim = %i.\n', sim);
    kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
    [theta1, accept(sim,1)] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept(sim,1), model, algo);
    theta1 = theta1(BurnIn+1:M+BurnIn,:);

    f_T = volatility_t_gas_mex(theta1, y);
    [y_H, eps_H] = predict_t_gas(theta1, y_T, f_T, H);

    ind_real = (sum(imag(y_H),2)==0);
    M_real = sum(ind_real); 
    fprintf('M_real = %i.\n',M_real)
    y_H = y_H(ind_real,:);
    theta1 = theta1(ind_real,:);  
    eps_H = eps_H(ind_real,:);

    [PL_H, ind] = sort(fn_PL(y_H));
    VaR_prelim(sim,1) = PL_H(round(p_bar*M_real));
    ES_prelim(sim,1) = mean(PL_H(round(1:p_bar*M)));   
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
end    


if plot_on
    Plot_hor_direct(y_H,y_T,VaR_prelim(sim,1),model,save_on);
end

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','cont1','summary1','accept','time_prelim')
end
% load(name)
% 
% % If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
% % kernel = @(a) posterior_t_garch_noS(a, data, S, L, hyper, GamMat);
% % kernel = @(xx) posterior_t_garch_noS_mex(xx, data, S, GamMat);
kernel = @(xx) posterior_t_gas_hyper_mex(xx, y, hyper, GamMat);
y_predict = @(draw) predict_t_gas_new(draw(:,1:d), y, H, draw(:,d+1:end));
tic
[draw_hl, VaR_est, ~, ~] = BigDraw(cont1.mit.N, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat, d);
time_bigdraw = toc;
 
if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','cont1','summary1','accept','time_prelim','draw_hl','VaR_est','time_bigdraw')
end

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

if (H > 10)
    cont2.mit.iter_max = 1;
else
    cont2.mit.iter_max = 5;
end
cont2.df.range = [1,10];
cont2.mit.Hmax = 1;
cont = cont2;

tic
% [pmit, CV_mix, CV, iter, pmit_pre, pmit_pre2, pmit_adapt]  = PMitISEM(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont2, GamMat);
[pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM2(draw0, lnk0, w0, kernel, fn_const_X, fn_input_X, partition, d, cont2, GamMat);
time_pmit(1,1) = time_pmit(1,1) + toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'pmit','CV_mix','CV','iter')
end