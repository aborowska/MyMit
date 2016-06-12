%% Initialisation
% clc
clear all
close all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));

x = (0:0.00001:50)' + 0.00001;
GamMat = gamma(x);

model = 'arch';
algo = 'Prelim';

% Data
data = csvread('GSPC_ret.csv');
data = 100*data;
ind_arch = find(data<=-5.5, 1, 'last' );
data = data(1:ind_arch,1);
data = data - mean(data);

T = length(data);
y_T = data(T);
S = var(data); % data variance for the variance targeting
         
mu_init = 0.03;
 
% Metropolis-Hastings for the preliminary
M = 10000; % number of draws for preliminary and IS computations
BurnIn = 1000;

N_sim = 20;
p_bar = 0.01;

H = 20; % forecast horizon
plot_on = true;
save_on = true;

% Control parameters for MitISEM (cont) and PMitiISEM (cont2)
cont1 = MitISEM_Control;
cont1.mit.dfnc = 5;
cont1.mit.N = 10000;


VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
accept = zeros(N_sim,1);
time_prelim = zeros(2,1);

kernel_init = @(a) - posterior_arch(a, data, S, true);
kernel = @(a) posterior_arch(a, data, S, true);
tic
[mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont1, GamMat);
time_prelim(1,1) = toc;

%% QERMit 1a.: 
tic
for sim = 1:N_sim
    fprintf('\nPrelim sim = %i.\n', sim);
    kernel = @(a) posterior_arch(a, data, S, true);

    [alpha1, accept(sim,1)] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept(sim,1), model, algo);
    alpha1 = alpha1(BurnIn+1:M+1000);

    eps1 = randn(M,H);
    y_H = predict_arch(alpha1, y_T, S, H, eps1);
    % get the preliminary VaR estimate as the 100th of the ascendingly sorted percentage loss values
    [PL_H, ind] = sort(fn_PL(y_H));
    VaR_prelim(sim,1) = PL_H(p_bar*M);
    ES_prelim(sim,1) = mean(PL_H(1:p_bar*M));    
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
end
time_prelim(2,1) = toc/N_sim;

if plot_on
    Plot_hor_direct(y_H,y_T,VaR_prelim(sim,1),model,save_on);
end

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','cont1','summary1','accept','time_prelim')
end

%% Choose the starting point (mu_hl) for the constuction of the approximaton
 % If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
kernel = @(xx) posterior_arch(xx, data, S, true);
y_predict = @(draw) predict_arch(draw(:,1), y_T, S, H, draw(:,2:end));  
% cont1.mit.N =10000; % the desired number of high-loss draws         
tic
[draw_hl, VaR_est, ~, ~] = BigDraw(cont1.mit.N, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat);
time_bigdraw = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','cont1','summary1','accept','time_prelim','draw_hl','VaR_est','time_bigdraw')
end