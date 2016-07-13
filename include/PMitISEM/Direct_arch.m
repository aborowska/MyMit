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
algo = 'Direct';

% Direct
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
H = 10; % forecast horizon
% d = H+1; % dimension of theta
% partition = [1,3:H+1];

plot_on = true;
save_on = false;

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
accept_direct = zeros(N_sim,1);

kernel_init = @(a) - posterior_arch(a, data, S, true);
% kernel_init = @(a) - posterior_arch_noS(a, data, true);
cont_direct = MitISEM_Control;
cont_direct.mit.dfnc = 5;

tic
[mu, Sigma] = fn_initopt(kernel_init, mu_init);
mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',cont_direct.mit.dfnc);
time_direct(1,1)= toc;

tic
for sim = 1:N_sim
    fprintf('\nDirect sim = %i.\n', sim);
    kernel = @(a) posterior_arch(a, data, S, true);
    % kernel = @(a) posterior_arch_noS(a, data, true);
    [alpha_direct, accept_direct(sim,1)] = Mit_MH(M+BurnIn, kernel, mit_direct, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept_direct(sim,1), model, algo);
    alpha_direct = alpha_direct(BurnIn+1:M+BurnIn,:);

    eps_direct = randn(M,H);
    
    y_direct = predict_arch(alpha_direct, y_T, S, H, eps_direct);
%     y_direct = predict_arch_noS(alpha_direct, y_T, S, H, eps_direct);

    PL_direct = sort(fn_PL(y_direct));
    VaR_direct(sim,1) = PL_direct(p_bar*M);
    ES_direct(sim,1) = mean(PL_direct(1:p_bar*M));
    
    fprintf('Direct 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_direct(sim,1), model, algo);
end
time_direct(2,1) = toc/N_sim;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','mit_direct','accept_direct','time_direct')
end

