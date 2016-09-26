%% Initialisation
% clc
clear all
close all
s = RandStream('mt19937ar','Seed',1); % or: 0
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));

x = (0:0.00001:50)' + 0.00001;
GamMat = gamma(x);

model = 'arch_ML';
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
        
M = 10000; % number of draws for preliminary and IS computations

N_sim = 20;
p_bar = 0.01;
H = 10; % forecast horizon
% d = H+1; % dimension of theta
% partition = [1,3:H+1];

plot_on = true;
save_on = true;

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
time_direct = zeros(2,1);

kernel_init = @(a) - posterior_arch(a, data, S, true);
% kernel_init = @(a) - posterior_arch_noS(a, data, true);
mu_init = 0.03;


tic
theta_mle = fn_initopt(kernel_init, mu_init);
time_direct(1,1)= toc;

theta_direct = repmat(theta_mle, M, 1);


tic
for sim = 1:N_sim
    fprintf('\nDirect sim = %i.\n', sim);
    
    [y_direct, eps_direct] = predict_arch(theta_direct, y_T, S, H);
%     y_direct = predict_arch_noS(alpha_direct, y_T, S, H, eps_direct);
    
    ind_real = find(sum(imag(y_direct),2)==0);
    M_real = length(ind_real); 
    fprintf('M_real = %i.\n',M_real)
    y_direct = y_direct(ind_real,:);
    theta_direct = theta_direct(ind_real,:);  
    eps_direct = eps_direct(ind_real,:);

    [PL_direct, ind] = sort(fn_PL(y_direct));
    VaR_direct(sim,1) = PL_direct(round(p_bar*M_real));
    ES_direct(sim,1) = mean(PL_direct(round(1:p_bar*M_real)));   
     
    fprintf('Direct 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_direct(sim,1), model, algo);
end
time_direct(2,1) = toc/N_sim;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','theta_mle','time_direct')
end


%% Generate many high loss draws to initialise the HL density approximation
% If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
y_predict = @(draw) predict_arch(draw(:,1), y_T, S, H, draw(:,2:end));  
tic
[draw_hl, VaR_est, ~, ~] = BigDraw(10000, H, [], p_bar, theta_mle, [], y_predict, GamMat);
time_bigdraw = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','theta_mle','time_direct','draw_hl','VaR_est','time_bigdraw')
end