%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 't_garch2_noS_ML';
algo = 'Direct';

crisis = false;
recent = false;
old = true;
if crisis 
    y = csvread('GSPC_ret_updated.csv'); 
    results_path = 'results/PMitISEM/crisis/';
elseif recent
    y = csvread('GSPC_ret_updated_short.csv');
    results_path = 'results/PMitISEM/recent';
elseif old
    y = csvread('GSPC_ret_tgarch.csv');
    results_path = 'results/PMitISEM/old/';        
else
    y = csvread('GSPC_ret_updated_short_end.csv');
    results_path = results_path;    
end
y = 100*y;
data = y;

T = size(data,1);
y_T = data(T);
S = var(data);
 
M = 10000;
N_sim = 20;

plot_on = true;
save_on = true;

p_bar = 0.01;
H = 1;     % prediction horizon 


VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
time_direct = zeros(2,1);

hyper = 0.01; 
% kernel_init = @(a) - posterior_t_garch_noS_hyper_mex(a, data , S, GamMat, hyper);
kernel_init = @(a) posterior_t_garch_noS_hyper_init_mex(a, data , S, GamMat, hyper);
kernel = @(a) posterior_t_garch_noS_hyper_mex(a, data, S, GamMat, hyper);


% theta = [omega, alpha, beta, mu, nu]
if crisis
    mu_init = [0.008, 0.07, 0.9, 0.01, 6.2];
    %     mu_init = [0.02, 0.12, 0.85, 0.075, 6.2];
    tic
%     x = fminsearch(kernel_init,mu_init);
%     [mu, Sigma] = fn_initopt(kernel_init, x);
%     [mu, Sigma] = fn_initopt(kernel_init, mu);
%     [mu, Sigma] = fn_initopt(kernel_init, mu);
    theta_mle = fn_initopt(kernel_init, mu_init);
    time_direct(1,1) = toc;   
elseif recent
    % mu_init = [0.008, 0.07, 0.9, 0.01, 10];
%     mu_init = [0.009, 0.07, 0.9, 0.05, 11];    
    mu_init = [0.043, 0.17, 0.78, 0.08, 6.1];    
    tic
    theta_mle  = fn_initopt(kernel_init, mu_init);
    time_direct(1,1) = toc; 
elseif old
% % %     mu_init = [0.009, 0.07, 0.9, 0.05, 11];
% %     mu_init = [0.009, 0.06, 0.9, 0.05, 11];
%     mu_init = [0.006, 0.065, 0.92, 0.048, 10.0];    
    mu_init = [0.008, 0.07, 0.92, 0.048, 10.0];    
    tic
    theta_mle = fn_initopt(kernel_init, mu_init);
    Sigma = Sigma/T;
    time_direct(1,1) = toc;
else
%     mu_init = [0.043, 0.17, 0.78, 0.08, 6.1];    
    mu_init = [0.05, 0.15, 0.78, 0.09, 6.1];    
    tic
    theta_mle = fn_initopt(kernel_init, mu_init);
    time_direct(1,1) = toc;     
end
mu_init = [0.009, 0.07, 0.9, 0.05, 11];
d = size(mu_init,2);


h_mle = volatility_t_garch_noS_mex(theta_mle, y, S);
theta_direct = repmat(theta_mle, M, 1);
h_direct = repmat(h_mle, M, 1);

tic
for sim = 1:N_sim
    fprintf('\nDirect sim = %i.\n', sim);
    [y_direct, eps_direct] = predict_t_garch_noS(theta_direct, y_T, h_direct, H);

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
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','theta_mle','h_mle','time_direct')
end


%% Generate many high loss draws to initialise the HL density approximation
% If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
y_predict = @(draw) predict_t_garch_new_noS(draw(:,1:d), y, S, H, draw(:,d+1:end));

tic
[draw_hl, VaR_est, ~, ~] = BigDraw(10000, H, [], p_bar, theta_mle, [], y_predict, GamMat, 5);
time_bigdraw = toc;

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','theta_mle','h_mle','time_direct','draw_hl','VaR_est','time_bigdraw')
end