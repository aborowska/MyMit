clear all
close all

%% Initialization
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s); 

addpath(genpath('include/'));
plot_on = true;
save_on = false;

x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);

model = 'WN';
algo = 'Direct';

% Artificial, white noise data 
T = 10000;
y = randn(T,1); 
y = y - mean(y);

% sigma is the VARIANCE of the error term, i.e. y_t ~ NID(0, sigma)
sigma_init = 0.9;

% Control parameters for MitISEM (cont) and PMitiISEM (cont2)
cont_direct = MitISEM_Control;
cont_direct.mit.dfnc = 1;

% hyper parameters for the prior for sigma2(inv. gamma)
a = 1; % if a == 0, then the flat prior is used; if a == 1, then the conjugate prior (inv. gamma)
b = 1; 
% posterior parameters:
a_post = a + T/2;
b_post = b + sum(y)/2;

sim = 1;
N_sim = 20;

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
accept_direct = zeros(N_sim,1);
time_direct = zeros(2,1);

% Metropolis-Hastings for the parameters
M = 10000; % number of draws for preliminary and IS computations
BurnIn = 1000;

H = 10; % forecast horizon
p_bar = 0.01;


kernel_init = @(x) - posterior_WN(x, y, a, b, true);
tic
[mu, Sigma] = fn_initopt(kernel_init, sigma_init);
mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',cont_direct.mit.dfnc);
time_direct(1,1)= toc;

tic
for sim = 1:N_sim
    fprintf('\nDirect sim = %i.\n', sim);
    kernel = @(x) posterior_WN(x, y, a, b, true);
    [sigma_direct, accept_direct(sim,1)] = Mit_MH(M+BurnIn, kernel, mit_direct, GamMat);
    fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept_direct(sim,1), model, algo);
    sigma_direct = sigma_direct(BurnIn+1:M+BurnIn,:);

    eps_direct = randn(M,H);
    y_direct = bsxfun(@times,eps_direct,sqrt(sigma_direct)); 

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
