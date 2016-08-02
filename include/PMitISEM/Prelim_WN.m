clear all
close all

%% Initialization
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s); 

addpath(genpath('include/'));
plot_on = false;
save_on = false;

x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);

model = 'WN';
algo = 'Prelim';

% Artificial, white noise data 
T = 10000;
y = randn(T,1); 
y = y - mean(y);

% sigma is the VARIANCE of the error term, i.e. y_t ~ NID(0, sigma)
sigma_init = 0.9;

% Control parameters for MitISEM
cont1 = MitISEM_Control;
cont1.mit.dfnc = 5;
cont1.mit.Hmax = 2;
cont1.mit.N = 10000;
 
% hyper parameters for the prior for sigma2(inv. gamma)
a = 1; % if a == 0, then the flat prior is used; if a == 1, then the conjugate prior (inv. gamma)
b = 1; 
% posterior parameters:
a_post = a + T/2;
b_post = b + sum(y)/2;

sim = 1;
N_sim = 20;

% Metropolis-Hastings for the parameters
M = 10000; % number of draws for preliminary and IS computations
BurnIn = 1000;

H = 100; % forecast horizon
p_bar = 0.01;
% d = H+1; % dimension of theta

VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
RNE_prelim = zeros(N_sim,1);
accept = zeros(N_sim,1);
time_prelim = zeros(2,1);

% Construct the approximation to the parameter posterior
kernel_init = @(x) - posterior_WN(x, y, a, b, true);
kernel = @(x) posterior_WN(x, y, a, b, true);
tic
[mit1, summary1] = MitISEM_new(kernel_init, kernel, sigma_init, cont1, GamMat);
time_prelim(1,1) = toc;


%% Preliminary VaR
tic
for sim = 1:N_sim 
    fprintf('\nVaR prelim iter: %d\n',sim)
    kernel = @(x) posterior_WN(x, y, a, b, true);
    [sigma1, accept(sim,1)] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    fprintf('(%s) MH acceptance rate: %4.2f. \n', model, accept(sim,1));
    sigma1 = sigma1((BurnIn+1):(M+BurnIn));

    %% Future disturbances
    eps_H = randn(M,H); % --> if future disturbances drawn from the target then their weights are 1

    % the returns coresponding to disturbances: y = sqrt(sigma)*eps
    y_H = bsxfun(@times,eps_H,sqrt(sigma1)); 

    % preliminary VaR
    PL_H_ind = fn_PL(y_H);
    PL_H = sort(PL_H_ind);
    VaR_prelim(sim,1) = PL_H(p_bar*M);  
    ES_prelim(sim,1) = mean(PL_H(1:p_bar*M));   
    
    ind_prelim = double((PL_H_ind < VaR_prelim(sim,1)));
    RNE_prelim(sim,1) = fn_RNE(ind_prelim, 'MH',[],'Q');
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
end
time_prelim(2,1) = toc/N_sim;

if plot_on
    Plot_hor_direct(y_H,y(end), VaR_prelim(sim,1),model, save_on);
end

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','cont1','summary1','accept','time_prelim','RNE_prelim')
end

%% Generate many high loss draws to initialise the HL density approximation
% If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
kernel = @(xx) posterior_WN(xx, y, a, b, true);
y_predict = @(draw) bsxfun(@times,draw(:,2:end),sqrt(draw(:,1))); 
% cont1.mit.N = 10000; % the desired number of high-loss draws         
tic
[draw_hl, VaR_est, ~, ~] = BigDraw(cont1.mit.N, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat);
time_bigdraw = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_prelim','ES_prelim','mit1','cont1','summary1','accept','time_prelim','draw_hl','VaR_est','time_bigdraw','RNE_prelim')
end