% horizon=[10,20,40,100,250];
% 
% for H = horizon
clear all
close all

%% Initialization
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s); 

addpath(genpath('include/'));
plot_on = false;
save_on = true;

x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);

model = 'WN_ML';
algo = 'Direct';

estimation = 'true'; % 'true' or 'mle'

% Artificial, white noise data 
T = 10000;
sigma2_true = 1;
y = sqrt(sigma2_true)*randn(T,1); 
y = y - mean(y);

% sigma2 is the VARIANCE of the error term, i.e. y_t ~ NID(0, sigma2)
sigma2_init = 0.9;

sim = 1;
N_sim = 20;

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
time_direct = zeros(2,1);

M = 10000; % number of draws for preliminary and IS computations

H = 40; % forecast horizon
p_bar = 0.01;

% loglik = @(ss,yy) - 0.5*sum(log(2*pi) + log(ss) + (y.^2)./ss);
kernel_init = @(ss) 0.5*(T*log(2*pi) + T*log(ss) + sum(y.^2)./ss);

if strcmp(estimation, 'mle');
    tic
    sigma2_mle = fn_initopt(kernel_init, sigma2_init);
    time_direct(1,1) = toc;  
    sigma2_used = sigma2_mle;
else
    time_direct(1,1) = 0;  
    sigma2_used = sigma2_true; 
end  

if plot_on
    figure(99) 
    set(gcf,'defaulttextinterpreter','latex');  
    xx = 0.5:0.01:1.5;
    yy = kernel_init(xx);
    plot(xx,-yy);
    xlabel('$$\sigma^2$$','FontSize', 12) % x-axis label
    ylabel('Loglikelihood','FontSize', 12) % y-axis label
    plotTickLatex2D('FontSize',12);
end

tic
for sim = 1:N_sim
    fprintf('\nDirect sim = %i.\n', sim); 
    y_direct = sqrt(sigma2_used)*randn(M,H);
    PL_direct = sort(fn_PL(y_direct));
    VaR_direct(sim,1) = PL_direct(p_bar*M);
    ES_direct(sim,1) = mean(PL_direct(1:p_bar*M));
    
    fprintf('Direct 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_direct(sim,1), model, algo);
end
time_direct(2,1) = toc/N_sim;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',estimation,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','sigma2_used','estimation','time_direct')
end

%% Generate many high loss draws to initialise the HL density approximation
% If we want many draws (to obtain a better approximation) better use BigDraw function (memory considerations)
y_predict = @(draw) bsxfun(@times,draw(:,2:end),sqrt(draw(:,1))); 
tic
[draw_hl, VaR_est, ~, ~] = BigDraw(10000, H, [], p_bar, sigma2_used, [], y_predict, GamMat);
time_bigdraw = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',estimation,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_direct','ES_direct','sigma2_used','time_direct','draw_hl','VaR_est','time_bigdraw','estimation')
end

% end