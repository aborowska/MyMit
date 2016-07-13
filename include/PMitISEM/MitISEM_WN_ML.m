clear all
close all

%% Initialization
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s); 

addpath(genpath('include/'));
plot_on = true;
save_on = true;

x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);

model = 'WN_ML';
algo = 'MitISEM';

% Artificial, white noise data 
T = 10000;
y = randn(T,1); 
y = y - mean(y);

% sigma is the VARIANCE of the error term, i.e. y_t ~ NID(0, sigma)
sigma_init = 0.9;

% Control parameters for MitISEM (cont) and PMitiISEM (cont2)
cont2 = MitISEM_Control;
cont2.mit.dfnc = 5;
cont2.mit.N = 10000;
cont2.mit.iter_max = 5;
cont2.df.range = [1,10];

sim = 1;
N_sim = 20;

VaR_mit = zeros(N_sim,1);
ES_mit = zeros(N_sim,1);
time_mit = zeros(2,1);

M = 10000; % number of draws for preliminary and IS computations

H = 10; % forecast horizon
p_bar = 0.01;

%% BIG DRAW
name =  ['results/PMitISEM/',model,'_Direct_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name);
draw_hl = draw_hl(:,2:end);
[N,d] = size(draw_hl);
w_hl = ones(N,1);

%% Standard MitISEM
[mu_hl, Sigma_hl] = fn_muSigma(draw_hl, w_hl);

% cont2.mit.N = 10000;
cont2.mit.Hmax = 10;
cont = cont2;

mit_hl.mu = mu_hl;
mit_hl.Sigma = Sigma_hl;
mit_hl.df = cont2.mit.dfnc;
mit_hl.p = 1;
% mu_init = mu_hl;
% mit_init = mit_hl;

% kernel_init = @(x) - posterior_debug_hl(x, y, 0, 0, mean(VaR_prelim), true); 
% kernel = @(x) posterior_debug_hl(x, y, 0, 0, mean(VaR_prelim), true); 

kernel_init = @(xx) -  MLtarget_WN_hl(xx, sigma2_mle, mean(VaR_direct));
kernel = @(xx)  MLtarget_WN_hl(xx, sigma2_mle, mean(VaR_direct));
 

tic
if (H <= 100)
    [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);
else
    [mit2, summary2] = MitISEM_new(mit_hl, kernel, mu_hl, cont2, GamMat);
end
time_mit(1,1) = time_mit(1,1) + toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'mit2','summary2')
end

% kernel_init = @(x) - posterior_debug(x, y, a, b, true);
% [mu, Sigma] = fn_initopt(kernel_init, sigma_init);
% mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',5);

%% VaR with standard MitISEM
tic
for sim = 1:N_sim 
    fprintf('\nVaR IS iter: %d\n',sim)

    draw_opt = rmvgt2(M, mit2.mu, mit2.Sigma, mit2.df, mit2.p); 

    %% IS weights
    kernel = @(xx) - 0.5*(log(2*pi) + xx.^2);
    lnk_opt = sum(kernel(draw_opt),2);
    lnd_opt = dmvgt(draw_opt, mit2, true, GamMat);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % IS VaR estimation
    y_opt = sigma2_mle.*draw_opt;  
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_mit(sim,1) = IS_estim(1,1);
    ES_mit(sim,1) = IS_estim(1,2);   
  
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_mit(sim,1), model, algo);  
end
time_mit(2,1) = toc/N_sim;


if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_mit','ES_mit','mit2','summary2','time_mit')
end


labels_in = {'prelim','mitisem'};
Boxplot_PMitISEM(VaR_prelim, VaR_mit, ES_prelim, ES_mit, model, algo, H, N_sim, true, labels_in);


y2 = bsxfun(@times,draw2(:,2:d),sqrt(draw2(:,1)));  
Plot_hor_pmit(y2, y(end), mean(VaR_prelim),model,algo,save_on)