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
algo = 'PMitISEM';

% Artificial, white noise data 
T = 10000;
y = randn(T,1); 
y = y - mean(y);

% sigma is the VARIANCE of the error term, i.e. y_t ~ NID(0, sigma)
sigma_init = 0.9;

% Control parameters for MitISEM (cont) and PMitiISEM (cont2)
cont2 = MitISEM_Control;
cont2.mit.N = 10000;

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

H = 20; % forecast horizon
p_bar = 0.01;
% d = H+1; % dimension of theta

VaR_pmit = zeros(N_sim,1);
ES_pmit = zeros(N_sim,1);
RNE_pmit = zeros(N_sim,1);
time_pmit = zeros(2,1);

% PRELIM & BIG DRAW
name =  ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name);

% WEIGHTS to initialise PMitISEM
% future disturbances are generated from the target thus have weights 1
% log kernel evaluation - only for the parameter draws
kernel = @(xx) posterior_WN(xx, y, a, b, true);
lnk_hl = kernel(draw_hl(:,1)); 
% log candidate evaluation
lnd_hl = dmvgt(draw_hl(:,1), mit1, true, GamMat);
% weights
w_hl = lnk_hl - lnd_hl;
w_hl = exp(w_hl - max(w_hl));

%% PMIT ALL
% patial mixture - based on multidimensional stuctures with 4 mit fields
% mit_struct = struct('mu',[],'Sigma',[],'df',[],'p',[]);
% pmit_struct = struct('mu',cell(1,S),'Sigma',cell(1,S),'df',cell(1,S),'p',cell(1,S));

partition = [1,3:H+1];
SS = length(partition);
d = H+1;


fn_input_X = @(xx) WN_input_X3(xx);
fn_const_X = @(xx) WN_const_X3(xx);
 

kernel = @(x) posterior_WN_hl(x, y, a, b, mean(VaR_prelim), true); 
CV_old = cont2.mit.CV_old;
CV_tol = cont2.mit.CV_tol;

draw0 = draw_hl;
w0 = w_hl;
lnk0 = lnk_hl; %kernel(draw0);
% clear draw_hl w_hl lnk_hl lnd_hl

cont2.mit.iter_max = 1;
if (H == 100)
    cont2.df.range = [10,20];
    cont2.mit.dfnc = 13; % <---------
else
    cont2.df.range = [10,20];
    cont2.mit.dfnc = 15; % <---------
end
cont2.mit.Hmax = 1;
cont = cont2;

tic
[pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM3(draw0, lnk0, w0, kernel, fn_const_X, fn_input_X, partition, d, cont2, GamMat);
time_pmit(1,1) = toc;

 if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'pmit','CV_mix','CV','iter')
 end
% load(name);

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
     
    sigma1 = rmvgt2(M/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
    eps1 = randn(M/2,H); 
    draw1 = [sigma1, eps1];
    input_X_1 = fn_input_X(draw1);
    [lnd1, input_X_1] = fn_dpmit3(input_X_1, pmit, partition, fn_const_X, true, GamMat);        

%     [draw_pmit, ~, input_X_pmit] = fn_p_rmvgt3(M/2, pmit, d, partition, [], fn_const_X, fn_input_X);         
    [draw_pmit, lnd_pmit, input_X_pmit] = fn_p_rmvgt_dpmit3(M/2, pmit,  d, SS, partition, fn_const_X, fn_input_X, GamMat);
    
    draw_opt = [draw1; draw_pmit];
    
    kernel = @(x) posterior_WN(x, y, a, b, true);
    lnk_opt = kernel(draw_opt(:,1));

    kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
    eps_pdf = sum(kernel(draw_opt(:,2:d)),2);
    lnk_opt = lnk_opt + eps_pdf;

    % optimal weights 
    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,1), mit1, true, GamMat));
%     exp_lnd2 = fn_dpmit3(input_X, pmit, partition, fn_const_X, true, GamMat);        
    exp_lnd2 = [lnd1; lnd_pmit];        
    exp_lnd2 = 0.5*exp(exp_lnd2);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % IS VaR estimation
%     y_opt = bsxfun(@times,draw_opt(:,2:d),sqrt(draw_opt(:,1))); 
    y_opt = [input_X_1.y_cum; input_X_pmit.y_cum];
%     ind_opt = (fn_PL(y_opt) <= mean(VaR_prelim));
%     RNE_pmit(sim,1) = fn_RNE(ind_opt, 'IS', w_opt);     
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_pmit(sim,1) = IS_estim(1,1);
    ES_pmit(sim,1) = IS_estim(1,2);   
  
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_pmit(sim,1), model, algo);  
end
time_pmit(2,1) = toc/N_sim;

VaR_step2 = VaR_pmit;
ES_step2 = ES_pmit;

% VaR_pmit = VaR_step2;
% ES_pmit = ES_step2;

VaR_step2_up = VaR_pmit;
ES_step2_up = ES_pmit;

% time_pmit(1,1) = time_pmit(1,1) + time_step2_up;

y_pmit = bsxfun(@times,draw_pmit(:,2:d),sqrt(draw_pmit(:,1)));  
PL_pmit = fn_PL(y_pmit);
pmit_eff = sum(PL_pmit <= mean(VaR_prelim))/(M/2);

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit','pmit_eff','RNE_pmit')
end

if plot_on
    Boxplot_PMitISEM(VaR_prelim,VaR_pmit,ES_prelim,ES_pmit,model,algo,H,N_sim,true);

    y_pmit = bsxfun(@times,draw_pmit(:,2:d),sqrt(draw_pmit(:,1)));  
    Plot_hor_pmit(y_pmit, y(end), mean(VaR_prelim),model,algo,save_on)

    Beta = Plot_beta(pmit,model,H,save_on,2); % the last parmeter: version==2 ==> plot only the second beta coefficient
end