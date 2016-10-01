clear all
close all

% horizon = [10,20,40, 100, 250];
% 
% 
% for H = horizon
% 
% close all

%% Initialization
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s); 

addpath(genpath('include/'));
plot_on = true;
save_on = true;

x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);

model = 'WN_ML';
algo = 'PMitISEM';
estimation = 'mle'; % 'true' or 'mle'

% Artificial, white noise data 
T = 10000;
sigma2_true = 1;
y = sqrt(sigma2_true)*randn(T,1); 
y = y - mean(y);

% Control parameters for PMitiISEM (cont2)
cont2 = MitISEM_Control;
cont2.mit.dfnc = 5;
cont2.mit.N = 10000;
cont2.mit.iter_max = 5;
cont2.df.range = [5,20];

 
sim = 1;
N_sim = 20;

M = 10000; % number of draws for IS computations

H = 250; % forecast horizon
p_bar = 0.01;
tail_shift = 1;

VaR_pmit = zeros(N_sim,1);
ES_pmit = zeros(N_sim,1);
time_pmit = zeros(2,1);

name =  ['results/PMitISEM/',model,'_Direct_',estimation,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
% if (H > 40)
    %% BIG DRAW
    load(name);
    draw_hl = draw_hl(:,2:H+1);
    N = size(draw_hl,1);
    w_hl = ones(N,1);
    % kernel = @(eps) -0.5*(H*log(2*pi) + H*log(sigma2_used) + sum(eps.^2,2)./sigma2_used);
    kernel = @(xx) -0.5*(H*log(2*pi) + H*log(1) + sum(xx.^2,2)./1);
    lnk_hl = kernel(draw_hl);

    draw0 = draw_hl;
    w0 = w_hl;
    lnk0 = lnk_hl; %kernel(draw0);
% else
%     load(name,'estimation','VaR_direct','ES_direct','sigma2_used')
%     mu_init = -ones(1,H);
%     kernel_init = @(xx) - MLtarget_WN_hl(xx, sigma2_used, tail_shift*mean(VaR_direct));
% end


%% PMIT ALL
% patial mixture - based on multidimensional stuctures with 4 mit fields
% mit_struct = struct('mu',[],'Sigma',[],'df',[],'p',[]);
% pmit_struct = struct('mu',cell(1,S),'Sigma',cell(1,S),'df',cell(1,S),'p',cell(1,S));

partition = 1:H;
fn_const_X = @(eps) WN_ML_const_X(eps,sigma2_used);
fn_input_X = @(xx) xx;

CV_old = cont2.mit.CV_old;
CV_tol = cont2.mit.CV_tol;
d = H;

cont2.mit.iter_max = 1;
% cont2.df.range = [1,20];
if (H < 250)
    cont2.mit.Hmax = 10;
else
    cont2.mit.Hmax = 2;    
end
cont2.mit.dfnc = 10;
cont = cont2;

kernel = @(xx) MLtarget_WN_hl(xx, sigma2_used, tail_shift*mean(VaR_direct));

tic
% if (H > 40)
    % [pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM_debug(draw0, lnk0, w0, kernel, fn_const_X, partition, d, cont, GamMat);
     [pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM2(draw0, lnk0, w0, kernel, fn_const_X, fn_input_X, partition, d, cont2, GamMat);
% else
%     [pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM2(kernel_init, mu_init, [], kernel, fn_const_X, fn_input_X, partition, d, cont2, GamMat);    
% end
time_pmit(1,1) = toc;

 if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',estimation,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','pmit','CV_mix','CV','iter')
 end
% load(name);

%% VaR with PMit
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s); 
pmit = pmit_step2;

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s); 
pmit = pmit_step2_up;

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s); 
pmit = pmit_step3;

tic
for sim = 1:N_sim 
    fprintf('\nVaR IS iter: %d\n',sim);
    
    draw_pmit = fn_p_rmvgt(M, pmit, H, partition, [], fn_const_X);  
    y_opt = sqrt(sigma2_used).*draw_pmit;  
    PL = fn_PL(y_opt);
%     kernel = @(xx) - 0.5*(log(2*pi) + log(sigma2_used) + (xx.^2)/sigma2_used);
    kernel = @(xx) - 0.5*(log(2*pi) + log(1) + (xx.^2)/1);
    lnk_opt = sum(kernel(draw_pmit),2);
    lnd_opt = fn_dpmit(draw_pmit, pmit, partition, fn_const_X, true, GamMat);
    w_opt = exp(lnk_opt - lnd_opt)/M;
    [PL, ind] = sort(PL);         
    w_opt = w_opt(ind,:);
    cum_w = cumsum(w_opt);
%     ind_var = min(find(cum_w >= p_bar))-1; 
    ind_var = min(find(cum_w > p_bar))-1; 
%     VaR_mit(sim,1) = PL(ind_var);
    VaR_pmit(sim,1)  = (PL(ind_var+1) + PL(ind_var))/2; % intrapolate
    ES = (w_opt(1:ind_var)/sum(w_opt(1:ind_var))).*PL(1:ind_var);
    ES_pmit(sim,1) = sum(ES(isfinite(ES)));       
    
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_pmit(sim,1), model, algo);  
end
time_pmit(2,1) = toc/N_sim;

VaR_step2 = VaR_pmit;
ES_step2 = ES_pmit;

VaR_step2_up = VaR_pmit;
ES_step2_up = ES_pmit;

VaR_step3 = VaR_pmit;
ES_ste3 = ES_pmit;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',estimation,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit')
end

if plot_on
    labels_in = {'naive','pmit'};
    algo_est = [algo,'_',estimation];
    Boxplot_PMitISEM(VaR_direct,VaR_pmit,ES_direct,ES_pmit,model,algo_est,H,N_sim,save_on, labels_in);

    y_pmit = sqrt(sigma2_used).*draw_pmit;  
    if (H>1)
        Plot_hor_pmit(y_pmit, y(end), tail_shift*mean(VaR_direct),model,algo,save_on)
    end
    
    Beta = Plot_beta(pmit,model,H,save_on,2); % the last parmeter: version==2 ==> plot only the second beta coefficient
end

% end