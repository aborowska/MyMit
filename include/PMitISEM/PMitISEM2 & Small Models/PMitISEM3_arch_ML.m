%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
 
x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 'arch_ML';
algo = 'PMitISEM';

data = csvread('GSPC_ret.csv');
data = 100*data;
ind_arch = find(data<=-5.5, 1, 'last' );
data = data(1:ind_arch,1);
data = data - mean(data);

T = length(data);
y_T = data(T);
S = var(data); % data variance for the variance targeting
        
p_bar = 0.01;
H = 100;

M = 10000;
N_sim = 20;
sim = 1;

plot_on = true;
save_on = false;

% Control parameters for PMitISEM
cont2 = MitISEM_Control;

VaR_pmit = zeros(N_sim,1);
ES_pmit = zeros(N_sim,1);
time_pmit = zeros(2,1);

%% PRELIM & BIG DRAW
name =  ['results/PMitISEM/',model,'_Direct_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name);

theta_mat = repmat(theta_mle,M,1);

draw_hl = draw_hl(:,2:end);
[N,d] = size(draw_hl);
w_hl = ones(N,1);
kernel = @(xx) -0.5*(H*log(2*pi) + H*log(1) + sum(xx.^2,2)./1);
lnk_hl = kernel(draw_hl);

%% PMitISEM
partition = 1:H;
d = H;

fn_const_X = @(xx) arch_ML_const_X3(xx,S);
fn_input_X = @(xx) arch_ML_input_X3(xx, theta_mle, y_T);

kernel_init = @(xx) - MLtarget_arch_hl(xx, theta_mle, y_T, S, mean(VaR_direct));
kernel = @(xx) MLtarget_arch_hl(xx, theta_mle, S, y_T, mean(VaR_direct));


CV_old = cont2.mit.CV_old;
CV_tol = cont2.mit.CV_tol;

draw0 = draw_hl;
w0 = w_hl;
lnk0 = lnk_hl; %kernel(draw0);

 if (H == 250)
    cont2.mit.Hmax = 1;
else
    cont2.mit.Hmax = 1;    
end

if (H >= 100)
    cont2.mit.dfnc = 20; 
    cont2.df.range = [15,25];
% elseif (H == 100)
%     cont2.mit.dfnc = 15; 
%     cont2.df.range = [5,20];
else
    cont2.mit.dfnc = 10; 
    cont2.df.range = [5,15];
end
cont = cont2;

tic
[pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM3(draw0, lnk0, w0, kernel, fn_const_X, fn_input_X, partition, d, cont2, GamMat);
time_pmit(1,1) = toc;

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','pmit','CV_mix','CV','iter')
end



%% VaR with PMit

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
pmit = pmit_step2;

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
pmit = pmit_step2_up;

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
pmit = pmit_step3;

tic
for sim = 1:N_sim
    fprintf('\nVaR IS iter: %d\n',sim)
    [draw_pmit, lnd_pmit, input_X_pmit] = fn_p_rmvgt_dpmit3(M, pmit,  d, SS, partition, fn_const_X, fn_input_X, GamMat);
    y_pmit = input_X_pmit.y_cum;
    kernel = @(xx) - 0.5*(log(2*pi) + log(1) + (xx.^2)/1);
    lnk_pmit = sum(kernel(draw_pmit),2);
    PL_pmit = fn_PL(y_pmit);
    
    w_pmit = exp(lnk_pmit - lnd_pmit)/M;
    [PL, ind] = sort(PL_pmit);         
    w_pmit = w_pmit(ind,:);
    cum_w = cumsum(w_pmit);
    ind_var = min(find(cum_w >= p_bar))-1; 
    VaR_pmit(sim,1) = PL(ind_var);
    ES = (w_pmit(1:ind_var)/sum(w_pmit(1:ind_var))).*PL(1:ind_var);
    ES_pmit(sim,1) = sum(ES(isfinite(ES)));    
    
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_pmit(sim,1), model, algo); 
end
time_pmit(2,1) = toc/N_sim;

VaR_step2 = VaR_pmit;
ES_step2 = ES_pmit;

% VaR_pmit = VaR_step2;s
% ES_pmit = ES_step2;

VaR_step2_up = VaR_pmit;
ES_step2_up = ES_pmit;

VaR_step3 = VaR_pmit;
ES_step3 = ES_pmit;

% VaR_pmit = VaR_step2_up;
% ES_pmit = ES_step2_up;

% time_pmit(1,1) = time_pmit(1,1) + time_step2_up;


if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit')
end


y_pmit = predict_arch(theta_mat, y_T, S, H, draw_pmit);
PL_pmit = fn_PL(y_pmit);
pmit_eff = sum(PL_pmit <= mean(VaR_direct))/M;


if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','pmit','CV_mix','CV','iter','VaR_pmit','ES_pmit','time_pmit','pmit_eff')
end

if plot_on
    labels_in = {'naive','pmit'};
    [VaR_outlier, ES_outlier] = Boxplot_PMitISEM(VaR_direct, VaR_pmit, ES_direct, ES_pmit, model, algo, H, N_sim, true, labels_in);

    [y_pmit, ~] = predict_t_gas(theta_mat, y_T, f_mat, H, draw_pmit);
    Plot_hor_pmit(y_pmit, y_T, mean(VaR_direct), model, algo, save_on)

    [Beta, Sigma, nu] = Plot_beta(pmit,model,H,save_on,2); % the last parmeter: version==2 ==> plot only the second beta coefficient
end