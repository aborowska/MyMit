%% Initialisation
clear all
close all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

model = 't_garch2_noS_ML';
algo = 'MitISEM';

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

M = 10000;
BurnIn = 1000;
N_sim = 20;


% theta = [mu, omega, A, B, nu]
mu_init = [0, 0.01, 0.1, 0.89, 8];
d = size(mu_init,2);

plot_on = false;
save_on = true;

cont2 = MitISEM_Control;
% cont2.mit.dfnc = 5;
cont2.mit.dfnc = 10; % <--!!

cont2.mit.Hmax = 10;
cont2.df.range = [5, 15];

p_bar = 0.01;
H = 10;     % prediction horizon 

VaR_mit = zeros(N_sim,1);
ES_mit = zeros(N_sim,1);
time_mit = zeros(2,1);

%% PRELIM & BIG DRAW
name =  [results_path,model,'_Direct_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name);

theta_mat = repmat(theta_mle,M,1);
h_mat = repmat(h_mle,M,1);    
nu_mat = theta_mat(:,5);


draw_hl = draw_hl(:,6:end);
[N,d] = size(draw_hl);
w_hl = ones(N,1);


[mu_hl, Sigma_hl] = fn_muSigma(draw_hl, w_hl);

mit_hl.mu = mu_hl;  
mit_hl.Sigma = Sigma_hl;
mit_hl.df = cont2.mit.dfnc;
mit_hl.p = 1;
% mu_init = mu_hl;

if (H < 40)
    cont2.mit.Hmax = 10;
    cont2.mit.iter_max = 10;
else
    cont2.mit.Hmax = 1;  % <<<<<< !!
    cont2.mit.iter_max = 1;
end
% cont = cont2;

kernel_init = @(xx) - MLtarget_t_garch_noS_hl(xx, theta_mle, h_mle, y_T, mean(VaR_direct));
kernel = @(xx) MLtarget_t_garch_noS_hl(xx, theta_mle, h_mle, y_T, mean(VaR_direct));

tic
if (H < 20)
    [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);
else
    [mit2, summary2] = MitISEM_new(mit_hl, kernel, mu_hl, cont2, GamMat);
end
time_mit(1,1) = toc;


if (plot_on && (H == 1))
    figure(100) 
    set(gcf,'defaulttextinterpreter','latex');
    xx = (-5:0.01:0);
    MitISEM_plot(mit2, 2, xx , [], GamMat);
    hold on
    yy = exp(kernel(xx'));
    yy = yy/max(yy);
    plot(xx,yy,'r')
    hold off
    xlabel('$$\varepsilon_{1}$$','FontSize', 12) % x-axis label
    ylabel('Approximation to the tail by $$q_{2}$$','FontSize', 12) % y-axis label
    plotTickLatex2D('FontSize',12);
end 

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','mit2','summary2')
end


%% QERMit 2:  MONTE CARLO VaR_mit and ES_mit (and their NSEs) ESTIMATION 
% use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
% to estiamte VaR and ES for theta and y (or alpha in eps)
tic    
for sim = 1:N_sim
    fprintf('VaR IS sim = %i.\n', sim);
    draw_mit = rmvgt2(M, mit2.mu, mit2.Sigma, mit2.df, mit2.p); 
    y_mit = predict_t_garch_noS(theta_mat, y_T, h_mat, H, draw_mit);
    PL_mit = fn_PL(y_mit);
    
    lnk_opt = duvt(draw_mit, nu_mat, H, true); %log density
    lnd_opt = dmvgt(draw_mit, mit2, true, GamMat);
    w_opt = exp(lnk_opt - lnd_opt)/M;
    [PL, ind] = sort(PL_mit);         
    w_opt = w_opt(ind,:);
    cum_w = cumsum(w_opt);
    ind_var = min(find(cum_w >= p_bar))-1; 
    VaR_mit(sim,1) = PL(ind_var);
    ES = (w_opt(1:ind_var)/sum(w_opt(1:ind_var))).*PL(1:ind_var);
    ES_mit(sim,1) = sum(ES(isfinite(ES)));    
    
    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_mit(sim,1), model, algo);  
end
time_mit(2,1) = toc/N_sim;

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','mit2','summary2','VaR_mit','ES_mit','time_mit')
end

y_mit = predict_t_garch_noS(theta_mat, y_T, h_mat, H, draw_mit);
PL_mit = fn_PL(y_mit);
mit_eff = sum(PL_mit <= mean(VaR_direct))/M;

if save_on
    name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'cont2','mit2','summary2','VaR_mit','ES_mit','time_mit','mit_eff')
end

if plot_on
    labels_in = {'naive','mit'};
    Boxplot_PMitISEM(VaR_direct, VaR_mit, ES_direct, ES_mit, model, algo, H, N_sim, true, labels_in);
end