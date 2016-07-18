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
estimation = 'mle'; % 'true' or 'mle'

% Artificial, white noise data 
T = 10000;
sigma2_true = 1;
y = sqrt(sigma2_true)*randn(T,1); 
y = y - mean(y);

% Control parameters for MitISEM
cont2 = MitISEM_Control;
cont2.mit.dfnc = 5;
cont2.mit.N = 10000;
cont2.mit.iter_max = 5;
cont2.df.range = [1,20];

sim = 1;
N_sim = 20;

VaR_mit = zeros(N_sim,1);
ES_mit = zeros(N_sim,1);
time_mit = zeros(2,1);

M = 10000; % number of draws for preliminary and IS computations

H = 1; % forecast horizon
p_bar = 0.01;
tail_shift = 1;

%% BIG DRAW
name =  ['results/PMitISEM/',model,'_Direct_',estimation,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
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
  
kernel_init = @(xx) - MLtarget_WN_hl(xx, sigma2_used, tail_shift*mean(VaR_direct));
kernel = @(xx) MLtarget_WN_hl(xx, sigma2_used, tail_shift*mean(VaR_direct));

tic
if (H < 100)
    [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);
else
    [mit2, summary2] = MitISEM_new(mit_hl, kernel, mu_hl, cont2, GamMat);
end
time_mit(1,1) = toc;

if (plot_on && (H == 1))
    figure(100) 
    set(gcf,'defaulttextinterpreter','latex');
    MitISEM_plot(mit2, 2, (-4:0.01:-1) , [], GamMat);
    xlabel('$$\sigma^2$$','FontSize', 12) % x-axis label
    ylabel('Approximation to the tail by $$q_{2}$$','FontSize', 12) % y-axis label
    plotTickLatex2D('FontSize',12);
end 

if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'mit2','summary2')
end

%% VaR with standard MitISEM
tic
for sim = 1:N_sim 
    fprintf('\nVaR IS iter: %d\n',sim)

    draw_opt = rmvgt2(M, mit2.mu, mit2.Sigma, mit2.df, mit2.p); 
    y_opt = sqrt(sigma2_used).*draw_opt;  
    PL = fn_PL(y_opt);
    kernel = @(xx) - 0.5*(log(2*pi) + log(sigma2_used)+ (xx.^2)/sigma2_used);

    lnk_opt = sum(kernel(draw_opt),2);
    lnd_opt = dmvgt(draw_opt, mit2, true, GamMat);
    w_opt = exp(lnk_opt - lnd_opt)/M;
    [PL, ind] = sort(PL);         
    w_opt = w_opt(ind,:);
    cum_w = cumsum(w_opt);
%     ind_var = min(find(cum_w >= p_bar))-1; 
    ind_var = min(find(cum_w > p_bar))-1; 
%     VaR_mit(sim,1) = PL(ind_var);
    VaR_mit(sim,1)  = (PL(ind_var+1) + PL(ind_var))/2; % intrapolate
    ES = (w_opt(1:ind_var)/sum(w_opt(1:ind_var))).*PL(1:ind_var);
    ES_mit(sim,1) = sum(ES(isfinite(ES)));     

    fprintf('IS 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_mit(sim,1), model, algo);  
end
time_mit(2,1) = toc/N_sim;


if save_on
    name = ['results/PMitISEM/',model,'_',algo,'_',estimation,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    save(name,'VaR_mit','ES_mit','mit2','summary2','time_mit')
end


labels_in = {'Direct','MitISEM'};
algo = [algo,'_',estimation];  
Boxplot_PMitISEM(VaR_direct, VaR_mit, ES_direct, ES_mit, model, algo, H, N_sim, true, labels_in);