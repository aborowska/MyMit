%% Initialisation
clear all
addpath(genpath('include/'));

% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s); 

v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014b)')
    v_new = 1;
else
    v_new = 0;
end

algo = 'Direct';
model = 't_garch';
p_bar = 0.01;
hp = 1;

x = (0:0.00001:100)'+0.00001;
GamMat = gamma(x);

plot_on = false;
print_on = false;
save_on = true;

data = csvread('GSPC_ret_tgarch.csv');
data = 100*data;
T = size(data,1);
y_T = data(T);
S = var(data);

M = 10000;
N_sim = 100;

% L = true;
% hyper = 1;
% theta = [alpha, beta, mu, nu]
% mu_init = [0.065, 0.93, 0.05, 8];
mu_init = [0.065 0.93 0.048 8.4];


% mu = [0.065 0.93 0.048 8.4];
% Sigma = [8.18505622977260e-05 -9.24538074257755e-05 1.09593583441952e-05 0.000599696689595632 -9.24538074257755e-05 0.000107858596235363 -1.14723061801121e-05 -0.000355774488666053 1.09593583441952e-05 -1.14723061801121e-05 0.000243547476400608 -0.000444051744593874 0.000599696689595632 -0.000355774488666053 -0.000444051744593874 0.790171013159574];
cont.mit.dfnc = 5;
% df = 10 ==> Acceptance rate 76
% df = 1 ==> Acceptance rate 51

% kernel_init = @(a) - posterior_t_garch(a, data, S, L, hyper, GamMat);
kernel_init = @(a) - posterior_t_garch_mex(a, data, S, GamMat);
% kernel = @(a) posterior_t_garch(a, data, S, L, hyper, GamMat);
kernel = @(a) posterior_t_garch_mex(a, data, S, GamMat);
[mu, Sigma] = fn_initopt(kernel_init, mu_init);
mit_direct = struct('mu',mu,'Sigma',Sigma,'df',cont.mit.dfnc,'p',1);

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
accept_direct = zeros(N_sim,1);


% [theta_direct, accept_direct] = Mit_MH(M+1000, kernel, mit_direct, GamMat);
% fprintf('Direct MH acceptance rate: %4.2f. \n',accept_direct);
% theta_direct = theta_direct(1001:M+1000,:);

for sim = 1:N_sim   
    fprintf('Direct sim = %i.\n', sim);
    [theta_direct, accept_direct(sim,:)] = Mit_MH(M+1000, kernel, mit_direct, GamMat);
    fprintf('MH acceptance rate: %6.4f (%s, %s). \n', accept_direct(sim,:), model, algo);
    theta_direct = theta_direct(1001:M+1000,:);

%     eps_direct = zeros(M, hp);
%     for hh = 1:hp
%        eps_direct(:,hh) = trnd(theta_direct(:,4)); % ERRORS ARE iid T!!
%     end
    
    h_T_direct = volatility_t_garch_mex(theta_direct, data, S);
    [y_direct, eps_direct] = predict_t_garch(theta_direct, y_T, S, h_T_direct, hp);

    ind_real = (imag(y_direct)==0);
    M_real = sum(ind_real); 
    y_direct = y_direct(ind_real,:);
    eps_direct = eps_direct(ind_real,:);

    PL_direct = sort(fn_PL(y_direct));
    VaR_direct(sim,1) = PL_direct(round(p_bar*M_real));
    ES_direct(sim,1) = mean(PL_direct(round(1:p_bar*M_real)));
  
    fprintf('p_bar = %4.2f, VaR_direct = %4.5f. \n', p_bar, VaR_direct(sim,1))
    fprintf('p_bar = %4.2f, NSE VaR_direct = %4.5f. \n', p_bar, std(VaR_direct(VaR_direct<0,1)))
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_direct(sim,1), model, algo);  
end

mean_VaR_direct = mean(VaR_direct);
mean_ES_direct = mean(ES_direct);

NSE_VaR_direct = std(VaR_direct);      
NSE_ES_direct = std(ES_direct);        


fprintf('Direct VAR (mean) estimate: %6.4f. \n',mean_VaR_direct);
fprintf('Direct ES (mean) estimate: %6.4f. \n',mean_ES_direct);

fprintf('Direct NSE VaR estimate: %6.4f. \n',NSE_VaR_direct);
fprintf('Direct NSE ES estimate: %6.4f. \n',NSE_ES_direct);

% model = 't_garch';
% s='direct';
% gen_out
if plot_on    
    figure(8)
    if v_new
        set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    else
        set(gcf,'units','normalized','outerposition',[0 0 0.5 0.75]);
    end
    set(gcf,'defaulttextinterpreter','latex');
    hold on
    plot(PL_direct)
    pos =  max(find(PL_direct <= mean(VaR_direct)));
    scatter(pos, mean(VaR_direct),'MarkerEdgeColor','red','MarkerFaceColor','red')
    hold off
%     title('Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$.')
    ylim([-10,10])
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    if print_on
        name = 'figures/t_garch_predict_direct.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
end

if save_on
    save(['results/t_garch_',algo,'.mat'], 'mu_init', 'mit_direct', 'accept_direct','theta_direct', 'eps_direct', 'h_T_direct', 'y_direct',  'p_bar', 'M', 'N_sim', 'VaR_direct', 'ES_direct');
    gen_out2;
end  

% clear GamMat x
% save(['results/t_garch_direct_',num2str(hp),'.mat']);