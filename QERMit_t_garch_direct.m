%% Initialisation
clear all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x = (0:0.00001:100)'+0.00001;
GamMat = gamma(x);

plot_on = false;
print_on = false;

data = csvread('GSPC_ret_tgarch.csv');
data = 100*data;
T = size(data,1);
y_T = data(T);
S = var(data);

M = 10000;
N_sim = 20;

L = true;
hyper = 1;
% theta = [alpha, beta, mu, nu]
% mu_init = [0.065, 0.93, 0.05, 8];
mu_init = [0.065 0.93 0.048 8.4];

% mu = [0.065 0.93 0.048 8.4];
% Sigma = [8.18505622977260e-05 -9.24538074257755e-05 1.09593583441952e-05 0.000599696689595632 -9.24538074257755e-05 0.000107858596235363 -1.14723061801121e-05 -0.000355774488666053 1.09593583441952e-05 -1.14723061801121e-05 0.000243547476400608 -0.000444051744593874 0.000599696689595632 -0.000355774488666053 -0.000444051744593874 0.790171013159574];
df = 10;
% df = 10 ==> Acceptance rate 76
% df = 1 ==> Acceptance rate 51

kernel_init = @(a) - posterior_t_garch(a, data, S, L, hyper, GamMat);
[mu, Sigma] = fn_initopt(kernel_init, mu_init);
mit_direct = struct('mu',mu,'Sigma',Sigma,'df',df,'p',1);
kernel = @(a) posterior_t_garch(a, data, S, L, hyper, GamMat);


VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);

hp = 1;
p_bar = 0.01;

[theta_direct, accept_direct] = Mit_MH(M+1000, kernel, mit_direct, GamMat);
fprintf('Direct MH acceptance rate: %4.2f. \n',accept_direct);
theta_direct = theta_direct(1001:M+1000,:);

for sim = 2:N_sim   
    [theta_direct, accept_direct] = Mit_MH(M+1000, kernel, mit_direct, GamMat);
    fprintf('Direct MH acceptance rate: %4.2f. \n',accept_direct);
    theta_direct = theta_direct(1001:M+1000,:);

    eps_direct = zeros(M, hp);
    for hh = 1:hp
       eps_direct(:,hh) = trnd(theta_direct(:,4)); % ERRORS ARE iid T!!
    end
    draw_direct = [theta_direct, eps_direct];

    h_T_direct = volatility_t_garch(draw_direct(:,1:4), data, S);
    [y_direct, ~] = predict_t_garch(draw_direct(:,1:4), y_T, S, h_T_direct, hp, draw_direct(:,5:5+hp-1));
    PL_direct = sort(fn_PL(y_direct));
    VaR_direct(sim,1) = PL_direct(p_bar*M);
    ES_direct(sim,1) = mean(PL_direct(1:p_bar*M));
    
     
    if plot_on    
        figure(8)
        set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
        set(gcf,'defaulttextinterpreter','latex');
        hold on
        plot(PL_direct)
        pos =  max(find(PL_direct<=VaR_direct(sim,1)));
        scatter(pos, VaR_direct(sim,1),'MarkerFaceColor','red')
        hold off
    %     subplot(3,1,1)
    %     subplot(3,1,2)
    %     subplot(3,1,3)
        title('Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$.')
        set(gca,'TickLabelInterpreter','latex')

        if print_on
            name = 'figures/t_garch_predict_direct.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    end
    fprintf('Direct VAR estimate: %6.4f. \n',VaR_direct(sim,1));
    fprintf('Direct ES estimate: %6.4f. \n',ES_direct(sim,1));
end

mean_VaR_direct = mean(VaR_direct);
mean_ES_direct = mean(ES_direct);

NSE_VaR_direct = std(VaR_direct);      
NSE_ES_direct = std(ES_direct);        


fprintf('Direct VAR (mean) estimate: %6.4f. \n',mean_VaR_direct);
fprintf('Direct ES (mean) estimate: %6.4f. \n',mean_ES_direct);

fprintf('Direct NSE VaR estimate: %6.4f. \n',NSE_VaR_direct);
fprintf('Direct NSE ES estimate: %6.4f. \n',NSE_ES_direct);

model = 't_garch';
s='direct';
gen_out


clear GamMat x
save(['results/t_garch_direct_',num2str(hp),'.mat']);