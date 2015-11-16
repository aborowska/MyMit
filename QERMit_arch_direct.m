%% Initialisation
% clc
clear all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));


x = (0:0.00001:50)'+0.00001;
GamMat = gamma(x);


data = csvread('GSPC_ret.csv');
% data = csvread('GSPC_ret_tgarch.csv');
data = 100*data;

% QERMit ARCH 
ind_arch = find(data<=-5.5, 1, 'last' );
data = data(1:ind_arch,1);
data = data - mean(data);

T = length(data);
y_T = data(T);
S = var(data);

p_bar = 0.01; % p_bar = 1-alpha, 100alpha% VaR

mu_init = 0.03;
mu_hl = [0.15, -3];
M = 10000;
N_sim = 20;
df = 10;
plot_on = false;
print_on  = false;


L = true;
kernel_init = @(a) - posterior_arch(a, data, S, L);
[mu, Sigma] = fn_initopt(kernel_init, mu_init);
mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',df);
kernel = @(a) posterior_arch(a, data, S, L);
VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);

for sim = 1:N_sim
    [alpha_direct, accept_direct] = Mit_MH(M, kernel, mit_direct, GamMat);
    fprintf('Direct MH acceptance rate: %4.2f. \n',accept_direct);
    draw_direct = alpha_direct;

    h_direct = S*(1-draw_direct(:,1)) + draw_direct(:,1)*(y_T^2);
    y_direct = randn(M,1).*sqrt(h_direct);
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
            name = 'figures/arch_predict_direct.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    end
    fprintf('Direct VAR estimate: %6.4f. \n',VaR_direct(sim,1));
    fprintf('Direct ES estimate: %6.4f. \n',ES_direct(sim,1));
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

mean_VaR_direct = mean(VaR_direct);
mean_ES_direct = mean(ES_direct);

NSE_VaR_direct = std(VaR_direct);      
NSE_ES_direct = std(ES_direct);      

fprintf('Direct VAR (mean) estimate: %6.4f. \n',mean_VaR_direct);
fprintf('Direct VAR (mean) estimate: %6.4f. \n',mean_ES_direct);

fprintf('Direct NSE VaR estimate: %6.4f. \n',NSE_VaR_direct);
fprintf('Direct NSE ES estimate: %6.4f. \n',NSE_ES_direct);

model = 'arch';
s = 'direct';
hp = 1;
gen_out;

clear GamMat x fig
save(['results/arch_direct_',num2str(hp),'_df_',num2str(df),'.mat']);