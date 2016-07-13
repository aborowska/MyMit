%% Initialisation
% clc
clear all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));

v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014b)')
    v_new = 1;
else
    v_new = 0;
end

x = (0:0.00001:50)'+0.00001;
GamMat = gamma(x);

algo = 'Direct';
model = 'arch';
s = 'direct';
hp = 1;

plot_on = false;
print_on = false;
save_on = true;

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
N_sim = 100;
cont.mit.dfnc = 1;
plot_on = true;
print_on  = false;


kernel_init = @(a) - posterior_arch(a, data, S, true);
[mu, Sigma] = fn_initopt(kernel_init, mu_init);
mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',cont.mit.dfnc);
kernel = @(a) posterior_arch(a, data, S, true);

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
accept_direct = zeros(N_sim,1);


for sim = 1:N_sim
    fprintf('Direct sim = %i.\n', sim);
    [alpha_direct, accept_direct(sim,1)] = Mit_MH(M, kernel, mit_direct, GamMat);
    fprintf('Direct MH acceptance rate: %4.2f. \n',accept_direct(sim,1));
    draw_direct = alpha_direct;

    eps_direct = randn(M,1);
    h_direct = S*(1-draw_direct(:,1)) + draw_direct(:,1)*(y_T^2);
    y_direct = eps_direct.*sqrt(h_direct);
    
    ind_real = (imag(y_direct)==0);
    M_real = sum(ind_real); 
    y_direct = y_direct(ind_real,:);
    eps_direct = eps_direct(ind_real,:);

    PL_direct = sort(fn_PL(y_direct));
    VaR_direct(sim,1) = PL_direct(round(p_bar*M_real));
    ES_direct(sim,1) = mean(PL_direct(round(1:p_bar*M_real)));
    
    if (plot_on && (sim == N_sim))
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
%         title('Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$.')
        ylim([-15,15])
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
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

mean_VaR_direct = mean(VaR_direct);
mean_ES_direct = mean(ES_direct);

NSE_VaR_direct = std(VaR_direct);      
NSE_ES_direct = std(ES_direct);      

fprintf('Direct VaR (mean) estimate: %6.4f. \n',mean_VaR_direct);
fprintf('Direct VaR (mean) estimate: %6.4f. \n',mean_ES_direct);

fprintf('Direct NSE VaR estimate: %6.4f. \n',NSE_VaR_direct);
fprintf('Direct NSE ES estimate: %6.4f. \n',NSE_ES_direct);

if save_on
    gen_out2;
    save(['results/arch_',algo,'.mat'], 'mu_init', 'mit_direct', 'accept_direct','draw_direct', 'eps_direct', 'y_direct',  'p_bar', 'M', 'N_sim', 'VaR_direct', 'ES_direct');
end