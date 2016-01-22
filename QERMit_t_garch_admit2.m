%% Initialisation
clear all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x);

data = csvread('GSPC_ret_tgarch.csv');
data = 100*data;

T = size(data,1);
y_T = data(T);
S = var(data);

M = 10000;
N_sim = 20;

plot_on = false;
print_on = false;

AdMit_Control
cont.IS.opt = true;
cont.IS.scale = [0.8, 1.0, 1.2]; 
cont.IS.perc = [0.25, 0.30, 0.35];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1a.:
L = true;
hyper = 1;
% theta = [alpha, beta, mu, nu]
mu_init = [0.03, 0.9, 0.03, 6];

posterior_t_garch(mu_init, data, S, L, hyper, GamMat)

kernel_init = @(a) - posterior_t_garch(a, data, S, L, hyper, GamMat);
kernel = @(a) posterior_t_garch(a, data, S, L, hyper, GamMat);

[mit1, summary1] = AdMit(kernel_init, kernel, mu_init, cont, GamMat);
% save('results/t_garch_admit.mat','mit1','summary1');

%% >>>> not necessary when NSE estimation
% % for QERMit 2.
% [draw1, lnk1, ind_red1] = fn_rmvgt_robust(M, mit1, kernel);
% % draw1 = [draw1, trnd(draw1(:,4))]; % ERRORS ARE T!!
% save('results/t_garch_admit_draw1.mat','draw1', 'lnk1', 'ind_red1');


%% Figure draws
if plot_on
    figure(1)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'defaulttextinterpreter','latex');

    subplot(2,2,1)
    [f,xi] = ksdensity(draw1(:,1));
    plot(xi,f);
    title('$$\alpha$$')
    plotTickLatex2D;

    subplot(2,2,2)
    [f,xi] = ksdensity(draw1(:,2));
    plot(xi,f);
    title('$$\beta$$')
    plotTickLatex2D;

    subplot(2,2,3)
    [f,xi] = ksdensity(draw1(:,3));
    plot(xi,f);
    title('$$\mu$$')
    plotTickLatex2D;

    subplot(2,2,4)
    [f,xi] = ksdensity(draw1(:,4));
    plot(xi,f);
    title('$$\nu$$')
    plotTickLatex2D;

    suptitle('Smoothed kernel estimates')

    if print_on
        name = 'figures/t_garch_admit_draw1_admit.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.

%% QERMit 1b.:
% generate set opf draws of theta using independence MH with
% candiate from AdMit; then simulate returns based on the draw of theta 
[theta, accept] = Mit_MH(M+1000, kernel, mit1, GamMat);
fprintf('(AdMit) MH acceptance rate: %4.2f. \n',accept);
% save('results/t_garch_admit_theta.mat','theta', 'accept');
theta = theta(1001:M+1000,:);
% acf = [autocorr(theta(:,1),1), autocorr(theta(:,2),1), autocorr(theta(:,3),1), autocorr(theta(:,4),1)];
% acf = acf(2,:);

%% High loss, 10 days horizon
% approximate the high loss distribution of (theta,eps*) where eps*={eps_T+1,...,eps_T+hp}
p_bar = 0.01;
hp = 1; % prediction horizon 

h_T = volatility_t_garch(theta, data, S);
[y_hp, eps_hp] = predict_t_garch(theta, y_T, S, h_T, hp);
% get the preliminary 10-day-ahead 99% VaR estimate as the 100th of the ascendingly sorted percentage loss values
[PL_hp, ind] = sort(fn_PL(y_hp));
eps_hl = eps_hp(ind,:); 
eps_hl = eps_hl(1:p_bar*M,:);
theta_hl = theta(ind,:);
theta_hl = theta_hl(1:p_bar*M,:);

VaR_prelim = PL_hp(p_bar*M);
fprintf('hp = %i, y_T = %4.2f, VaR_prelim = %4.5f. \n', hp, y_T, VaR_prelim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1c.:
% get mit approximation of the conditional joint density of
% parameters and future returns given the returns are below VaR_prelim
% approximation of the joint high loss distribution
% here: not future returns but future disturbances  (varepsilons)

%% High loss density approximation
L = true;
kernel_init = @(a) - posterior_t_garch_hl(a, data, S, VaR_prelim, L, hyper, GamMat);
kernel = @(a) posterior_t_garch_hl(a, data, S, VaR_prelim, L, hyper, GamMat);
% theta = [alpha, beta, mu, nu, eps_T1, ..., eps_Thp] <-- size 14 
mu_init_hl = [0.07, 0.93, 0.05, 8.5, -4*ones(1, hp)];
% mu_init_hl = [0.0628    0.9311    0.0491    8.2365   -2.9356];
% mu_hl = median([theta_hl, eps_hl],1);
% kernel(mu_hl)

%% IS estimation of the initial mixture component 
draw_hl = [theta_hl, eps_hl];
% lnk = fn_lkereval(kernel, dupa);
lnk = kernel(draw_hl);
ind = find(lnk~=-Inf);
lnk = lnk(ind,:);
draw_hl = draw_hl(ind,:);
% lnd = dmvgt(theta,mit_init,true);
% w = fn_ISwgts(lnk, lnd, norm);
w = lnk/sum(lnk);
[mu_hl, Sigma_hl] = fn_muSigma(draw_hl, w);

mit_hl.mu = mu_hl;
mit_hl.Sigma = Sigma_hl;
mit_hl.df = cont.dfnc;
mit_hl.p = 1;

cont2 = cont;
cont2.IS.opt = true; 
cont2.Hmax = 2;
cont2.dfnc = 5;
cont2.IS.scale = [1.0]; 
cont2.IS.perc = [0.25];

[mit2, summary2] = AdMit(mit_hl, kernel, mu_init_hl, cont2, GamMat);
% [mit2, summary] = AdMit(kernel_init, kernel, mu_hl, cont, GamMat);
save(['results/t_garch_admit_hl_',num2str(hp),'.mat'],'mit2','summary2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 2:
% use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
% to estiamte VaR and ES for theta and y (or alpha in eps)

plot_on = false;
print_on = false;
VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONTE CARLO NSE AND RNE ESTIMATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sim = 5:N_sim
    fprintf('NSE sim = %i.\n', sim);
    kernel = @(a) posterior_t_garch(a, data, S, L, hyper, GamMat);
    [draw1, lnk1, ind_red1] = fn_rmvgt_robust(M/2, mit1, kernel);
    eps1 = zeros(M/2, hp);
    for hh = 1:hp
       eps1(:,hh) = trnd(draw1(:,4)); % ERRORS ARE iid T!!
    end
    draw1_eps1 = [draw1(1:M/2,:), eps1];
    
    kernel = @(a) posterior_t_garch_hl(a, data, S, VaR_prelim, L, hyper, GamMat);
    [draw2, lnk2, ind_red2] = fn_rmvgt_robust(M/2, mit2, kernel);
    % col 1-4 parameters, col 5:5+hp-1 eps

    draw_opt = [draw1_eps1; draw2];

    %% IS weights
    L = true;
    kernel = @(a) posterior_t_garch_whole(a, data, S, L, hyper, GamMat);
    lnk = kernel(draw_opt);

    % exp_lnd1 = 0.5*normpdf(draw_opt(:,2)).*dmvgt(draw_opt(:,1:4),mit1,false);
    % ep = duvt(draw_opt(:,5:5+hp-1), nu, hp, false);  % false: exp dnesitydraw_hl
    exp_lnd1 = 0.5*sum(duvt(draw_opt(:,5:5+hp-1), draw_opt(:,4), hp, false),2).*dmvgt(draw_opt(:,1:4), mit1, false, GamMat);
    exp_lnd2 = 0.5*dmvgt(draw_opt, mit2, false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd = log(exp_lnd);

    w_opt = fn_ISwgts(lnk, lnd, false);

    %% VaR and ES IS estimates
    % % h_opt = S*(1-draw_opt(:,1)) + draw_opt(:,1)*(y_T^2);
    % % y_opt = draw_opt(:,2).*sqrt(h_opt);
    h_T = volatility_t_garch(draw_opt(:,1:4), data, S);
    [y_opt, ~] = predict_t_garch(draw_opt(:,1:4), y_T, S, h_T, hp, draw_opt(:,5:5+hp-1));
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);
    
    PL_opt= fn_PL(y_opt);

    if plot_on    
        figure(8)
        set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
        set(gcf,'defaulttextinterpreter','latex');
        hold on
        plot(sort(PL_opt))
        pos =  max(find(sort(PL_opt)<=VaR_IS(sim,1)));
        scatter(pos, VaR_IS(sim,1),'MarkerFaceColor','red')
        hold off
    %     subplot(3,1,1)
    %     subplot(3,1,2)
    %     subplot(3,1,3)
        title('Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$.')
        set(gca,'TickLabelInterpreter','latex')

        if print_on
            name = 'figures/t_garch_predict_admit.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    end
    fprintf('(AdMit) IS VAR estimate: %6.4f. \n',VaR_IS(sim,1));
    fprintf('(AdMit) IS ES estimate: %6.4f. \n',ES_IS(sim,1));
    save(['results/t_garch_admit_',num2str(hp),'_sim.mat'],'draw_opt','lnk','lnk2','VaR_IS','ES_IS');

end

mean_VaR_IS = mean(VaR_IS(VaR_IS<0));
mean_ES_IS = mean(ES_IS(ES_IS<0));

NSE_VaR_IS = std(VaR_IS(VaR_IS<0));              
NSE_ES_IS = std(ES_IS(ES_IS<0));                

fprintf('AdMit IS VAR (mean) estimate: %6.4f. \n',mean_VaR_IS);
fprintf('AdMit IS VAR (mean) estimate: %6.4f. \n',mean_VaR_IS);

fprintf('AdMit IS NSE VaR estimate: %6.4f. \n',NSE_VaR_IS);
fprintf('AdMit IS NSE ES estimate: %6.4f. \n',NSE_ES_IS);

cont = cont2;
model = 't_garch';
s='admit';
gen_out

clear GamMat x fig
save(['results/t_garch_admit_',num2str(hp),'.mat']);


x = (0:0.00001:50)'+0.00001;
GamMat = gamma(x);
 