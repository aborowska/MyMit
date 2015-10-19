%% Initialisation
clear all
addpath('include/');
addpath('include/MEX/');

% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

% data = csvread('GSPC_ret_tgarch.csv');
% data = 100*data;
y = csvread('GSPC_ret.csv');
y = 100*y;
ind_arch = find(y<=-5.5, 1, 'last' );
y = y(1:ind_arch,1);
y = y - mean(y);
data=y;

T = size(y,1);


T = size(data,1);
y_T = data(T);
S = var(data);

% M = 10000;
% N_sim = 20;
M = 10000;
N_sim = 10;
SS = 10;

model = 't_garch';
    
save_on = false;
plot_on = false;
print_on = false;

if plot_on
	figure(1)
	set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
    xx = linspace(1998,2008,T);
    plot(xx, data)
%     set(gca,'XTickLabel',{'1998',' ','1999',' ','2000',' '})
	title('S$$\&$$P 500 log-returns $$y$$')
%    plotTickLatex2D;
	set(gca,'TickLabelInterpreter','latex')
	if print_on
		name = 'figures/t_garch_data.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end

MitISEM_Control
cont.mit.dfnc = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1a.:
L = true;
hyper = 1;
% theta = [alpha, beta, mu, nu]
mu_init = [0.03, 0.9, 0.03, 6];

posterior_t_garch(mu_init, data, S, L, hyper, GamMat)

kernel_init = @(a) - posterior_t_garch(a, data, S, L, hyper, GamMat);
kernel = @(a) posterior_t_garch(a, data, S, L, hyper, GamMat);

[mit1, summary1] = MitISEM(kernel_init, kernel, mu_init, cont, GamMat);
if save_on
    save('results/t_garch_mitisem.mat','mit1','summary1');
end

%% >>>> not necessary when NSE estimation
% % for QERMit 2.
% [draw1, lnk1, ind_red1] = fn_rmvgt_robust(M, mit1, kernel);
% % draw1 = [draw1, trnd(draw1(:,4))]; % ERRORS ARE T!!
% save('results/t_garch_mitisem_draw1.mat','draw1', 'lnk1', 'ind_red1');

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
        name = 'figures/t_garch_mitisem_draw1.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.

%% QERMit 1b.:
% generate set opf draws of theta using independence MH with
% candiate from MitISEM; then simulate returns based on the draw of theta 

VaR_prelim = zeros(SS,1);
accept = zeros(SS,1);
for sim = 1:SS  
    [theta, accept(sim,1)] = Mit_MH(M+1000, kernel, mit1, GamMat);
    fprintf('(MitISEM) MH acceptance rate: %6.4f. \n',accept(sim));
    if save_on
        save('results/t_garch_mitisem_theta.mat','theta', 'accept');
    end

    theta = theta(1001:M+1000,:);
    % acf = [autocorr(theta(:,1),1), autocorr(theta(:,2),1), autocorr(theta(:,3),1), autocorr(theta(:,4),1)];
    % acf = acf(2,:);

    %% High loss, 10 days horizon
    % approximate the high loss distribution of (theta,eps*) where eps*={eps_T+1,...,eps_T+hp}
    p_bar = 0.01*sim;
    hp = 1; % prediction horizon 

    h_T = volatility_t_garch(theta, data, S);
    [y_hp, eps_hp] = predict_t_garch(theta, y_T, S, h_T, hp);
    % get the preliminary 10-day-ahead 99% VaR estimate as the 100th of the ascendingly sorted percentage loss values
    [PL_hp, ind] = sort(fn_PL(y_hp));
    eps_hl = eps_hp(ind,:); 
    eps_hl = eps_hl(1:p_bar*M,:);
    theta_hl = theta(ind,:);
    theta_hl = theta_hl(1:p_bar*M,:);

    VaR_prelim(sim) = PL_hp(round(p_bar*M));
    fprintf('p_bar = %4.2f, y_T = %4.2f, VaR_prelim = %4.5f. \n', p_bar, y_T, VaR_prelim(sim))
end

figure(158)
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
plot(0.01:0.01:0.1,VaR_prelim)
xlabel('$$p_{bar}$$')
ylabel('$$VaR_{prelim}(p_{bar})$$')
% set(ax,'XTickLabel',{'0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08','0.09','0.10'});
% ax.XTickLabel = {'0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08','0.09','0.10'};
plotTickLatex2D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1c.:
% get mit approximation of the conditional joint density of
% parameters and future returns given the returns are below VaR_prelim
% approximation of the joint high loss distribution
% here: not future returns but future disturbances  (varepsilons)
p_bar = 0.01
VaR_prel = VaR_prelim;
VaR_prelim = VaR_prelim(1,1);

%% High loss density approximation
L = true;
kernel_init = @(a) - posterior_t_garch_hl(a, data, S, VaR_prelim, L, hyper, GamMat);
kernel = @(a) posterior_t_garch_hl(a, data, S, VaR_prelim, L, hyper, GamMat);
% theta = [alpha, beta, mu, nu, eps_T1, ..., eps_Thp] <-- size 14 
mu_init_hl = [0.07, 0.93, 0.05, 8.5, -1.5*ones(1, hp)];
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
mit_hl.df = 5;
mit_hl.p = 1;

cont2 = cont;  
cont2.mit.Hmax = 2;
[mit2, summary2] = MitISEM(mit_hl, kernel, mu_init_hl, cont2, GamMat);
%[mit2, summary] = MitISEM(kernel_init, kernel, mu_hl, cont, GamMat);
if save_on
    save(['results/t_garch_mitisem_hl_',num2str(hp),'.mat'],'mit2','summary2');
end

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
for sim = 1:N_sim
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
% lnk_opt = lnk1 + duvt(eps1, draw1(:,4), 1, true);
% lnk_opt = [lnk_opt; lnk2];
 
    % exp_lnd1 = 0.5*normpdf(draw_opt(:,2)).*dmvgt(draw_opt(:,1:4),mit1,false);
    % ep = duvt(draw_opt(:,5:5+hp-1), nu, hp, false);  % false: exp dnesity
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
% [PL_opt_h1, ind] = sort(PL_opt);
% w_opt_h1 = w_opt(ind)/sum(w_opt);
% cum_w = cumsum(w_opt_h1);
% 
% lnd_opt_h1 = lnd(ind);
% lnk_opt_h1 = lnk(ind);
% 
% 
% figure(11)
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% set(gcf,'defaulttextinterpreter','latex');
% subplot(2,2,1); plot(cum_w); hold on; plot(0.01*ones(M),'r'); hold off; title('cum w');
% subplot(2,2,2); plot(lnk_opt_h1);   title('lnk opt h1');
% subplot(2,2,3); plot(w_opt_h1);     title('w opt h1');
% subplot(2,2,4); plot(lnd_opt_h1);   title('lnd opt h1');


    if plot_on    
        figure(8)
        set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
        set(gcf,'defaulttextinterpreter','latex');
        hold on
        plot(PL_opt_h1)
        pos = max(find(sort(PL_opt)<=VaR_IS(sim,1)));
        scatter(pos, VaR_IS(sim,1),'MarkerEdgeColor','red','MarkerFaceColor','red')
        pos = max(find(PL_opt_h1<=VaR_prelim));
        scatter(pos, VaR_prelim,'MarkerEdgeColor','green','MarkerFaceColor','green')
       
        hold off
    %     subplot(3,1,1)
    %     subplot(3,1,2)
    %     subplot(3,1,3)
        title(['Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$. Model: ',model,'.'])
        set(gca,'TickLabelInterpreter','latex')

        if print_on
            name = 'figures/t_garch_predict.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    end
    if save_on
        save(['results/t_garch_mitisem_',num2str(hp),'_sim.mat'],'draw_opt','lnk','lnk2','VaR_IS','ES_IS');
    end
    
    fprintf('MitISEM IS VAR estimate: %6.4f. \n',VaR_IS(sim,1));
    fprintf('MitISEM IS ES estimate: %6.4f. \n',ES_IS(sim,1));

end

mean_VaR_IS = mean(VaR_IS(VaR_IS<0));
mean_ES_IS = mean(ES_IS(ES_IS<0));

NSE_VaR_IS = std(VaR_IS(VaR_IS<0));              
NSE_ES_IS = std(ES_IS(ES_IS<0));                
           

fprintf('MitISEM IS VAR (mean) estimate: %6.4f. \n',mean_VaR_IS);
fprintf('MitISEM IS ES (mean) estimate: %6.4f. \n',mean_ES_IS);

fprintf('MitISEM IS NSE VaR estimate: %6.4f. \n',NSE_VaR_IS);
fprintf('MitISEM IS NSE ES estimate: %6.4f. \n',NSE_ES_IS);


if save_on
    cont = cont2;
    
    s='mitisem';
    gen_out
        
    clear GamMat x fig
    save(['results/t_garch_mitisem_',num2str(hp),'.mat']);
end