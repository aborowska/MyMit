%% Initialisation
% clc
clear all
close all
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

x = (0:0.00001:50)' + 0.00001;
GamMat = gamma(x);

model = 'arch';
data = csvread('GSPC_ret.csv');
data = 100*data;

% QERMit ARCH 
ind_arch = find(data<=-5.5, 1, 'last' );
data = data(1:ind_arch,1);
data = data - mean(data);

T = length(data);
y_T = data(T);
S = var(data);
f_stdev = @(aa) sqrt(S+(y_T^2-S)*aa); % then the one-day-ahead forecast given alpha is y_T1 = f_stdev(alpha).*randn
        
mu_init = 0.03;
% mu_hl = [0.15, -3]; % <-- for p_bar = 0.01 
M = 10000;
N_sim1 = 10;
N_sim = 10;

hp_on = true; % explicitly approximate the high PROFIT density, i.e. the one for the COMPLIMENT of the high LOSS region
if hp_on
    model = [model,'_hp'];
end

plot_on = false;
print_on  = true;
plot_on2 = true;
save_on = true;

MitISEM_Control
cont.mit.dfnc = 5;
cont.mit.N = 1000;

cont2 = cont;
% cont2.mit.Hmax = 10;
% cont2.mit.dfnc = 5;
% cont2.mit.CV_tol = 0.1;
cont2.df.range = [1, 50];
if hp_on
    cont3 = cont;
    cont3.mit.Hmax = 5;
    cont3.mit.N = 10000;
    cont3.mit.CV_tol = 0.01;
end


% p_bar = 0.05; % p_bar = 1-alpha, 100alpha% VaR
P_bars = [0.01, 0.02, 0.05, 0.1, 0.5];

    VaR_prelim = zeros(N_sim1,1);
    VaR_prelim_MC = zeros(N_sim1,length(P_bars));
    ES_prelim = zeros(N_sim1,length(P_bars));
    accept = zeros(N_sim1,length(P_bars));
    
    VaR_IS = zeros(N_sim,length(P_bars));
    ES_IS = zeros(N_sim,length(P_bars));
    hl_w = zeros(N_sim,length(P_bars)); % Sum of weights for high losses
    hp_w = zeros(N_sim,length(P_bars)); % Sum of weights for high profits

for p_bar = P_bars
    fprintf('\np_bar: %4.2f\n',p_bar);
    %% QERMit 1a.: 
    kernel_init = @(a) - posterior_arch(a, data, S, true);
    kernel = @(a) posterior_arch(a, data, S, true);
%     [mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);

    for sim = 1:N_sim1
        [mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);

        %% QERMit 1b.:
        % generate set of draws of alpha using independence MH with candidate from MitISEM; 
        % then simulate returns from normal with variance based on the draw of alpha 

        [alpha1, accept(sim,P_bars==p_bar)] = Mit_MH(M+1000, kernel, mit1, GamMat);
        fprintf('(%s) MH acceptance rate: %4.2f. \n', model, accept(sim,P_bars==p_bar));
        alpha1 = alpha1(1001:M+1000);

        stdev = f_stdev(alpha1);
        eps1 = randn(M,1);
        y_T1 = stdev.*eps1;

        % get the preliminary VaR estimate as the 100th of the ascendingly sorted percentage loss values
        [PL_T1, ind] = sort(fn_PL(y_T1));
        VaR_prelim(sim,1) = PL_T1(p_bar*M); % VaR_prelim = 0; VaR_prelim = Inf;
        ES_prelim(sim,P_bars==p_bar) = mean(PL_T1(1:p_bar*M));    
        fprintf('(%s) Preliminary 100*%4.2f%% VaR estimate: %6.4f. \n', model, p_bar, VaR_prelim(sim,1));
    end

    % take one value of VaR_prelim to construct mit2
    VaR_prelim_MC(:,P_bars==p_bar) = VaR_prelim;
%     VaR_prelim = VaR_prelim_MC(N_sim1,1);       % the last one
    VaR_prelim = mean(VaR_prelim);              % the mean
    
    
    arch_plot0; % The posterior density and the approximation to the posterior density

    %% high loss
    % get mit approximation of the conditional joint density of
    % parameters and future returns given the returns are below VaR_prelim
    % approximation of the joint high loss distribution
    % here: not future returns but future disturbances  (varepsilons)
    
    arch_plot1; % Fig. 6.1: the high loss subspace
    
    %% QERMit 1c.: 
    kernel_init = @(a) - posterior_arch_hl(a, data, S, VaR_prelim, true);
    kernel = @(a) posterior_arch_hl(a, data, S, VaR_prelim, true);

    % Choose the starting point (mu_hl) for the constuction of the approximaton 
    % to the high loss (hl) region density
    alpha1_hl = alpha1(ind); 
    eps_hl = eps1(ind); 
    draw_hl = [alpha1(ind), eps1(ind)];
    mu_hl = draw_hl(max(find(PL_T1 < VaR_prelim)),:);    
   
    % [mit2, summary2] = MitISEM(kernel_init, kernel, mu_hl, cont2, GamMat);
    [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);

    arch_plot2; %  The high loss density and the approximation to the high loss density
    
%% High profit    
    if hp_on 
        kernel_init = @(x) - posterior_arch_hp(x, data, S, VaR_prelim, true);
        kernel = @(x) posterior_arch_hp(x, data, S, VaR_prelim, true);

        mu_hp = mean(draw_hl((PL_T1 > VaR_prelim),:));
        % fn_PL(f_stdev(mu_hp(1,1))*mu_hl(1,2))
        [mit3, summary3] = MitISEM_new(kernel_init, kernel, mu_hp, cont3, GamMat);
        plot_HighProfit;
    end

    %% QERMit 2: 
    % use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
    % to estimate VaR and ES for theta and y (or alpha in eps)

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MONTE CARLO VaR_IS and ES_IS (and their NSEs) ESTIMATION 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for sim = 1:N_sim
        if hp_on
            kernel = @(x) posterior_arch_hp(x, data, S, VaR_prelim, true);
            [draw1, lnk1, ~] = fn_rmvgt_robust(M, mit3, kernel); % Robust drawing from the multivariate mixture of t        
        else
            kernel = @(a) posterior_arch(a, data, S, true);
            [draw1, lnk1, ~] = fn_rmvgt_robust(M, mit1, kernel);
            eps1 = randn(M,1);
            draw1 = [draw1, eps1];
            lnk1 = lnk1 - 0.5*(log(2*pi) + eps1.^2);
        end
        kernel = @(a) posterior_arch_hl(a, data, S, VaR_prelim, true);
        [draw2, lnk2, ~] = fn_rmvgt_robust(M, mit2, kernel);
        
        draw_opt = [draw1; draw2];
        
        arch_plot3; %  Draws from the optimal importance density, the joint density and the approximation to the optimal posterior

        %% IS weights
    %     kernel = @(a) posterior_arch_whole(a, data, S, true);
    %     lnk = kernel(draw_opt);
        lnk = [lnk1; lnk2];
        if hp_on      
            exp_lnd1 = 0.5*dmvgt(draw_opt,mit3,false, GamMat);
        else
            exp_lnd1 = 0.5*normpdf(draw_opt(:,2)).*dmvgt(draw_opt(:,1), mit1, false, GamMat);
        end
        exp_lnd2 = 0.5*dmvgt(draw_opt, mit2, false, GamMat);
        exp_lnd = exp_lnd1 + exp_lnd2;
        lnd = log(exp_lnd);
        w_opt =  fn_ISwgts(lnk, lnd, false);

        hl_w(sim,P_bars==p_bar) = sum( w_opt(f_stdev(draw_opt(:,1)).*draw_opt(:,2)<VaR_prelim,:)/sum(w_opt) );    
        hp_w(sim,P_bars==p_bar) = sum( w_opt(f_stdev(draw_opt(:,1)).*draw_opt(:,2)>VaR_prelim,:)/sum(w_opt) );   
        
        fprintf('Sum of weights for high losses: %6.4f and for high profits: %6.4f.\n', hl_w(sim,P_bars==p_bar), hp_w(sim,P_bars==p_bar));

        %% VaR and ES IS estimates
        y_opt = f_stdev(draw_opt(:,1)).*draw_opt(:,2);
        dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
        IS_estim = fn_PL(dens, 1);
        VaR_IS(sim,P_bars==p_bar) = IS_estim(1,1);
        ES_IS(sim,P_bars==p_bar) = IS_estim(1,2);
        PL_opt = fn_PL(y_opt);
    % [PL_opt_h1, ind] = sort(PL_opt);
    % w_opt_h1 = w_opt(ind)/sum(w_opt);
    % cum_w = cumsum(w_opt_h1);
    % lnd_opt_h1 = lnd(ind);
    % lnk_opt_h1 = lnk(ind);  
    % 
    % figure(11)
    % subplot(2,2,1); plot(cum_w); hold on; plot(0.01*ones(M),'r'); hold off; title('cum w');
    % subplot(2,2,2); plot(lnk_opt_h1);   title('lnk opt h1');
    % subplot(2,2,3); plot(w_opt_h1);     title('w opt h1');
    % subplot(2,2,4); plot(lnd_opt_h1);   title('lnd opt h1');

        arch_plot4; % The sorted future profit/losses 

        fprintf('(%s) IS 100*%4.2f%% VAR estimate: %6.4f. \n', model, p_bar, VaR_IS(sim,P_bars==p_bar));
        fprintf('(%s) IS 100*%4.2f%%ES estimate: %6.4f. \n', model, p_bar, ES_IS(sim,P_bars==p_bar));  
    end

    mean_VaR_prelim = mean(VaR_prelim_MC(:,P_bars==p_bar));
    mean_ES_prelim = mean(ES_prelim(:,P_bars==p_bar));

    NSE_VaR_prelim = std(VaR_prelim_MC(:,P_bars==p_bar));
    NSE_ES_prelim = std(ES_prelim(:,P_bars==p_bar));
    
    mean_VaR_IS = mean(VaR_IS(:,P_bars==p_bar));
    mean_ES_IS = mean(ES_IS(:,P_bars==p_bar));

    NSE_VaR_IS = std(VaR_IS(:,P_bars==p_bar));
    NSE_ES_IS = std(ES_IS(:,P_bars==p_bar));

    fprintf('(%s) 100*%4.2f%% VaR prelim (mean) estimate: %6.4f. \n', model, p_bar, mean_VaR_prelim);
    fprintf('(%s) NSE VaR prelim: %6.4f. \n', model, NSE_VaR_prelim);
    fprintf('(%s) VaR prelim: [%6.4f, %6.4f]. \n \n', model, mean_VaR_prelim - NSE_VaR_prelim, mean_VaR_prelim + NSE_VaR_prelim);

    fprintf('(%s) 100*%4.2f%% VaR IS (mean) estimate: %6.4f. \n', model, p_bar, mean_VaR_IS);
    fprintf('(%s) NSE VaR IS estimate: %6.4f. \n', model, NSE_VaR_IS);
    fprintf('(%s) VaR: [%6.4f, %6.4f]. \n \n', model, mean_VaR_IS - NSE_VaR_IS, mean_VaR_IS + NSE_VaR_IS);

    fprintf('(%s) 100*%4.2f%% ES prelim (mean) estimate: %6.4f. \n', model, p_bar, mean_ES_prelim);
    fprintf('(%s) NSE ES prelim: %6.4f. \n', model, NSE_ES_prelim);
    fprintf('(%s) ES prelim: [%6.4f, %6.4f]. \n \n', model, mean_ES_prelim - NSE_ES_prelim, mean_ES_prelim + NSE_ES_prelim);

    fprintf('(%s) 100*%4.2f%% ES IS (mean) estimate: %6.4f. \n', model, p_bar, mean_ES_IS);
    fprintf('(%s) NSE ES IS estimate: %6.4f. \n', model, NSE_ES_IS);
    fprintf('(%s) ES: [%6.4f, %6.4f]. \n', model, mean_ES_IS - NSE_ES_IS, mean_ES_IS + NSE_ES_IS);
    
    if plot_on2
        figure(390+100*p_bar)
        set(gcf, 'visible', 'off');
%         set(gcf,'defaulttextinterpreter','latex');
        boxplot([VaR_prelim_MC(:,P_bars==p_bar), VaR_IS(:,P_bars==p_bar)],'labels',{'VaR_prelim MC','VaR_IS'})        
        title(['(',model,' M = ',num2str(M),') ','100*', num2str(p_bar),'% VaR estimates: prelim and IS.'])
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        if print_on
            name = ['figures/(',model,')', num2str(p_bar),'_VaR_box_',num2str(M),'.png'];
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(3900+100*p_bar)
        set(gcf, 'visible', 'off');
%         set(gcf,'defaulttextinterpreter','latex');
        hold on; 
        bar(VaR_IS(:,P_bars==p_bar),'FaceColor',[0 0.4470 0.7410], 'EdgeColor','w'); 
        plot(0:(N_sim+1), (mean_VaR_prelim - NSE_VaR_prelim)*ones(N_sim+2,1),'r--'); 
        plot(0:(N_sim+1), (mean_VaR_prelim + NSE_VaR_prelim)*ones(N_sim+2,1),'r--'); 
        plot(0:(N_sim+1), mean_VaR_prelim*ones(N_sim+2,1),'r'); 
        hold off;
        title(['(',model,' M = ',num2str(M),') ','100*', num2str(p_bar),'% VaR IS estimates and the mean VaR prelim (+/- NSE VaR prelim).'])
    
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        if print_on
            name = ['figures/(',model,')', num2str(p_bar),'_VaR_bar_',num2str(M),'.png'];
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    end
    if save_on
        gen_out2;
    end
end

% set(1:n, 'visible', 'on');
set(390+100*P_bars, 'visible', 'on');
set(3900+100*P_bars, 'visible', 'on');

% model = 'arch';
% s = 'mitisem';
% hp = 1;
% gen_out;
% 
% clear GamMat x fig
% save(['results/arch_mitisem_',num2str(hp),'.mat']);
% 
% x = (0:0.00001:50)'+0.00001;
% GamMat = gamma(x);
