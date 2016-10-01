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

algo = 'AdMit';
model = 'arch';

plot_on = false;
print_on  = false;
plot_on2 = true;
save_on = false;

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
N_sim = 100;

plot_on = false;
print_on  = false;
plot_on2 = true;
save_on = false;

AdMit_Control
cont.resampl_on = false;
% cont.IS.opt = true;
cont.IS.opt = false;
cont.IS.scale = [1, 0.25, 4]; 
cont.IS.perc = [0.25, 0.30, 0.35];
cont.CV_tol = 0.1;
cont.Hmax = 10;
cont.dfnc = 5;


cont2 = cont; 
cont2.Ns = 1e4;
cont2.Hmax = 10;
cont2.CV_tol = 0.01;
cont2.dfnc = 5;
cont2.IS.scale = [  0.35  ]; 
cont2.IS.perc = [0.3 0.25   ]; 


% p_bar = 0.05; % p_bar = 1-alpha, 100alpha% VaR
% P_bars = [0.01, 0.05, 0.1, 0.5];
P_bars = 0.01;

VaR_prelim = zeros(N_sim,1);
VaR_prelim_MC = zeros(N_sim,length(P_bars));
ES_prelim = zeros(N_sim,length(P_bars));
accept = zeros(N_sim,length(P_bars));

VaR_IS = zeros(N_sim,length(P_bars));
ES_IS = zeros(N_sim,length(P_bars));

for p_bar = P_bars
    fprintf('\np_bar: %4.2f\n',p_bar);
    %% QERMit 1a.: 
    kernel_init = @(a) - posterior_arch(a, data, S, true);
    kernel = @(a) posterior_arch(a, data, S, true);
    [mit1, summary1] = AdMit_old(kernel_init, kernel, mu_init, cont, GamMat);
    
    if save_on
        save(['results/arch_mit1_',algo,'.mat'],'mit1','summary1','cont','mu_init','p_bar');
    end
    
    for sim = 1:N_sim
          fprintf('\nPrelim sim = %i.\n', sim);
%         [mit1, summary1] = AdMit(kernel_init, kernel, mu_init, cont, GamMat);

        %% QERMit 1b.:
        % generate set of draws of alpha using independence MH with candidate from MitISEM; 
        % then simulate returns from normal with variance based on the draw of alpha 

        [alpha1, accept(sim,P_bars==p_bar)] = Mit_MH(M+1000, kernel, mit1, GamMat);
        fprintf('MH acceptance rate: %4.2f (%s, %s). \n', accept(sim,P_bars==p_bar), model, algo);
        alpha1 = alpha1(1001:M+1000);

        stdev = f_stdev(alpha1);
        eps1 = randn(M,1);
        y_T1 = stdev.*eps1;

        % get the preliminary VaR estimate as the 100th of the ascendingly sorted percentage loss values
        [PL_T1, ind] = sort(fn_PL(y_T1));
        VaR_prelim(sim,1) = PL_T1(p_bar*M); % VaR_prelim = 0; VaR_prelim = Inf;
        ES_prelim(sim,P_bars==p_bar) = mean(PL_T1(1:p_bar*M));    
        fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(sim,1), model, algo);
    end

    % take one value of VaR_prelim to construct mit2
    VaR_prelim_MC(:,P_bars==p_bar) = VaR_prelim;
%     VaR_prelim = VaR_prelim_MC(N_sim,1);       % the last one
    VaR_prelim = mean(VaR_prelim);              % the mean
 
    if save_on
        save(['results/arch_prelim_',algo,'.mat'],'mit1','accept','alpha1', 'eps1', 'y_T1', 'summary1',...
            'cont','p_bar','M','N_sim','VaR_prelim_MC','VaR_prelim','ES_prelim','ind');
    end
    
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
   
    [mit2, summary2] = AdMit_old(kernel_init, kernel, mu_hl, cont2, GamMat);
    
    if save_on
        save(['results/arch_mit2_',algo,'.mat'],'mit2','summary2','cont2','mu_hl','p_bar','VaR_prelim');
    end
    
    arch_plot2; %  The high loss density and the approximation to the high loss density
    
    %% QERMit 2: 
    % use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
    % to estimate VaR and ES for theta and y (or alpha in eps)

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MONTE CARLO VaR_IS and ES_IS (and their NSEs) ESTIMATION 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for sim = 1:N_sim
        resampl_on = false;
        fprintf('NSE sim = %i.\n', sim);    
        
%         kernel = @(a) posterior_arch(a, data, S, true);
%         [draw1, lnk1, ~] = fn_rmvgt_robust(M, mit1, kernel, resampl_on);
        draw1 = rmvgt2(M/2, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
        eps1 =  randn(M/2,1);
        draw1_eps1 = [draw1, eps1];

%         lnk1 = lnk1 - 0.5*(log(2*pi) + eps1.^2);
  
%         kernel = @(a) posterior_arch_hl(a, data, S, Inf, true);
%         [draw2, lnk2, ~] = fn_rmvgt_robust(M, mit2, kernel, resampl_on);
        
        draw2 = rmvgt2(M/2, mit2.mu, mit2.Sigma, mit2.df, mit2.p); 
        draw_opt = [draw1_eps1; draw2];
       
        arch_plot3; %  Draws from the optimal importance density, the joint density and the approximation to the optimal posterior

        %% IS weights
    %     kernel = @(a) posterior_arch_whole(a, data, S, true);
    %     lnk = kernel(draw_opt);
%         lnk = [lnk1; lnk2];
        
        kernel = @(a) posterior_arch(a, data, S, true);
        lnk = kernel(draw_opt(:,1));
        eps_pdf = normpdf(draw_opt(:,2));
        lnk = lnk + log(eps_pdf);
        
        
        exp_lnd1 = 0.5*normpdf(draw_opt(:,2)).*dmvgt(draw_opt(:,1), mit1, false, GamMat);
        exp_lnd2 = 0.5*dmvgt(draw_opt, mit2, false, GamMat);
        exp_lnd = exp_lnd1 + exp_lnd2;
        lnd = log(exp_lnd);
        w_opt = fn_ISwgts(lnk, lnd, false);

        %% VaR and ES IS estimates
        y_opt = f_stdev(draw_opt(:,1)).*draw_opt(:,2);
        dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
        IS_estim = fn_PL(dens, 1);
        VaR_IS(sim,P_bars==p_bar) = IS_estim(1,1);
        ES_IS(sim,P_bars==p_bar) = IS_estim(1,2);
        PL_opt = fn_PL(y_opt);

        arch_plot4; % The sorted future profit/losses 

        fprintf('IS 100*%4.2f%% VAR estimate: %6.4f (%s, %s). \n', p_bar, VaR_IS(sim,P_bars==p_bar), model, algo);
        fprintf('IS 100*%4.2f%% ES estimate: %6.4f (%s, %s). \n', p_bar, ES_IS(sim,P_bars==p_bar), model, algo);  
    end
    
    if save_on
        save(['results/arch_IS_',algo,'.mat'],'mit1','mit2','accept','draw_opt', 'y_opt', 'lnk','lnd','summary1','summary2',...
            'cont','cont2','p_bar','M','N_sim','VaR_prelim','VaR_IS','ES_IS');
    end
    
    mean_VaR_prelim = mean(VaR_prelim_MC(:,P_bars==p_bar));
    mean_ES_prelim = mean(ES_prelim(:,P_bars==p_bar));

    NSE_VaR_prelim = std(VaR_prelim_MC(:,P_bars==p_bar));
    NSE_ES_prelim = std(ES_prelim(:,P_bars==p_bar));
    
    mean_VaR_IS = mean(VaR_IS(:,P_bars==p_bar));
    mean_ES_IS = mean(ES_IS(:,P_bars==p_bar));

    NSE_VaR_IS = std(VaR_IS(:,P_bars==p_bar));
    NSE_ES_IS = std(ES_IS(:,P_bars==p_bar));


    fprintf('100*%4.2f%% VaR prelim (mean) estimate: %6.4f (%s, %s). \n', p_bar, mean_VaR_prelim, model, algo);
    fprintf('NSE VaR prelim: %6.4f (%s, %s). \n', NSE_VaR_prelim, model, algo);
    fprintf('VaR prelim: [%6.4f, %6.4f] (%s, %s). \n \n', mean_VaR_prelim - NSE_VaR_prelim, mean_VaR_prelim + NSE_VaR_prelim, model, algo);

    fprintf('100*%4.2f%% VaR IS (mean) estimate: %6.4f (%s, %s). \n',  p_bar, mean_VaR_IS, model, algo);
    fprintf('NSE VaR IS estimate: %6.4f (%s, %s). \n', NSE_VaR_IS, model, algo);
    fprintf('VaR: [%6.4f, %6.4f] (%s, %s). \n \n', mean_VaR_IS - NSE_VaR_IS, mean_VaR_IS + NSE_VaR_IS, model, algo);

    fprintf('100*%4.2f%% ES prelim (mean) estimate: %6.4f (%s, %s). \n', p_bar, mean_ES_prelim, model, algo);
    fprintf('NSE ES prelim: %6.4f (%s, %s). \n',NSE_ES_prelim, model, algo);
    fprintf('ES prelim: [%6.4f, %6.4f] (%s, %s). \n \n', mean_ES_prelim - NSE_ES_prelim, mean_ES_prelim + NSE_ES_prelim, model, algo);

    fprintf('100*%4.2f%% ES IS (mean) estimate: %6.4f (%s, %s). \n', p_bar, mean_ES_IS, model, algo);
    fprintf('NSE ES IS estimate: %6.4f (%s, %s). \n', NSE_ES_IS, model, algo);
    fprintf('ES: [%6.4f, %6.4f] (%s, %s). \n',  mean_ES_IS - NSE_ES_IS, mean_ES_IS + NSE_ES_IS, model, algo);
    
    if plot_on2
        figure(390+100*p_bar)
%         set(gcf, 'visible', 'off');
        if v_new
            set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
        else
            set(gcf,'units','normalized','outerposition',[0 0 0.5 0.75]);
        end
        set(gcf,'defaulttextinterpreter','latex');
        title(['100*', num2str(p_bar),'\% VaR estimates: prelim and IS (',strrep(model,'_','\_'),', ',algo,', M = ',num2str(M),', N\_sim = ', num2str(N_sim),').'])  
        title(['100*', num2str(p_bar),'% VaR estimates: prelim and IS (',model,', ',algo,', M = ',num2str(M),', N\_sim = ', num2str(N_sim),').'])
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        if print_on
            name = ['figures/(',model,'_',algo,')', num2str(p_bar),'_VaR_box_',num2str(M),'.png'];
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%
%         
%         figure(3900+100*p_bar)
% %         set(gcf, 'visible', 'off');
%         if v_new
%             set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
%         else
%             set(gcf,'units','normalized','outerposition',[0 0 0.5 0.75]);
%         end
%         set(gcf,'defaulttextinterpreter','latex');
%         hold on; 
%         bar(VaR_IS(:,P_bars==p_bar),'FaceColor',[0 0.4470 0.7410], 'EdgeColor','w'); 
%         plot(0:(N_sim+1), (mean_VaR_prelim - NSE_VaR_prelim)*ones(N_sim+2,1),'r--'); 
%         plot(0:(N_sim+1), (mean_VaR_prelim + NSE_VaR_prelim)*ones(N_sim+2,1),'r--'); 
%         plot(0:(N_sim+1), mean_VaR_prelim*ones(N_sim+2,1),'r'); 
%         hold off;
%         title(['100*', num2str(p_bar),'\% VaR IS estimates and the mean VaR prelim (+/- NSE VaR prelim) (',strrep(model,'_','\_'),', ',algo,', M = ',num2str(M),', N\_sim = ', num2str(N_sim),').'])    
%         if v_new
%             set(gca,'TickLabelInterpreter','latex')
%         else
%             plotTickLatex2D;
%         end
%         if print_on
%             name = ['figures/(',model,'_',algo,')', num2str(p_bar),'_VaR_bar_',num2str(M),'.png'];
%             fig = gcf;
%             fig.PaperPositionMode = 'auto';
%             print(name,'-dpng','-r0')
%         end
    end
    if save_on
        gen_out2;
    end
end
