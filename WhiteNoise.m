% All computational functions are in the folder "include"
% exept for MitISEM which is in the main "MyMit" folder.
% MitISEM is in 2 versions, the one with "_new" is currently being used.
% Plot functions are in include/Plots and are produced when plot_on = true.
% Prior and posterior functions are in incude/Models and they follow the
% idea from the R packages AdMit and MitISEM.

%% Initialisation
% clc
clear all
close all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));
% addpath('include/Debug/');


v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014b)')
    v_new = 1;
else
    v_new = 0;
end

plot_on = false; % whether the plots are genereated
print_on = false;
plot_on2 = true;
save_on = false;

% Tabularised values of the gamma function to speed up the copmutaions
x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);

model = 'WN';
% Artificial, white noise data 
T = 10000;
y = randn(T,1); 
y = y - mean(y);

% sigma is the VARIANCE of the error term, i.e. y_t ~ NID(0, sigma)
% should be sigma2, but I skip the power for readibility
sigma_init = 0.9;
M = 10000;
N_sim = 20; % number of MC replications

H = 10; % prediction horizon

% Control parameters for  MitISEM 
MitISEM_Control
cont.mit.dfnc = 5;
cont.mit.N = 10000;
cont.resmpl_on = false;

% Different control parameters might be used for high loss approximation
cont2 = cont;
cont2.mit.dfnc = 1;
cont2.mit.Hmax = 10;


% hyper parameters for the prior for sigma2(inv. gamma)
a = 1; % if a == 0, then the flat prior is used; if a == 1, then the conjugate prior (inv. gamma)
b = 1; 

% Posterior, theoretical moments
a_post = a + T/2;
b_post = b + sum(y.^2)/2;
mean_post = b_post/(a_post-1);
var_post = (b_post.^2)/(((a_post-1)^2)*(a_post-2));
std_post = sqrt(var_post);

% P_bars = [0.01, 0.05, 0.1, 0.5];
P_bars = [0.01];

% p_bar = 0.010; % quantile alpha for VaR (100*alpha% VaR)
% [ususally 0.01 or 0.05; 0.5 or 0.99 for debugging]

VaR_prelim = zeros(N_sim,1);
VaR_prelim_MC = zeros(N_sim,length(P_bars));
ES_prelim = zeros(N_sim,length(P_bars));
accept = zeros(N_sim,length(P_bars));
CV = zeros(N_sim,length(P_bars));

VaR_IS = zeros(N_sim,length(P_bars));
ES_IS = zeros(N_sim,length(P_bars));

for p_bar = P_bars
    fprintf('\np_bar: %4.2f\n',p_bar);

    % logkernel 
    kernel_init = @(x) - posterior_debug(x, y, a, b, true);
    kernel = @(x) posterior_debug(x, y, a, b, true);

    %% MC for VaR_prelim
    % If the MC compuations of VaR_prelim are supposed to be based on the same
    % mixture of t (so that we disregard the noise in the construction
    % process), then comment out the following line ...
    [mit1, summary1] = MitISEM_new(kernel_init, kernel, sigma_init, cont, GamMat);

    for sim = 1:N_sim
        
        % ... and comment the following line:
%         [mit1, summary1] = MitISEM_new(kernel_init, kernel, sigma_init, cont, GamMat);
%         CV(sim,P_bars==p_bar) = summary1.CV(end);
        
        % draw from posterior
        resampl_on = true; % ? true? false?
        [draw1, lnk1, ~] = fn_rmvgt_robust(M, mit1, kernel, resampl_on);
        % the moments of draw1 can be copmared with the theoretical moments 

        % Metropolis-Hastings to compute VaR_prelim
        [sigma1, accept(sim,P_bars==p_bar)] = Mit_MH(M+1000, kernel, mit1, GamMat);
        fprintf('(%s) MH acceptance rate: %4.2f. \n', model, accept(sim,P_bars==p_bar));
        sigma1 = sigma1(1001:M+1000);

        % Future logreturns
        eps1 = randn(M,H);
        y_T1 = repmat(sqrt(sigma1),1,H).*eps1; 

        % fn_PL computes (among others) the profit-loss value of the given returns
        [PL_T1, ind] = sort(fn_PL(y_T1));
        VaR_prelim(sim,1) = PL_T1(p_bar*M); % VaR_prelim = 0; VaR_prelim = Inf;
        ES_prelim(sim,P_bars==p_bar) = mean(PL_T1(1:p_bar*M));    
        fprintf('(%s) Preliminary 100*%4.2f%% VaR estimate: %6.4f. \n', model, p_bar, VaR_prelim(sim,1));
     end

    % take one value of VaR_prelim to construct mit2
    VaR_prelim_MC(:,P_bars==p_bar) = VaR_prelim;
%     VaR_prelim = VaR_prelim_MC(N_sim,1);       % the last one
    % VaR_prelim = 0;
    VaR_prelim = mean(VaR_prelim);              % the mean

    wn_plot0;  % Plot the approximation to the posterior

    %% Choose the starting point (mu_hl) for the constuction of the approximation 
    % to the high loss (hl) region density
    sigma1_hl = sigma1(ind); 
    eps_hl = eps1(ind,:); 
    draw_hl = [sigma1_hl, eps_hl];
    % Take as mu_hl the last draw which leads to the PL lower than VaR_prelim 
    % (to ensure that the restictions in the numerical optimisation for the initial coponent are satisfied)
    mu_hl = draw_hl(max(find(PL_T1 < VaR_prelim)),:);

    %% Might be based on weighted average of the "negative draws"
    % sigma1_hl = sigma1_hl(1:p_bar*M);
    % eps1_hl = eps1_hl(1:p_bar*M);
    % lnk_hl = kernel(sigma1_hl) - 0.5*(log(2*pi) + eps_hl.^2);
    % w_hl = lnk_hl - max(lnk_hl);
    % w_hl = exp(w_hl);
    % w_hl = w_hl/sum(w_hl);
    % [mu_hl, Sigma_hl] = fn_muSigma(draw_hl, w_hl);

    %% Alternatively: take the inital component as given (e.g. when numerical optimisation crashes)
    % (not based on mode and Hessian of the target)
    % use draws leading to very negative losses
    % mit_hl.mu = mu_hl;
    % mit_hl.Sigma = Sigma_hl;
    % mit_hl.df = 5;
    % mit_hl.p = 1;

    %% High loss  
    kernel_init = @(x) - posterior_debug_hl(x, y, a, b, VaR_prelim, true);
    kernel = @(x) posterior_debug_hl(x, y, a, b, VaR_prelim, true);

    % [mit2, summary2] = MitISEM(mit_hl, kernel, mu_hl, cont2, GamMat);
    [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);
  
    wn_plot1; %  Plots for high loss density and its approximation by mit2

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MONTE CARLO VaR_IS and ES_IS (and their NSEs) ESTIMATION 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for sim = 1:N_sim
        resampl_on = false;
        % DRAWS FROM THE OPTIMAL CANDIDATE
        % M draws from the whole space and M draws from the high loss subspace
        % ==> 2M draws in total corresponding to the optimal 50-50 candidate

        % draw sigma2 from the posterior approximation (mit1) 
        % lnk1 is the correponding log kernel evaluation
           
        kernel = @(x) posterior_debug(x, y, a, b, true);
        [draw1, lnk1, ~] = fn_rmvgt_robust(M, mit1, kernel, resampl_on); % Robust drawing from the multivariate mixture of t
        % draw future disturbance for the standard normal distribution
        eps1 = randn(M,H);
        draw1 = [draw1, eps1]; % "Combined" draw (parameter + disturbance)
        % the logkernel evaluations for the draw from "the whole": corrected by the corresponding logkernel value for epsilon)        
        lnk1 = lnk1 - 0.5*(H*log(2*pi) + sum(eps1.^2,2));

        % Draw from the high loss approximation (mit2)
        kernel = @(x) posterior_debug_hl(x, y, a, b, Inf, true);
        [draw2, lnk2, ~] = fn_rmvgt_robust(M, mit2, kernel, resampl_on);

        % "Optimal" draw = from the optimal 50-50 candidate
        draw_opt = [draw1; draw2];
        
        wn_plot2;  % Approximation to the optimal posterior density

        % IMPORTANCE WEIGHTS
        % Use the computed logkernel evaluations: combine the logkernel evaluations
        lnk_opt = [lnk1; lnk2];

        % formula (10) from the QERMit paper: evaluation on the optimal candidate
        
        % 0.5*q_1(sigma2)*p(eps) [independent distubances]

%        exp_lnd1 = 0.5*normpdf(draw_opt(:,2)).*dmvgt(draw_opt(:,1),mit1,false, GamMat);
        exp_lnd1 = 0.5*exp(-0.5*(H*log(2*pi) + sum(draw_opt(:,2:H+1).^2,2)) + dmvgt(draw_opt(:,1), mit1, true, GamMat));
        
        % 0.5*q_2(sigma2,eps) 
        exp_lnd2 = 0.5*dmvgt(draw_opt,mit2,false, GamMat);
        exp_lnd = exp_lnd1 + exp_lnd2;
        lnd_opt = log(exp_lnd); % take log to comute the importance weights in fn_ISwgts
        w_opt =  fn_ISwgts(lnk_opt, lnd_opt, false); % false - not normalised --> will be in fn_PL function
        
        % VaR_IS ESTIMATION
        y_opt = repmat(sqrt(draw_opt(:,1)),1,H).*draw_opt(:,2:H+1); % the correponding logreturns
        dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
        IS_estim = fn_PL(dens, 1); %  computed for the case when (L == 1) 
        VaR_IS(sim,P_bars==p_bar) = IS_estim(1,1);
        ES_IS(sim,P_bars==p_bar) = IS_estim(1,2);

        fprintf('(%s) IS 100*%4.2f%% VAR estimate: %6.4f. \n', model, p_bar, VaR_IS(sim,P_bars==p_bar));
        fprintf('(%s) IS 100*%4.2f%%ES estimate: %6.4f. \n', model, p_bar, ES_IS(sim,P_bars==p_bar));  

    %     PL_opt = fn_PL(y_opt);
    %     [PL_opt_h1, ind] = sort(PL_opt);
    %     hold on
    %     plot(PL_opt_h1)
    %     pos =  max(find(PL_opt_h1<=VaR_IS(sim,1)));
    %     scatter(pos, VaR_IS(sim,1),'MarkerFaceColor','red')
    %     hold off

    end
    
%     if plot_on2
%         figure(2000+100*p_bar)
%         hold on
%         xx = 0.45:0.01:0.55;
%         yy = 0.55:-0.01:0.45;
%         plot(xx,yy,'r')
%         scatter(hl_w,hp_w)
%         plot(xx,xx,'g')
%         hold off
%         ylabel({'sum of weights for high profits'});
%         xlabel({'sum of weights for high losses'});
%     end
    
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
        figure(180+100*p_bar)
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
        
        figure(1800+100*p_bar)
        set(gcf, 'visible', 'off');
%         set(gcf,'defaulttextinterpreter','latex');
        hold on; 
        bar(VaR_IS(:,P_bars==p_bar),'FaceColor',[0 0.4470 0.7410], 'EdgeColor','w'); 
        plot(0:(N_sim+1), (mean_VaR_prelim - NSE_VaR_prelim)*ones(N_sim+2,1),'r--'); 
        plot(0:(N_sim+1), (mean_VaR_prelim + NSE_VaR_prelim)*ones(N_sim+2,1),'r--'); 
        plot(0:(N_sim+1), mean_VaR_prelim*ones(N_sim+2,1),'r'); 
        hold off;
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        title(['(',model,' M = ',num2str(M),') ','100*', num2str(p_bar),'% VaR IS estimates and the mean VaR prelim (+/- NSE VaR prelim).'])

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
set(180+100*P_bars, 'visible', 'on');
set(1800+100*P_bars, 'visible', 'on');
