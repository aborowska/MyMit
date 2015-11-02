% All computational functions are in the folder "include"
% exept for MitISEM which is in the main "MyMit" folder.
% MitISEM is in 2 versions, the one with "_new" is currently being used.
% Plot functions are in include/Plots and are produced when plot_on = true.
% Prior and posterior functions are in incude/Models and they follow the
% idea from the R packages AdMit and MitISEM.

%% Initialisation
% clc
clear all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));
% addpath('include/Debug/');

plot_on = false; % whether the plots are genereated

% Tabularised values of the gamma function to speed up the copmutaions
x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);

% Artificial, white noise data 
T = 10000;
y = randn(T,1); 
y = y - mean(y);

% sigma is the VARIANCE of the error term, i.e. y_t ~ NID(0, sigma)
% should be sigma2, but I skip the power for readibility
sigma_init = 0.9;
M = 1000;
N_sim = 10; % number of MC replications for everything

p_bar = 0.050; % quantile alpha for VaR (100*alpha% VaR)
% [ususally 0.01 or 0.05; 0.5 or 0.99 for debugging]
% p_bar2 = 0.02; % quantile for "somewhat higher value of alpha" to compute
% VaR_prelim 

% Control parameters for  MitISEM 
MitISEM_Control
cont.mit.dfnc = 5;
cont2 = cont;

% hyper parameters for the prior for sigma2(inv. gamma)
a = 1; % if a == 0 then the flat prior is used
b = 1; 

% Posterior, theoretical moments
a_post = a + T/2;
b_post = b + sum(y.^2)/2;
mean_post = b_post/(a_post-1);
var_post = (b_post.^2)/(((a_post-1)^2)*(a_post-2));
std_post = sqrt(var_post);
        
% MC for VaR_prelim
VaR_prelim = zeros(N_sim,1);
ES_prelim = zeros(N_sim,1);
accept = zeros(N_sim,1);

VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);
hl_w = zeros(N_sim,1); % Sum of weights for high losses
hp_w = zeros(N_sim,1); % Sum of weights for high profits

for sim = 1:N_sim
    % logkernel
    kernel_init = @(x) - posterior_debug(x, y, a, b, true);
    kernel = @(x) posterior_debug(x, y, a, b, true);
 
    [mit1, summary1] = MitISEM_new(kernel_init, kernel, sigma_init, cont, GamMat);

    % draw from posterior
    [draw1, lnk1, ~] = fn_rmvgt_robust(M, mit1, kernel);
    % the moments of draw1 can be copmared with the theoretical moments 
    
    % Metropolis-Hastings to compute VaR_prelim
    [sigma1, accept(sim,1)] = Mit_MH(M+1000, kernel, mit1, GamMat);
    fprintf('(df = %i) MH acceptance rate: %4.2f. \n',cont.mit.dfnc, accept(sim,1));
    sigma1 = sigma1(1001:M+1000);

    % Future logreturns
    eps1 = randn(M,1);
    y_T1 = sqrt(sigma1).*eps1; 

    % fn_PL computes (among others) the profit-loss value of the given
    % returns
    [PL_T1, ind] = sort(fn_PL(y_T1));
    VaR_prelim(sim) = PL_T1(p_bar*M); % VaR_prelim = 0; VaR_prelim = Inf;
    ES_prelim(sim) = mean(PL_T1(1:p_bar*M));    
    fprintf('(df = %i) Preliminary 100*%4.2f%% VaR estimate: %6.4f. \n', cont.mit.dfnc, p_bar, VaR_prelim(sim,1));


    sigma1_hl = sigma1(ind); sigma1_hl = sigma1_hl(1:p_bar*M); 
    eps_hl = eps1(ind); eps_hl = eps_hl(1:p_bar*M);
    draw_hl = [sigma1_hl, eps_hl];

    % Take as mu_hl the draws "before" the draw which leads to VaR_prelim 
    %(to ensure that the restictions in the numerical optimaisation for the initial coponent are satisfied)
    mu_hl = draw_hl(end-1,:);
    PL_mu_hl = fn_PL(sqrt(mu_hl(1,1))*mu_hl(1,2)); % to see what the correspondinf PL value is



    %% High loss 
    % kernel_init = @(x) - posterior_debug_hl(x, y, a, b, VaR_constr, true);
    % kernel = @(x) posterior_debug_hl(x, y, a, b, VaR_constr, true);
    kernel_init = @(x) - posterior_debug_hl(x, y, a, b, VaR_prelim(sim,1), true);
    kernel = @(x) posterior_debug_hl(x, y, a, b, VaR_prelim(sim,1), true);

    % [mit2, summary2] = MitISEM(mit_hl, kernel, mu_hl, cont2, GamMat);
    [mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);

    % DRAWS FROM THE OPTIMAL CANDIDATE
    % M draws from the whole space and M draws from the high loss subspace
    % ==> 2M draws in total corresponding to the optimal 50-50 candidate
    
    % draw sigma2 from the posterior approximation (mit1) 
    % lnk1 is the correponding log kernel evaluation
    kernel = @(x) posterior_debug(x, y, a, b, true);
    [draw1, lnk1, ~] = fn_rmvgt_robust(M, mit1, kernel); % Robust drawing from the multivariate mixture of t
    
    % draw future disturbance for the standard normal distribution
    eps1 = randn(M,1);
    draw1 = [draw1,eps1]; % "Combined" draw (parameter + disturbance)
    
    % Draw from the high loss approximation (mit2)
    kernel = @(x) posterior_debug_hl(x, y, a, b, VaR_prelim(sim,1), true);
    [draw2, lnk2, ~] = fn_rmvgt_robust(M, mit2, kernel);
    
    % "Optimal" draw = from the optimal 50-50 candidate
    draw_opt = [draw1; draw2];
        
    % IMPORTANCE WEIGHTS
    % Use the computed logkernel evaluations
    % (for the draw from "the whole" we need to compute the corresponding logkernel value for epsilon)
    lnk_opt = lnk1 - 0.5*(log(2*pi) + eps1.^2);
    lnk_opt = [lnk_opt; lnk2]; % combine the logkernel evaluations

    % formula (10) from the QERMit paper: evaluation on the optimal candidate
    % 0.5*q_1(sigma2)*p(eps) [independent distubances]
%     exp_lnd1 = 0.5*normpdf(draw_opt(:,2)).*dmvgt(draw_opt(:,1),mit1,false, GamMat);
    exp_lnd1 = 0.5*exp(-0.5*(log(2*pi) + draw_opt(:,2).^2) + dmvgt(draw_opt(:,1), mit1, true, GamMat));

    % 0.5*q_2(sigma2,eps) 
    exp_lnd2 = 0.5*dmvgt(draw_opt,mit2,false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd); % take log to comute the importance weights in fn_ISwgts
    w_opt =  fn_ISwgts(lnk_opt, lnd_opt, false); % false - not normalised --> will be in fn_PL function
    
    if (p_bar == 0.5) % then the theoretical VaR is 0
        hl_w(sim,1) = sum( w_opt(draw_opt(:,2)<0,:)/sum(w_opt) );    
        hp_w(sim,1) = sum( w_opt(draw_opt(:,2)>0,:)/sum(w_opt) );    
    else
        hl_w(sim,1) = sum( w_opt(sqrt(draw_opt(:,1)).*draw_opt(:,2)<VaR_prelim(sim,1),:)/sum(w_opt) );    
        hp_w(sim,1) = sum( w_opt(sqrt(draw_opt(:,1)).*draw_opt(:,2)>VaR_prelim(sim,1),:)/sum(w_opt) );  
    end     
    fprintf('Sum of weights for high losses: %6.4f and for high profits: %6.4f.\n', hl_w(sim,1), hp_w(sim,1));
    
    % VaR_IS ESTIMATION
    y_opt = sqrt(draw_opt(:,1)).*draw_opt(:,2); % the correponding logreturns
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1); %  computed for the case when (L == 1) 
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);
    
    fprintf('(WN) IS 100*%4.2f%% VAR estimate: %6.4f. \n', p_bar, VaR_IS(sim,1));
    fprintf('(WN) IS 100*%4.2f%%ES estimate: %6.4f. \n', p_bar, ES_IS(sim,1));  
end

figure(2000)
hold on
xx = 0.45:0.01:0.55;
yy = 0.55:-0.01:0.45;
plot(xx,yy,'r')
scatter(hl_w,hp_w)
plot(xx,xx,'g')
hold off
ylabel({'sum of weights for high profits'});
xlabel({'sum of weights for high losses'});

mean_VaR_IS = mean(VaR_IS);
mean_ES_IS = mean(ES_IS);

NSE_VaR_IS = std(VaR_IS);
NSE_ES_IS = std(ES_IS);

fprintf('(WN) 100*%4.2f%% VaR prelim (mean) estimate: %6.4f. \n', p_bar, mean(VaR_prelim));
fprintf('(WN) NSE VaR prelim: %6.4f. \n',std((VaR_prelim)));
fprintf('(WN) VaR prelim: [%6.4f, %6.4f]. \n \n', mean(VaR_prelim)-std(VaR_prelim), mean(VaR_prelim)+std(VaR_prelim));

fprintf('(WN) 100*%4.2f%% VaR IS (mean) estimate: %6.4f. \n', p_bar, mean_VaR_IS);
fprintf('(WN) NSE VaR IS estimate: %6.4f. \n',NSE_VaR_IS);
fprintf('(WN) VaR: [%6.4f, %6.4f]. \n \n',mean_VaR_IS-NSE_VaR_IS,mean_VaR_IS+NSE_VaR_IS);

fprintf('(WN) 100*%4.2f%% ES prelim (mean) estimate: %6.4f. \n', p_bar, mean(ES_prelim));
fprintf('(WN) NSE ES prelim: %6.4f. \n',std((ES_prelim)));
fprintf('(WN) ES prelim: [%6.4f, %6.4f]. \n \n', mean(ES_prelim)-std(ES_prelim), mean(ES_prelim)+std(ES_prelim));

fprintf('(WN) 100*%4.2f%% ES IS (mean) estimate: %6.4f. \n', p_bar, mean_ES_IS);
fprintf('(WN) NSE ES IS estimate: %6.4f. \n',NSE_ES_IS);
fprintf('(WN) ES: [%6.4f, %6.4f]. \n',mean_ES_IS-NSE_ES_IS,mean_ES_IS+NSE_ES_IS);

figure(190)
boxplot([VaR_prelim, VaR_IS],'labels',{'VaR_prelim','VaR_IS'})

figure(200)
% hold on; bar(VaR_IS); plot(VaR_prelim*ones(N_sim,1),'r'); hold off;
scatter(VaR_prelim,VaR_IS)
ylabel({'VaR IS'});
xlabel({'VaR prelim'});
X = [ones(N_sim,1),VaR_prelim];
coef = inv(X'*X)*X'*(VaR_IS);
s2 = sum((VaR_IS-X*coef).^2)/(N_sim-2);
var_coef = s2*inv(X'*X);
tstat = coef./sqrt(diag(var_coef));
ttest = tinv(1-0.05/2,N_sim-2);
